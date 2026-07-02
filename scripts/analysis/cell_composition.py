#!/usr/bin/env python3
"""
This script is a version of test_comp_ana.ipynb that only computes cell cmposition analysis:

1. Reads integrated AnnData files for one or more segmentation methods.
2. Maps cells to anatomical regions using a region annotation parquet file.
3. Counts cells per sample, region, method, and cell type.
4. Calculates relative cell counts:
      cell_counts_rel = cell_count_reg_ct / cell_count_sum
   where cell_count_sum is the total number of cells in the selected regions
   for a given sample and method.
5. Saves cell_counts_rel line plots and summary tables.

Example:
    python cell_counts_rel_cli.py \
        --region-annot-path /path/to/aging_brain_regions.parquet \
        --adata-int-path /path/to/analysis/aging \
        --out-dir /path/to/output \
        --regions CTX \
        --cell-types cell type column from merged_all \
        --filter-undefined
"""

from __future__ import annotations

import argparse
import os
import warnings
from pathlib import Path
from typing import Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns
from geopandas import read_parquet
from joblib import Parallel, delayed
from tqdm import tqdm

from ma_package.spatial_mapping import map_points_to_regions_from_anndata

try:
    from cellseg_benchmark._constants import method_colors
except Exception:
    method_colors = {}


DEFAULT_REGIONS = ["BS","CTX","STR", "HIP", "fiber_tracts"]

DEFAULT_EXCLUDED_METHODS = ["SIS_DAPI_total_mrna", "Watershed_Merlin"]

DEFAULT_UNDEFINED_LABELS = [
    "undefined",
    "Undefined",
    "UNDEFINED",
    "unassigned",
    "Unassigned",
    "UNASSIGNED",
    "unknown",
    "Unknown",
    "UNKNOWN",
    "nan",
    "NaN",
    "NA",
    "None",
    "",
]


def configure_plotting() -> None:
    """Set plotting defaults used for saved figures."""
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
            "font.size": 10,
            "axes.labelsize": 11,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 8,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )
    sns.set_style("whitegrid")


def build_regions_by_slide(region_annot: pd.DataFrame) -> dict:
    """Convert the annotation GeoDataFrame into the mapping expected by ma_package."""
    regions_by_slide = {}
    for sample, gdf_sample in region_annot.groupby("sample"):
        inner = {}
        for region, gdf_region in gdf_sample.groupby("label"):
            inner[region] = list(gdf_region.geometry)
        regions_by_slide[sample] = inner
    return regions_by_slide


def process_method(
    method: str,
    adata_int_path: Path,
    regions_by_slide: dict,
    cell_type_col: str,
    volume_col: str,
    coord_key: str,
    slide_key: str,
    filter_undefined: bool,
    undefined_labels: set[str],
) -> pd.DataFrame:
    """Read one method's AnnData file, map cells to regions, and count cell types."""
    adata_path = adata_int_path / method / "adatas" / "adata_integrated.h5ad.gz"
    if not adata_path.exists():
        raise FileNotFoundError(f"Could not find AnnData file for method '{method}': {adata_path}")

    adata = sc.read_h5ad(adata_path)

    required_obs = {cell_type_col, slide_key, "age_months", volume_col}
    missing_obs = sorted(required_obs.difference(adata.obs.columns))
    if missing_obs:
        raise KeyError(
            f"Method '{method}' is missing required .obs columns: {missing_obs}. "
            f"Available columns include: {list(adata.obs.columns)[:30]}"
        )

    if coord_key not in adata.obsm:
        raise KeyError(
            f"Method '{method}' is missing .obsm['{coord_key}']. "
            f"Available .obsm keys: {list(adata.obsm.keys())}"
        )

    obs = adata.obs.copy()

    if filter_undefined:
        labels = obs[cell_type_col].astype("string")
        keep = labels.notna() & ~labels.isin(undefined_labels)
        obs = obs.loc[keep].copy()
        adata = adata[obs.index].copy()

    mapped = map_points_to_regions_from_anndata(
        adata,
        regions_by_slide=regions_by_slide,
        coord_key=coord_key,
        slide_key=slide_key,
        include_boundary=True,
        dissolve=False,
        index_kind="name",
        return_df=True,
    )

    df_all = pd.concat([v["df"] for v in mapped.values()], ignore_index=True)
    df_all = df_all.set_index("obs_id").reindex(adata.obs_names)

    obs["anatom_region"] = df_all["label"].values
    obs["region_poly_index"] = df_all["poly_index"].values

    # Keep only mapped cells.
    obs = obs.loc[obs["anatom_region"].notna()].copy()

    counts_samp_reg_ct = (
        obs.groupby([slide_key, "anatom_region", cell_type_col], observed=True)
        .agg(cell_count_reg_ct=("anatom_region", "size"), volume_sum=(volume_col, "sum"))
        .reset_index()
        .rename(columns={slide_key: "sample", cell_type_col: "cell_type_revised"})
    )

    counts_samp_reg = (
        obs.groupby([slide_key, "anatom_region"], observed=True)
        .agg(cell_count_reg=("anatom_region", "size"))
        .reset_index()
        .rename(columns={slide_key: "sample"})
    )

    age = obs[[slide_key, "age_months"]].drop_duplicates().rename(columns={slide_key: "sample"})

    merged = counts_samp_reg.merge(
        counts_samp_reg_ct, on=["sample", "anatom_region"]
    ).merge(age, on="sample")

    merged["proportion"] = merged["cell_count_reg_ct"] / merged["cell_count_reg"]
    merged["method"] = method

    return merged


def add_area_columns(merged_all: pd.DataFrame, region_annot: pd.DataFrame) -> pd.DataFrame:
    """Add region area and density/coverage columns to the count table."""
    region_annot = region_annot.copy()
    region_annot["area"] = region_annot["geometry"].area
    region_annot["area_mm2"] = region_annot["geometry"].area * 1e-6

    area_lookup = (
        region_annot.groupby(["sample", "label"], observed=True)
        .agg(area_sum=("area", "sum"), area_mm2_sum=("area_mm2", "sum"))
        .reset_index()
        .rename(columns={"label": "anatom_region"})
    )

    merged_all = merged_all.merge(area_lookup, on=["sample", "anatom_region"], how="left")

    merged_all["cell_counts_per_mm2"] = (
        merged_all["cell_count_reg_ct"] / merged_all["area_mm2_sum"]
    )
    merged_all["area_covered"] = merged_all["volume_sum"] * 1e-6 / merged_all["area_mm2_sum"]
    merged_all["method"] = merged_all["method"].astype("category")
    merged_all["cell_type_revised"] = merged_all["cell_type_revised"].astype("category")
    return merged_all


def subset_regions(
    df: pd.DataFrame,
    regions: Iterable[str] | None = None,
    celltypes: Iterable[str] | str | None = None,
    methods: Iterable[str] | None = None,
) -> pd.DataFrame:
    """Subset regions/cell types/methods and calculate relative cell counts."""
    if regions is None:
        regions = list(df["anatom_region"].dropna().unique())
    elif isinstance(regions, str):
        regions = [regions]

    if celltypes is None:
        celltypes = list(df["cell_type_revised"].dropna().unique())
    elif isinstance(celltypes, str):
        celltypes = [celltypes]

    if methods is None:
        methods = list(df["method"].dropna().unique())
    elif isinstance(methods, str):
        methods = [methods]

    mask = (
        df["anatom_region"].isin(regions)
        & df["cell_type_revised"].isin(celltypes)
        & df["method"].isin(methods)
    )

    df_tmp = df.loc[mask].copy()

    if df_tmp.empty:
        return df_tmp

    if hasattr(df_tmp["method"], "cat"):
        df_tmp["method"] = df_tmp["method"].cat.remove_unused_categories()
    if hasattr(df_tmp["cell_type_revised"], "cat"):
        df_tmp["cell_type_revised"] = df_tmp["cell_type_revised"].cat.remove_unused_categories()

    area_by_sample_method = (
        df_tmp.groupby(["sample", "anatom_region", "method"], observed=True)
        .agg(area_mm2_sum=("area_mm2_sum", "first"))
        .reset_index()
        .groupby(["sample", "method"], observed=True)
        .agg(area_mm2_sum=("area_mm2_sum", "sum"))
        .reset_index()
    )

    total_count_by_sample_method = (
        df_tmp.groupby(["sample", "anatom_region", "method"], observed=True)
        .agg(cell_count_sum=("cell_count_reg", "first"))
        .reset_index()
        .groupby(["sample", "method"], observed=True)
        .agg(cell_count_sum=("cell_count_sum", "sum"))
        .reset_index()
    )

    result = (
        df_tmp.groupby(["sample", "cell_type_revised", "method"], observed=True)
        .agg(cell_count_reg_ct=("cell_count_reg_ct", "sum"), volume_sum=("volume_sum", "sum"))
        .reset_index()
    )

    sample_age = df_tmp[["sample", "age_months"]].drop_duplicates()

    result = (
        result.merge(area_by_sample_method, on=["sample", "method"])
        .merge(total_count_by_sample_method, on=["sample", "method"])
        .merge(sample_age, on="sample")
    )

    result["cell_counts_per_mm2"] = result["cell_count_reg_ct"] / result["area_mm2_sum"]
    result["cell_counts_rel"] = result["cell_count_reg_ct"] / result["cell_count_sum"]
    return result


def save_cell_counts_rel_outputs(
    merged_all: pd.DataFrame,
    out_dir: Path,
    regions: list[str],
    cell_types: list[str],
    methods: list[str],
    formats: list[str],
) -> pd.DataFrame:
    """Save one cell_counts_rel plot per cell type and return combined plotting table."""
    all_plot_data = []

    for cell_type in cell_types:
        data = subset_regions(merged_all, regions=regions, celltypes=cell_type, methods=methods)
        if data.empty:
            warnings.warn(f"No data found for cell type '{cell_type}'. Skipping plot.")
            continue

        all_plot_data.append(data)
        safe_cell_type = str(cell_type).replace("/", "_").replace(" ", "_")

        data.to_csv(out_dir / f"cell_counts_rel_{safe_cell_type}_plot_data.csv", index=False)
        data.to_excel(out_dir / f"cell_counts_rel_{safe_cell_type}_plot_data.xlsx", index=False)

        g = sns.relplot(
            data=data,
            x="age_months",
            y="cell_composition",
            hue="method",
            kind="line",
            estimator="mean",
            palette=method_colors if method_colors else None,
            errorbar=None,
            height=6,
            aspect=1.6,
        )
        g.set_axis_labels("Age months", "Relative cell count")
        g.figure.suptitle(f"{cell_type}: relative cell counts", y=1.02)

        for fmt in formats:
            g.savefig(out_dir / f"cell_counts_rel_{safe_cell_type}.{fmt}", bbox_inches="tight")
        plt.close(g.figure)

    if not all_plot_data:
        raise ValueError("No plot data were generated. Check regions, cell types, and methods.")

    combined = pd.concat(all_plot_data, ignore_index=True)
    combined.to_csv(out_dir / "cell_counts_rel_plot_data_all_cell_types.csv", index=False)
    combined.to_excel(out_dir / "cell_counts_rel_plot_data_all_cell_types.xlsx", index=False)

    summary = (
        combined.groupby(["cell_type_revised", "method", "age_months"], observed=True)
        .agg(
            n_samples=("sample", "nunique"),
            mean_cell_counts_rel=("cell_counts_rel", "mean"),
            sd_cell_counts_rel=("cell_counts_rel", "std"),
            mean_cell_count=("cell_count_reg_ct", "mean"),
            mean_total_count_in_selected_regions=("cell_count_sum", "mean"),
            mean_cell_counts_per_mm2=("cell_counts_per_mm2", "mean"),
        )
        .reset_index()
    )
    summary.to_csv(out_dir / "cell_counts_rel_summary_table.csv", index=False)
    summary.to_excel(out_dir / "cell_counts_rel_summary_table.xlsx", index=False)

    return combined


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate cell_counts_rel plots and summary tables from cellseg-benchmark AnnData outputs."
    )
    parser.add_argument(
        "--region-annot-path",
        required=True,
        type=Path,
        help="Path to the anatomical region annotation parquet file, e.g. aging_brain_regions.parquet.",
    )
    parser.add_argument(
        "--adata-int-path",
        required=True,
        type=Path,
        help="Path to the directory containing one subfolder per segmentation method.",
    )
    parser.add_argument(
        "--out-dir",
        required=True,
        type=Path,
        help="Output folder for plots and summary tables.",
    )
    parser.add_argument(
        "--methods",
        nargs="*",
        default=None,
        help="Segmentation methods to include. Default: all method folders in --adata-int-path, after exclusions.",
    )
    parser.add_argument(
        "--exclude-methods",
        nargs="*",
        default=DEFAULT_EXCLUDED_METHODS,
        help="Method folder names to exclude when --methods is not provided.",
    )
    parser.add_argument(
        "--exclude-method-prefixes",
        nargs="*",
        default=["Negative"],
        help="Method name prefixes to exclude when --methods is not provided.",
    )
    parser.add_argument(
        "--regions",
        nargs="+",
        default=DEFAULT_REGIONS,
        help="Anatomical regions to combine for the cell_counts_rel plot.",
    )
    parser.add_argument(
        "--cell-type-col",
        default="cell_type_revised",
        help="Column in adata.obs containing cell type labels.",
    )
    parser.add_argument(
        "--volume-col",
        default="volume_final",
        help="Column in adata.obs containing cell/cell-segment volume values.",
    )
    parser.add_argument(
        "--coord-key",
        default="spatial_microns",
        help="Key in adata.obsm containing spatial coordinates.",
    )
    parser.add_argument(
        "--slide-key",
        default="sample",
        help="Column in adata.obs identifying slide/sample.",
    )
    parser.add_argument(
        "--filter-undefined",
        action="store_true",
        help="Filter out undefined/unassigned cell types before counting.",
    )
    parser.add_argument(
        "--undefined-labels",
        nargs="*",
        default=DEFAULT_UNDEFINED_LABELS,
        help="Cell type labels to remove when --filter-undefined is used.",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=1,
        help="Number of parallel jobs. Use -1 for all available cores.",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        default=["png", "pdf", "svg"],
        help="Image formats to save.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    configure_plotting()
    warnings.filterwarnings("ignore", category=FutureWarning)

    args.out_dir.mkdir(parents=True, exist_ok=True)

    region_annot = read_parquet(args.region_annot_path)
    regions_by_slide = build_regions_by_slide(region_annot)

    if args.methods is None:
        methods = [p.name for p in args.adata_int_path.iterdir() if p.is_dir()]
        methods = [
            m
            for m in methods
            if m not in set(args.exclude_methods)
            and not any(m.startswith(prefix) for prefix in args.exclude_method_prefixes)
        ]
        methods.sort()
    else:
        methods = args.methods

    if not methods:
        raise ValueError("No methods selected. Check --adata-int-path, --methods, and exclusions.")

    print(f"Selected {len(methods)} methods:")
    for method in methods:
        print(f"  - {method}")

    undefined_labels = set(args.undefined_labels)

    res_tables = Parallel(n_jobs=args.n_jobs)(
        delayed(process_method)(
            method=method,
            adata_int_path=args.adata_int_path,
            regions_by_slide=regions_by_slide,
            cell_type_col=args.cell_type_col,
            volume_col=args.volume_col,
            coord_key=args.coord_key,
            slide_key=args.slide_key,
            filter_undefined=args.filter_undefined,
            undefined_labels=undefined_labels,
        )
        for method in tqdm(methods, desc="Processing methods")
    )

    merged_all = pd.concat(res_tables, ignore_index=True)
    merged_all = add_area_columns(merged_all, region_annot)

    merged_all.to_csv(args.out_dir / "cell_counts_full_table.csv", index=False)
    merged_all.to_excel(args.out_dir / "cell_counts_full_table.xlsx", index=False)

    cell_types = (
        merged_all[args.cell_type_col]
        .dropna()
        .astype(str)
        .sort_values()
        .unique()
        .tolist()
    )

    save_cell_counts_rel_outputs(
        merged_all=merged_all,
        out_dir=args.out_dir,
        regions=args.regions,
        cell_types=cell_types,
        methods=methods,
        formats=args.formats,
    )

    print(f"Done. Outputs saved to: {args.out_dir}")


if __name__ == "__main__":
    main()
