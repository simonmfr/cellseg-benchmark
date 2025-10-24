import argparse
import logging
import re
from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from cellseg_benchmark.dea_utils import (
    add_ensembl_id,
    add_group_sample_counts,
    safe_sheet,
)

plt.rcParams["font.family"] = (
    "Arial" if "Arial" in [f.name for f in fm.fontManager.ttflist] else "sans-serif"
)
plt.rcParams["font.weight"] = "normal"


def get_args(test_args=None):  # noqa: D103
    p = argparse.ArgumentParser(description="Run single-cell-based wilxocon DE test")
    p.add_argument("cohort", help="Cohort name, e.g. foxf2")
    p.add_argument(
        "seg_method", help="Segmentation method, e.g. Cellpose_1_nuclei_model"
    )
    p.add_argument(
        "--sample_key", default="sample", help="obs column for donor/sample ID"
    )
    p.add_argument(
        "--subset_key",
        default="cell_type",
        help="obs column used to subset data (e.g. cell_type, cluster, region)",
    )
    p.add_argument(
        "--subset_values",
        nargs="+",
        default=None,
        help="Values of subset_key to process (default: all unique values)",
    )
    p.add_argument(
        "--condition_key",
        default="genotype",
        help="obs column for condition (e.g. genotype)",
    )
    p.add_argument("--ref", default="WT", help="Reference group (default: WT)")
    p.add_argument(
        "--test_groups",
        nargs="+",
        default=None,
        help="Groups to test vs reference (default: all groups in condition_key except --ref)",
    )
    p.add_argument(
        "--overwrite",
        type=lambda x: str(x).lower() in ["true", "1", "yes"],
        default=True,
        help="Overwrite existing result files (default: True)",
    )

    if test_args is not None:
        return p.parse_args(test_args)
    else:
        return p.parse_args()


if __name__ == "__main__":
    args = get_args()

    logger = logging.getLogger("dea")
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(
            logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s")
        )
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        logger.propagate = False

    base_path = Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")
    method_path = base_path / "analysis" / args.cohort / args.seg_method
    output_dir = method_path / "dea"
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Loading integrated AnnData...")
    adata = sc.read_h5ad(method_path / "adatas" / "adata_integrated.h5ad.gz")
    adata.obs["cell_type"] = adata.obs["cell_type_revised"]

    # re-group cell types
    adata.obs["cell_type"] = (
        adata.obs["cell_type"]
        .astype(str)
        .replace({
            "Tanycytes": "Ependymal",
            "Astroependymal": "Astrocytes",
            "Neurons-Glyc-Gaba": "Neurons-Other"
        })
        .astype("category")
    )

    if args.test_groups is None:
        groups = adata.obs[args.condition_key].dropna().unique().tolist()
        args.test_groups = [g for g in groups if g != args.ref]
        msg = (
            f"Test groups inferred from '{args.condition_key}': {', '.join(map(str, args.test_groups))}"
            if args.test_groups
            else f"No test groups found different from ref '{args.ref}' under '{args.condition_key}'."
        )
        logger.info(msg)

    # sample overview
    if args.subset_values is None:
        groups_to_process = sorted(adata.obs[args.subset_key].unique())
    else:
        groups_to_process = args.subset_values
    for cond, df in adata.obs.groupby(args.condition_key, observed=True):
        samples = df[args.sample_key].unique()
        logger.info(f"{cond}: {len(samples)} samples → {', '.join(samples)}")

    adata.obs["n_cells_sum"] = adata.obs.groupby(
        [args.subset_key, args.sample_key], observed=True
    )[args.sample_key].transform("count")
    adata.obs["volume_mean"] = adata.obs.groupby(
        [args.subset_key, args.sample_key], observed=True
    )["volume_final"].transform("mean")
    adata.obs["volume_sum"] = adata.obs.groupby(
        [args.subset_key, args.sample_key], observed=True
    )["volume_final"].transform("sum")

    adata.X = adata.layers["volume_log1p_norm"].copy()

    del adata.varm
    del adata.obsm
    del adata.obsp

    logger.info("Run DEA...")

    adata_dict = {}
    unique_groups = adata.obs[args.subset_key].unique().tolist()
    for key in unique_groups:
        tmp = adata[adata.obs[args.subset_key] == key].copy()
        adata_dict[key] = tmp
    adata_dict["all"] = adata.copy()

    results = []
    for subset, adata_tmp in adata_dict.items():
        sc.tl.rank_genes_groups(
            adata_tmp,
            groupby=args.condition_key,
            method="wilcoxon",
            reference=args.ref,
            corr_method="benjamini-hochberg",
        )

        valid_groups = set(
            adata_tmp.obs[args.condition_key].astype("category").cat.categories
        )
        take = [g for g in args.test_groups if g in valid_groups]

        dfs = []
        for g in take:
            df = sc.get.rank_genes_groups_df(adata_tmp, group=g).sort_values("pvals")
            df = df.copy()
            df["test_group"] = g
            dfs.append(df)

        if not dfs:
            continue

        df = pd.concat(dfs, ignore_index=True).rename(
            columns={
                "names": "gene",
                "logfoldchanges": "log2FC",
                "pvals": "PValue",
                "pvals_adj": "FDR",
            }
        )

        df["FC"] = np.power(2, df["log2FC"])
        df["method"] = "wilcoxon_scanpy"
        df["subset"] = subset
        df["ref"] = args.ref

        df = df[
            [
                "subset",
                "gene",
                "FC",
                "log2FC",
                "PValue",
                "FDR",
                "scores",
                "method",
                "test_group",
                "ref",
            ]
        ]
        results.append(df)

    df_all = pd.concat(results, ignore_index=True).sort_values(
        ["subset", "test_group", "PValue"]
    )

    logger.info("Format output table...")
    df_all["test"] = df_all["test_group"].astype(str) + "vs" + df_all["ref"].astype(str)
    df_all = add_group_sample_counts(
        df_all,
        adata_dict,
        condition_key=args.condition_key,
        sample_key=args.sample_key,
        ref=args.ref,
        test_groups=args.test_groups,
        subset_group="subset",
    )

    logger.info("Add ensembl gene ids...")
    df_all = add_ensembl_id(df_all, logger=logger)

    logger.info("Export table(s)...")
    subset_key_clean = re.sub(r"[_-]", "", args.subset_key)
    name = f"{args.cohort}-by-{subset_key_clean}"

    for (method, test), df_mt in df_all.groupby(["method", "test"]):
        df_mt = df_mt.loc[:, df_mt.notna().any()]  # drop all-NaN cols
        xlsx = output_dir / f"{name}_{method}_{test}.xlsx"
        if xlsx.exists() and not args.overwrite:
            logger.info(f"  Exists, skip: {xlsx}")
            continue
        used = set()
        with pd.ExcelWriter(xlsx, engine="xlsxwriter") as writer:
            for gid, g in df_mt.groupby("subset", sort=True):
                g.sort_values("PValue").to_excel(
                    writer, sheet_name=safe_sheet(gid, used), index=False
                )
        logger.info(f"  Wrote: {xlsx}")

    logger.info("Done.")
