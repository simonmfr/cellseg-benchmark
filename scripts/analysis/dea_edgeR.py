import argparse
import logging
import re
import warnings
from importlib.resources import files
from pathlib import Path

import anndata as ad
import anndata2ri
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
import scanpy as sc
from rpy2.rinterface_lib.embedded import RRuntimeError
from rpy2.robjects.conversion import localconverter

import cellseg_benchmark as csb
from cellseg_benchmark._constants import cell_type_colors
from cellseg_benchmark.adata_utils import plot_pseudobulk_pca
from cellseg_benchmark.dea_utils import (
    add_ensembl_id,
    add_group_sample_counts,
    prepare_adata_for_rpy2,
    pseudobulk_aggregate_and_filter,
    safe_sheet,
)

plt.rcParams["font.family"] = (
    "Arial" if "Arial" in [f.name for f in fm.fontManager.ttflist] else "sans-serif"
)
plt.rcParams["font.weight"] = "normal"

warnings.filterwarnings("ignore", message=".*Observation names are not unique*")
VALID_METHODS = {"LRT", "QL"}


def get_args(test_args=None):  # noqa: D103
    p = argparse.ArgumentParser(description="Run edgeR pseudobulk-based DEA")
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
    p.add_argument(
        "--batch_key",
        default="slide",
        help="optional batch key for inclusion as covariate (default: slide)",
    )
    p.add_argument("--ref", default="WT", help="Reference group (default: WT)")
    p.add_argument(
        "--test_groups",
        nargs="+",
        default=None,
        help="Groups to test vs reference (default: all groups in condition_key except --ref)",
    )
    p.add_argument(
        "--edger_methods",
        nargs="+",
        default=["LRT"],
        choices=sorted(VALID_METHODS),
        help="edgeR method(s) to run (default: LRT)",
    )
    p.add_argument("--min_cells", type=int, default=15, help="Minimum cells per donor")
    p.add_argument(
        "--replicates_per_patient",
        type=int,
        default=1,
        help="Number of pseudoreplicates per donor",
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

    # Logger setup
    logger = logging.getLogger("dea")
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(
            logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s")
        )
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        logger.propagate = False
    rcb.logger.handlers = logger.handlers
    rcb.consolewrite_print = lambda x: logger.debug(f"R: {x.strip()}")
    rcb.consolewrite_error = lambda x: (_ for _ in ()).throw(RRuntimeError(x.strip()))
    rcb.consolewrite_message = lambda x: logger.info(f"R: {x.strip()}")
    rcb.consolewrite_warn = lambda x: (
        logger.warning if x.lstrip().lower().startswith("warning") else logger.info
    )(f"R: {x.strip()}")
    setattr(rcb, "consolewrite_warnerror", rcb.consolewrite_warn)

    conv = ro.default_converter + ro.pandas2ri.converter + anndata2ri.converter
    # r_script = Path(__file__).resolve().parent / "cellseg_benchmark" / "dea_utils.r"
    # ro.r["source"](str(r_script))
    # edgeR_loop = ro.globalenv["edgeR_loop"]
    ro.r["source"](str(files(csb) / "dea_utils.r"))
    edgeR_loop = ro.globalenv["edgeR_loop"]

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
        .replace(
            {
                "Tanycytes": "Ependymal",
                "Astroependymal": "Astrocytes",
                "Neurons-Glyc-Gaba": "Neurons-Other",
            }
        )
        .astype("category")
    )

    # Clean up group names for R conversion
    adata.obs[args.subset_key] = [
        key.replace(" ", "_")
        .replace("-", "_")
        .replace("/", "_")
        .replace("+", "")
        .replace("(", "")
        .replace(")", "")
        if not pd.isna(key)
        else key
        for key in adata.obs[args.subset_key]
    ]
    for col in [args.condition_key, args.sample_key, args.subset_key]:
        adata.obs[col] = adata.obs[col].astype("category")
    if args.test_groups is None:
        groups = adata.obs[args.condition_key].dropna().unique().tolist()
        args.test_groups = [g for g in groups if g != args.ref]
        msg = (
            f"Test groups inferred from '{args.condition_key}': {', '.join(map(str, args.test_groups))}"
            if args.test_groups
            else f"No test groups found different from ref '{args.ref}' under '{args.condition_key}'."
        )
        logger.info(msg)

    # print sample overview
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

    adata.X = adata.layers["counts"].copy()
    assert np.issubdtype(adata.X.dtype, np.integer)

    obs_to_keep = [
        args.condition_key,
        args.subset_key,
        args.sample_key,
        "n_cells_sum",
        "volume_mean",
        "volume_sum",
    ]
    if args.batch_key and args.batch_key in adata.obs.columns:
        obs_to_keep.append(args.batch_key)
    else:
        logger.warning(
            f"Batch column '{args.batch_key}' was not found in adata.obs "
            "or is invalid. Continuing **without** a batch covariate."
        )
        args.batch_key = None

    logger.info("Run pseudobulking...")
    total = len(groups_to_process)
    logger.info(f"Processing {groups_to_process[0]} ({1}/{total})")
    adata_pb = pseudobulk_aggregate_and_filter(
        adata,
        subset_value=groups_to_process[0],
        sample_key=args.sample_key,
        subset_key=args.subset_key,
        obs_to_keep=obs_to_keep,
        min_cells=args.min_cells,
        replicates_per_patient=args.replicates_per_patient,
        logger=logger,
    )
    for i, group in enumerate(groups_to_process[1:], 2):
        logger.info(f"Processing {group} ({i}/{total})")
        adata_pb_i = pseudobulk_aggregate_and_filter(
            adata,
            subset_value=group,
            sample_key=args.sample_key,
            subset_key=args.subset_key,
            obs_to_keep=obs_to_keep,
            min_cells=args.min_cells,
            replicates_per_patient=args.replicates_per_patient,
            logger=logger,
        )
        adata_pb = ad.concat(
            [adata_pb, adata_pb_i], join="outer", label=None, index_unique=None
        )
    adata_pb.obs_names_make_unique()

    # Add whole brain: pseudobulk per sample across all groups
    adata_pb_i = pseudobulk_aggregate_and_filter(
        adata,
        subset_value=None,
        sample_key=args.sample_key,
        subset_key=args.subset_key,
        obs_to_keep=obs_to_keep,
        min_cells=args.min_cells,
        replicates_per_patient=args.replicates_per_patient,
        logger=logger,
    )
    adata_pb_i.obs[args.subset_key] = "all"
    adata_pb = ad.concat([adata_pb, adata_pb_i])
    del adata

    plot_pseudobulk_pca(adata_pb, args, output_dir, cell_type_colors, logger)

    logger.info("Run DEA...")
    adata_pb = prepare_adata_for_rpy2(adata_pb, key=args.subset_key)
    adatas_pb = {}
    unique_groups = adata_pb.obs[args.subset_key].unique().tolist()
    for key in unique_groups:
        tmp = adata_pb[adata_pb.obs[args.subset_key] == key].copy()
        adatas_pb[key] = tmp
    all_degs = {}
    for group_i in adatas_pb:
        adata_tmp = adatas_pb[group_i]
        with localconverter(conv):
            combined_results = edgeR_loop(
                adata=adata_tmp,
                group_i=group_i,
                edger_methods=args.edger_methods,
                test_groups=args.test_groups,
                ref_group=args.ref,
                condition_col=args.condition_key,
                batch_col=args.batch_key,
            )
        if combined_results is not ro.NULL:
            all_degs[group_i] = combined_results

    logger.info("Format output table...")
    collapsed_df = pd.concat(
        [pd.DataFrame(v).assign(subset=k) for k, v in all_degs.items()],
        ignore_index=True,
    ).drop(columns=["result_id"], errors="ignore")

    collapsed_df = add_group_sample_counts(
        collapsed_df,
        adatas_pb,
        condition_key=args.condition_key,
        sample_key=args.sample_key,
        ref=args.ref,
        test_groups=args.test_groups,
        subset_group="subset",
    )
    collapsed_df = collapsed_df.set_index("gene")
    collapsed_df.insert(collapsed_df.columns.get_loc("FC"), "gene", collapsed_df.index)

    order = [
        "subset",
        "gene",
        "FC",
        "logFC",
        "PValue",
        "FDR",
        "logCPM",
        "LR",
        "method",
        "test_group",
        "ref",
        "test",
    ]
    collapsed_df = collapsed_df[
        order + [c for c in collapsed_df.columns if c not in order]
    ]

    logger.info("Add ensembl gene ids...")
    collapsed_df = add_ensembl_id(collapsed_df, logger=logger)

    logger.info("Export table(s)...")
    subset_key_clean = re.sub(r"[_-]", "", args.subset_key)
    name = f"{args.cohort}-by-{subset_key_clean}_pseudobulk_edgeR"

    for (method, test), df_mt in collapsed_df.groupby(["method", "test"]):
        df_mt = df_mt.loc[:, df_mt.notna().any()]  # drop all-NaN cols
        xlsx = output_dir / f"{name}_{method.split('_', 1)[-1]}_{test}.xlsx"
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
