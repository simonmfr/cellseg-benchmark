import logging
import random
import re
from typing import Any, Optional, Sequence

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.sparse import issparse


def pseudobulk_aggregate_and_filter(
    adata: sc.AnnData,
    subset_value: Optional[str],
    sample_key: str = "sample",
    subset_key: str = "cell_type",
    min_cells: int = 15,
    obs_to_keep: Optional[Sequence[str]] = None,
    replicates_per_patient: int = 1,
    logger: Optional[logging.Logger] = None,
) -> sc.AnnData:
    """Aggregates gene expression from individual cells into donor-level pseudobulk profiles while filtering out donors with insufficient cell counts.

    If subset_value is None, processes all cells without filtering (by subset_key).

    By setting replicates_per_patient > 1, n pseudo-replicates are created per donor
    by randomly splitting that donor's cells into n roughly equal parts.

    Based on https://www.sc-best-practices.org/conditions/differential_gene_expression.html.
    """
    if obs_to_keep is None:
        obs_to_keep = []
    if logger is None:
        logger = logging.getLogger(__name__)

    # Subset
    if subset_value is not None:
        adata_cell_pop = adata[adata.obs[subset_key] == subset_value].copy()
        if adata_cell_pop.n_obs == 0:
            logger.warning(
                f"\tNo cells for '{subset_key} == {subset_value}'. Skipping."
            )
            return sc.AnnData(
                X=np.zeros((0, adata.n_vars)),
                var=adata.var.copy(),
                obs=pd.DataFrame(index=[]),
            )
    else:
        logger.info(f"Processing all data (no subsetting by {subset_key})")
        adata_cell_pop = adata.copy()

    if sample_key not in adata_cell_pop.obs.columns:
        raise KeyError(f"'{sample_key}' not found in adata.obs columns.")

    # Per-donor sizes & filtering
    size_by_donor = adata_cell_pop.obs.groupby(sample_key, observed=True).size()
    donors_kept = [d for d, n in size_by_donor.items() if n >= min_cells]
    donors_dropped = {d: int(n) for d, n in size_by_donor.items() if n < min_cells}

    if donors_dropped:
        logger.info(f"\tDropping samples with <{min_cells} cells: {donors_dropped}")

    if len(donors_kept) == 0:
        logger.warning(
            f"\tAll donors below min_cells={min_cells} in this group; nothing to aggregate."
        )
        return sc.AnnData(
            X=np.zeros((0, adata_cell_pop.n_vars)),
            var=adata_cell_pop.var.copy(),
            obs=pd.DataFrame(index=[]),
        )

    extra_obs = [c for c in obs_to_keep if c != sample_key]
    out_rows = []

    adata_cell_pop.obs[sample_key] = adata_cell_pop.obs[sample_key].astype(str)

    total = len(donors_kept)
    for j, donor in enumerate(donors_kept, 1):
        if j in (1, total):
            logger.info(f"\tProcessing donor {j}/{total}...")
        adata_donor = adata_cell_pop[adata_cell_pop.obs[sample_key] == donor]
        if adata_donor.n_obs == 0:
            logger.debug(f"\tDonor {donor} has 0 cells after subsetting; skipping.")
            continue

        idx = list(adata_donor.obs_names)
        random.shuffle(idx)
        splits = np.array_split(np.array(idx), replicates_per_patient)

        for rep_idx, rep_indices in enumerate(splits):
            if len(rep_indices) == 0:
                continue

            adata_rep = adata_donor[rep_indices]
            X = (
                adata_rep.X.toarray()
                if issparse(adata_rep.X)
                else np.asarray(adata_rep.X)
            )

            df_rep = pd.DataFrame(
                X, index=adata_rep.obs_names, columns=adata_rep.var_names
            )
            obs_cols = [sample_key] + extra_obs
            df_rep = df_rep.join(adata_rep.obs[obs_cols])

            agg_dict = {gene: "sum" for gene in adata_rep.var_names}
            for obs in extra_obs:
                agg_dict[obs] = "first"

            gb = df_rep.groupby(sample_key, observed=True).agg(agg_dict)
            row = gb.loc[donor].copy()
            row[sample_key] = donor
            row.name = f"donor_{donor}_{rep_idx}"
            out_rows.append(row)

    if not out_rows:
        logger.warning("\tNo rows produced after aggregation (nothing to write).")
        return sc.AnnData(
            X=np.zeros((0, adata_cell_pop.n_vars)),
            var=adata_cell_pop.var.copy(),
            obs=pd.DataFrame(index=[]),
        )

    df = pd.DataFrame(out_rows)
    gene_cols = list(adata_cell_pop.var_names)
    meta_cols = [c for c in df.columns if c not in gene_cols]
    df = df[gene_cols + meta_cols]

    adata_out = sc.AnnData(
        X=df[gene_cols].values,
        var=pd.DataFrame(index=gene_cols),
        obs=df[meta_cols].copy(),
    )
    return adata_out


def prepare_adata_for_rpy2(adata_, key="cell_type"):
    """Make AnnData rpy2 compatible."""
    if isinstance(adata_.X, np.ndarray):
        adata_.X = sp.csr_matrix(adata_.X.astype(np.float32))
    else:
        adata_.X = sp.csr_matrix(adata_.X)

    for key_name in [
        "condition_colors",
        key + "_colors",
        key + "_sample_colors",
        "sample_colors",
        "batch_colors",
        "pca",
    ]:
        if key_name in adata_.uns and isinstance(adata_.uns[key_name], np.ndarray):
            adata_.uns[key_name] = adata_.uns[key_name].tolist()

    for key_name, value in adata_.obsm.items():
        if isinstance(value, np.ndarray):
            adata_.obsm[key_name] = pd.DataFrame(value, index=adata_.obs.index)

    if "X_pca" in adata_.obsm:
        adata_.obsm["X_pca"] = sp.csc_matrix(adata_.obsm["X_pca"])

    for key_name, value in adata_.varm.items():
        if isinstance(value, np.ndarray):
            adata_.varm[key_name] = pd.DataFrame(value, index=adata_.var.index)

    for col in adata_.obs.filter(like="volume").columns:
        adata_.obs[col] = adata_.obs[col].astype(float)

    if hasattr(adata_, "layers"):
        del adata_.layers

    return adata_


def add_ensembl_id(
    df: pd.DataFrame,
    species: str = "mouse",
    mg: Optional[Any] = None,
    logger: Optional[logging.Logger] = None,
    *,
    out_col: str = "ensembl_id",
    max_renames_preview: int = 20,
) -> pd.DataFrame:
    """Map gene symbols to Ensembl IDs, works if genes are in index or 'gene' column."""
    import mygene

    mg = mg or mygene.MyGeneInfo()

    # detect source
    source_is_index = "gene" not in df.columns
    genes = (df.index if source_is_index else df["gene"]).astype(str)

    sym_to_ens, alias_renames = {}, {}
    direct_q, alias_q = set(), set()

    for scope in ("symbol", "alias"):
        need = [
            g for g in genes.unique() if g not in sym_to_ens and g not in alias_renames
        ]
        if not need:
            break
        for h in mg.querymany(
            need, scopes=scope, species=species, fields="ensembl.gene,symbol"
        ):
            ens = _first_ens(h.get("ensembl"))
            canon, q = h.get("symbol") or h.get("query"), h.get("query")
            if not (ens and canon):
                continue
            sym_to_ens.setdefault(canon, ens)
            if scope == "symbol" and q:
                direct_q.add(q)
            elif scope == "alias" and q and q != canon:
                alias_renames[q] = canon
                alias_q.add(q)

    # apply alias renames
    if alias_renames:
        if source_is_index:
            df.rename(index=alias_renames, inplace=True)
        else:
            df["gene"] = df["gene"].replace(alias_renames)
        genes = (df.index if source_is_index else df["gene"]).astype(str)

    df[out_col] = genes.map(sym_to_ens).fillna("NA")

    if logger:
        n_input, normalized = (
            len(genes.unique()),
            {alias_renames.get(s, s) for s in genes.unique()},
        )
        mapped = sum(s in sym_to_ens for s in normalized)
        logger.info(
            f"{n_input} unique symbols -> {mapped} mapped ({mapped / n_input * 100:.1f}%)"
        )
        if alias_renames and max_renames_preview:
            items = list(alias_renames.items())
            preview = ", ".join(f"{a}->{s}" for a, s in items[:max_renames_preview])
            more = len(items) - min(len(items), max_renames_preview)
            logger.info(
                f"Renamed aliases: {preview}{' (+' + str(more) + ' more)' if more else ''}"
            )

    return df


def _first_ens(x):
    """Helper to extract first Ensembl gene ID from dict or list of dicts."""
    if isinstance(x, dict):
        return x.get("gene")
    if isinstance(x, list):
        for i in x:
            if "gene" in i:
                return i["gene"]


def add_group_sample_counts(
    df, adatas_pb, *, condition_key, sample_key, ref, test_groups, subset_group="subset"
):
    """Add per-group sample counts to a DEA results DataFrame."""
    sg = subset_group
    counts = pd.concat(
        [
            (
                a.obs.groupby(condition_key, observed=True)[sample_key]
                .nunique()
                .rename("n_samples")
                .reset_index()
                .assign(**{sg: g})
            )
            for g, a in adatas_pb.items()
        ],
        ignore_index=True,
    )

    wide = (
        counts.pivot(index=sg, columns=condition_key, values="n_samples")
        .fillna(0)
        .astype(int)
    )

    out = df.merge(
        counts[counts[condition_key].isin(test_groups)].rename(
            columns={condition_key: "test_group", "n_samples": "test_group_n"}
        )[[sg, "test_group", "test_group_n"]],
        on=[sg, "test_group"],
        how="left",
    ).merge(
        wide[[ref]].rename(columns={ref: "ref_n"}),
        left_on=sg,
        right_index=True,
        how="left",
    )

    out[["test_group_n", "ref_n"]] = (
        out[["test_group_n", "ref_n"]].fillna(0).astype(int)
    )
    # out = out.set_index("gene"); out.index.name = None
    return out


def safe_sheet(s, used):
    """Helper to return unique, Excel-safe sheet name."""
    s = re.sub(r"[][:*?/\\]+", "_", str(s))[:31] or "Sheet"
    base, i = s, 2
    while s in used:
        s = (base[:28] + f"_{i}")[:31]
        i += 1
    used.add(s)
    return s
