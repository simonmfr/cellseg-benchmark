#!/usr/bin/env python
import argparse
import anndata as ad
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Convert a Ficture factor x gene posterior-count matrix to a "
        "SpatialData .zarr (table under .tables), or a standalone AnnData .h5ad."
    )
    parser.add_argument("input", help="model.posterior.count.tsv.gz")
    parser.add_argument("--h5ad", default=None, help="standalone AnnData .h5ad output")
    parser.add_argument("--zarr", default=None, help="SpatialData .zarr output")
    parser.add_argument("--table_name", default="ficture_factor")
    args = parser.parse_args()
    if not (args.h5ad or args.zarr):
        parser.error("provide --h5ad and/or --zarr")

    # rows = genes, cols = factor IDs; transpose so obs = factors, var = genes
    df = pd.read_csv(args.input, sep="\t", index_col=0).T
    adata = ad.AnnData(X=df.values.astype("float32"))
    adata.obs_names = [f"factor_{c}" for c in df.index]
    adata.var_names = df.columns

    if args.zarr:
        import spatialdata as sd

        sd.SpatialData(tables={args.table_name: adata}).write(args.zarr)
    if args.h5ad:
        adata.write_h5ad(args.h5ad)


if __name__ == "__main__":
    main()