#!/usr/bin/env python
import argparse
import pathlib
import shutil

import pandas as pd
import sopa
from spatialdata import SpatialData, read_zarr
from spatialdata.models import PointsModel

parser = argparse.ArgumentParser(
    description="ProSeg 3D on transcripts only, seeded by Baysor 3D molecule assignments."
)
parser.add_argument(
    "baysor_dir", help="Baysor_3D_denovo results folder (with molecules.parquet)."
)
parser.add_argument("save_path", help="Output folder.")
parser.add_argument("--crop", nargs=3, type=float, metavar=("XMIN", "YMIN", "WIDTH"))
parser.add_argument(
    "proseg_flags", nargs=argparse.REMAINDER, help="Additional flags to pass to proseg."
)
args = parser.parse_args()


def main():
    """Proseg 3D seeded by Baysor 3D assignments via sopa."""
    save_path = pathlib.Path(args.save_path)
    save_path.mkdir(parents=True, exist_ok=True)

    mol = pd.read_parquet(
        pathlib.Path(args.baysor_dir, "molecules.parquet"),
        columns=["x", "y", "z", "gene", "cell", "is_noise"],
    )
    if args.crop:
        x0, y0, w = args.crop
        mol = mol[mol.x.between(x0, x0 + w) & mol.y.between(y0, y0 + w)]
    assigned = ~mol.is_noise & ~mol.cell.astype(str).isin(["0", "", "nan", "None"])
    mol["baysor_cell"] = pd.factorize(mol.cell.where(assigned))[0] + 1
    # proseg 3.0.7 fits z tilt on the first 1M rows only
    mol = mol.drop(columns=["cell", "is_noise"]).sample(frac=1, random_state=0)
    mol = mol.reset_index(drop=True)
    print(f"{len(mol)} molecules, {mol.baysor_cell.max()} prior cells")

    sdata = SpatialData(
        points={
            "transcripts": PointsModel.parse(mol, feature_key="gene", npartitions=64)
        }
    )
    sdata.write(save_path / "sdata_tmp.zarr", overwrite=True)
    sdata = read_zarr(save_path / "sdata_tmp.zarr")

    sopa.make_transcript_patches(
        sdata, patch_width=None, prior_shapes_key="baysor_cell", unassigned_value=0
    )
    sopa.segmentation.proseg(
        sdata, delete_cache=False, command_line_suffix=" ".join(args.proseg_flags)
    )
    sopa.aggregate(
        sdata, gene_column="gene", aggregate_channels=False, min_transcripts=10
    )

    cache_dir = sopa.utils.get_cache_dir(sdata)
    del sdata["transcripts"]
    sdata.write(save_path / "sdata.zarr", overwrite=True)
    shutil.copytree(cache_dir, save_path / "sdata.zarr" / cache_dir.name)
    shutil.rmtree(save_path / "sdata_tmp.zarr")


if __name__ == "__main__":
    main()
