#!/usr/bin/env python
"""One-off cleanup: remove outdated per-cell ovrlpy and FICTURE stats from master sdatas and adatas.

Drops the obsm entries ``Ovrlpy_stats`` and ``ficture_*`` in place, from every table in
``samples/*/sdata_z3.zarr`` and every ``analysis/*/*/adatas/*.h5ad*``. Z integrity is now
computed after merging (notebooks/metrics_z_integrity.ipynb). Dry run unless ``--apply``.
"""

import argparse
import pathlib
import sys

import h5py
import zarr

sys.path.insert(0, str(pathlib.Path(__file__).parents[2]))

from cellseg_benchmark import BASE_PATH


def outdated(keys):
    """Names of outdated obsm entries."""
    return [k for k in keys if k == "Ovrlpy_stats" or k.startswith("ficture")]


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--apply", action="store_true", help="Delete the entries (default: only list them).")
args = parser.parse_args()
base = pathlib.Path(BASE_PATH)

for sdata_path in sorted(base.glob("samples/*/sdata_z3.zarr")):
    changed = False
    for obsm_path in sorted(sdata_path.glob("tables/*/obsm")):
        group = zarr.open_group(str(obsm_path), mode="r+" if args.apply else "r")
        for key in outdated(group.keys()):
            print(obsm_path.relative_to(base), key)
            if args.apply:
                del group[key]
                changed = True
    if changed and (sdata_path / ".zmetadata").exists():
        zarr.consolidate_metadata(str(sdata_path))

for adata_path in sorted(base.glob("analysis/*/*/adatas/*.h5ad*")):
    with h5py.File(adata_path, "r+" if args.apply else "r") as f:
        for key in outdated(f.get("obsm", {}).keys()):
            print(adata_path.relative_to(base), key)
            if args.apply:
                del f["obsm"][key]
