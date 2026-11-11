#!/usr/bin/env python
"""Cut 3DMERFISH thick samples into 10 µm slabs and register them as samples.

Each slab keeps the transcripts with z0 <= global_z < z0 + 10 (global_z in µm, 1 µm planes),
shifted to start at 0, in a MERSCOPE-style folder next to the parent. Slabs are appended to
sample_metadata.yaml, so the standard Baysor 3D and Proseg 3D generators pick them up.
"""
import pathlib

import pandas as pd
import yaml

BASE_PATH = (pathlib.Path(__file__).parents[2] / "data").resolve()
THICK = 10
SLABS = {"3DMERFISH_cortex": [30, 50, 70], "3DMERFISH_hypothalamus": [50, 100, 150]}

meta_file = BASE_PATH / "misc/sample_metadata.yaml"
meta = yaml.safe_load(open(meta_file))

for parent, z0s in SLABS.items():
    src = pathlib.Path(meta[parent]["path"])
    keys = {z0: f"{parent}_z{z0}" for z0 in z0s}
    out = {z0: src.parent / k for z0, k in keys.items()}
    for d in out.values():
        d.mkdir(exist_ok=True)
    for i, chunk in enumerate(pd.read_csv(src / "detected_transcripts.csv", chunksize=5_000_000)):
        for z0, d in out.items():
            s = chunk[(chunk.global_z >= z0) & (chunk.global_z < z0 + THICK)].assign(global_z=lambda c: c.global_z - z0)
            s.to_csv(d / "detected_transcripts.csv", mode="w" if i == 0 else "a", header=i == 0, index=False)
    with open(meta_file, "a") as f:
        for z0, k in keys.items():
            if k not in meta:
                f.write(yaml.safe_dump({k: {**meta[parent], "path": str(out[z0]), "thickness": THICK}}, sort_keys=False))
    print(parent, "->", list(keys.values()))
