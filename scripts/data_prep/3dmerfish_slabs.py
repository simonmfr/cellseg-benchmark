#!/usr/bin/env python
"""Cut 3DMERFISH thick samples into 10 µm slabs and register them as samples.

Each slab keeps the transcripts with z0 <= global_z < z0 + 10 (global_z in µm, 1 µm planes),
shifted to start at 0, in a MERSCOPE-style folder next to the parent. Slabs are appended to
sample_metadata.yaml, so the standard Baysor 3D and Proseg 3D generators pick them up.
"""
import pathlib
import pandas as pd
import yaml

meta_file = (pathlib.Path(__file__).parents[2] / "data").resolve() / "misc/sample_metadata.yaml"
meta = yaml.safe_load(open(meta_file))

for parent, z0s in {"3DMERFISH_cortex": [30, 50, 70], "3DMERFISH_hypothalamus": [50, 100, 150]}.items():
    src = pathlib.Path(meta[parent]["path"])
    out = {z0: src.parent / f"{parent}_z{z0}" for z0 in z0s}
    for d in out.values():
        d.mkdir(exist_ok=True)
    for i, chunk in enumerate(pd.read_csv(src / "detected_transcripts.csv", chunksize=5_000_000)):
        for z0, d in out.items():
            chunk[(chunk.global_z >= z0) & (chunk.global_z < z0 + 10)].assign(global_z=lambda c: c.global_z - z0).to_csv(
                d / "detected_transcripts.csv", mode="w" if i == 0 else "a", header=i == 0, index=False)
    with open(meta_file, "a") as f:
        if not meta_file.read_text().endswith("\n"):
            f.write("\n")
        for d in out.values():
            if d.name not in meta:
                f.write(yaml.safe_dump({d.name: {**meta[parent], "path": str(d), "thickness": 10}}, sort_keys=False))
    print(parent, "->", [d.name for d in out.values()])
