#!/usr/bin/env python
"""Cut 3DMERFISH thick samples into 10 µm slabs and register them as samples.

Existing slab folders are skipped. Each slab keeps the transcripts with z0 <= global_z < z0 + 10 (global_z in µm, 1 µm planes),
shifted to start at 0, in a MERSCOPE-style folder next to the parent. Slabs are appended to
sample_metadata.yaml, so the standard Baysor 3D and Proseg 3D generators pick them up.
If the parent has images (mosaic_DAPI_z{global_z}.tif), each slab gets links to its planes, named by the shifted z.
"""
import pathlib
import pandas as pd
import yaml

meta_file = (pathlib.Path(__file__).parents[2] / "data").resolve() / "misc/sample_metadata.yaml"
meta = yaml.safe_load(open(meta_file))

for parent, z0s in {"3DMERFISH_cortex": range(0, 100, 10), "3DMERFISH_hypothalamus": range(0, 200, 10)}.items():
    src = pathlib.Path(meta[parent]["path"])
    out = {z0: d for z0 in z0s if not (d := src.parent / f"{parent}_z{z0}").exists()}
    for d in out.values():
        d.mkdir()
    for i, chunk in enumerate(pd.read_csv(src / "detected_transcripts.csv", chunksize=5_000_000) if out else []):
        for z0, d in out.items():
            chunk[(chunk.global_z >= z0) & (chunk.global_z < z0 + 10)].assign(global_z=lambda c: c.global_z - z0).to_csv(
                d / "detected_transcripts.csv", mode="w" if i == 0 else "a", header=i == 0, index=False)
    for z0 in z0s if (src / "images").exists() else []:
        img = src.parent / f"{parent}_z{z0}" / "images"
        img.mkdir(exist_ok=True)
        names = {"micron_to_mosaic_pixel_transform.csv": "micron_to_mosaic_pixel_transform.csv"} | {f"mosaic_DAPI_z{i}.tif": f"mosaic_DAPI_z{z0 + i}.tif" for i in range(10)}
        for name, parent_name in names.items():
            if (src / "images" / parent_name).exists() and not (img / name).exists():
                (img / name).symlink_to(src / "images" / parent_name)
    with open(meta_file, "a") as f:
        if not meta_file.read_text().endswith("\n"):
            f.write("\n")
        for d in out.values():
            if d.name not in meta:
                f.write(yaml.safe_dump({d.name: {**meta[parent], "path": str(d), "thickness": 10}}, sort_keys=False))
    print(parent, "->", [d.name for d in out.values()])
