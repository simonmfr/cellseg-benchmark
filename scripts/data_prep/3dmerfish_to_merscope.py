#!/usr/bin/env python
"""Convert Fang et al. 2024 (eLife) 3D MERFISH transcripts to MERSCOPE-style detected_transcripts.csv.

barcode_id is the row of codebook.csv; global_x/y/z are already in µm (1 µm z planes).
Usage: 3dmerfish_to_merscope.py [out_dir]   (default: Fang_2024_eLife/merscope_format)
"""
import pathlib
import sys
import pandas as pd

FANG = (pathlib.Path(__file__).parents[2] / "data").resolve().parent / "merfish-public-data/Fang_2024_eLife"
SAMPLES = {"3DMERFISH_cortex": "mouse.cortex.242genes.100um", "3DMERFISH_hypothalamus": "mouse.hypothalamus.156genes.200um"}
out_dir = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else FANG / "merscope_format"

for key, raw in SAMPLES.items():
    genes = pd.read_csv(FANG / "data" / raw / "codebook.csv", usecols=["name"])["name"]
    (out_dir / key).mkdir(parents=True, exist_ok=True)
    reader = pd.read_csv(FANG / "data" / raw / "transcripts.csv", usecols=["barcode_id", "global_x", "global_y", "global_z"], chunksize=5_000_000)
    for i, t in enumerate(reader):
        t.assign(gene=genes.to_numpy()[t.barcode_id])[["global_x", "global_y", "global_z", "gene"]].to_csv(
            out_dir / key / "detected_transcripts.csv", mode="w" if i == 0 else "a", header=i == 0, index=False)
    print(key, "->", out_dir / key)
