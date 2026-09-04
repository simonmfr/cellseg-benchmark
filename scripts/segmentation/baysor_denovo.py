#!/usr/bin/env python
import argparse
import pathlib
import subprocess
import sys

import pandas as pd

parser = argparse.ArgumentParser(
    description="Compute Baysor segmentation without prior, using the native Baysor CLI."
)
parser.add_argument("data_path", help="Path to merfish output folder.")
parser.add_argument("sample", help="sample name.")
parser.add_argument("dimension", choices=["2D", "3D"], help="segmentation mode.")
parser.add_argument("--z_step", type=float, default=1.5, help="micron spacing of z planes.")
parser.add_argument("--z_start", type=float, default=1.5, help="micron position of plane 0.")
parser.add_argument("--keep_transcripts", action="store_true")
args = parser.parse_args()

# CLI built into the micromamba env by baysor_denovo_build.sh (see that file for setup).
BAYSOR = "/dss/dsshome1/0C/ra98gaq/micromamba/envs/baysor_denovo/bin/baysor"
COLS = ["global_x", "global_y", "global_z", "gene"]


def write_transcripts(data_path, out_csv, z_start, z_step):
    """MERSCOPE detected_transcripts.csv -> Baysor input (x, y, z, gene)."""
    src = pathlib.Path(data_path, "detected_transcripts.csv")
    kept, header = 0, True
    with open(out_csv, "w") as f:
        for c in pd.read_csv(src, usecols=COLS, chunksize=2_000_000):
            c = c[~c.gene.astype(str).str.lower().str.startswith("blank")]
            if c.empty:
                continue
            pd.DataFrame({
                "x": c.global_x.to_numpy(),
                "y": c.global_y.to_numpy(),
                "z": z_start + c.global_z.to_numpy() * z_step,
                "gene": c.gene.to_numpy(),
            }).to_csv(f, index=False, header=header)
            header = False
            kept += len(c)
    if not kept:
        sys.exit(f"no molecules read from {src}")
    print(f"{kept:,} molecules -> {out_csv}", flush=True)


def main(data_path, sample, dimension, z_start, z_step, keep_transcripts):
    """Baysor without prior segmentation, native CLI, parquet output."""
    path = pathlib.Path(
        "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples",
        sample,
        "results",
        f"Baysor_{dimension}_denovo",
    )
    path.mkdir(parents=True, exist_ok=True)
    config = pathlib.Path(__file__).parents[2] / "configs" / "baysor_denovo.toml"
    transcripts = path / "transcripts.csv"

    if not transcripts.exists():
        write_transcripts(data_path, transcripts, z_start, z_step)

    cmd = [BAYSOR, "run", "-c", str(config), "--output-style", "parquet", "-o", str(path)]
    if dimension == "2D":
        cmd.append("--force-2d")
    cmd.append(str(transcripts))
    subprocess.run(cmd, check=True)

    if not keep_transcripts:
        transcripts.unlink()


if __name__ == "__main__":
    main(
        args.data_path,
        args.sample,
        args.dimension,
        args.z_start,
        args.z_step,
        args.keep_transcripts,
    )
