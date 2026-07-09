#!/usr/bin/env python
import argparse
import os
import pathlib
import yaml

parser = argparse.ArgumentParser(
    description="Prepare scripts for vpt pipeline. Single staining."
)
parser.add_argument("staining", help="Name of staining (e.g. nuclei or PolyT).")
args = parser.parse_args()
BASE_PATH = pathlib.Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")
REPO_PATH = "$HOME/gitrepos/cellseg-benchmark"
SBATCH_DIR = BASE_PATH / "misc/sbatches/sbatch_vpt_2D_simple"

with open(BASE_PATH / "misc/sample_metadata.yaml") as f:
    data = yaml.safe_load(f)

EXPERIMENT_JSON_PATH = f"{REPO_PATH}/configs/vpt_2D_{args.staining}.json"
RUN_LOG_PATH = f"{REPO_PATH}/scripts/sbatch_utils/run_log.sh"

SBATCH_DIR.mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    res_path = f"{BASE_PATH}/samples/{key}/results/vpt_2D_DAPI_{args.staining}"
    for dire in os.listdir(value["path"]):
        if dire.endswith(".vzg"):
            vzg_path = pathlib.Path(value["path"], dire)
    f = open(SBATCH_DIR / f"{key}_{args.staining}.sbatch", "w")
    f.write(f"""#!/bin/bash
#SBATCH -p lrz-hgx-h100-94x4,lrz-hgx-a100-80x4,lrz-dgx-a100-80x8
#SBATCH -t 1-12:00:00
#SBATCH --mem=600G
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=40
#SBATCH -J vpt2D_{key}_{args.staining}
#SBATCH -o {BASE_PATH}/misc/logs/outputs/vpt2D_{key}_{args.staining}.out
#SBATCH -e {BASE_PATH}/misc/logs/errors/vpt2D_{key}_{args.staining}.err
#SBATCH --container-image="{BASE_PATH}/misc/enroot_images/vpt.sqsh"

set -euo pipefail

git -C {REPO_PATH} pull -q
source {RUN_LOG_PATH}

RES_PATH="{res_path}"
EXPERIMENT_JSON="{EXPERIMENT_JSON_PATH}"
INPUT_IMAGES="{pathlib.Path(value["path"], "images")}"
INPUT_TRANSFORM="{pathlib.Path(value["path"], "images/micron_to_mosaic_pixel_transform.csv")}"
INPUT_TRANSCRIPTS="{pathlib.Path(value["path"], "detected_transcripts.csv")}"
VZG_PATH="{vzg_path}"

BOUNDARIES="{pathlib.Path(res_path, "analysis_outputs/cellpose2_micron_space.parquet")}"
OUT_ANALYSIS="{pathlib.Path(res_path, "analysis_outputs")}"
OUT_CBG="{pathlib.Path(res_path, "analysis_outputs/cell_by_gene.csv")}"
OUT_META="{pathlib.Path(res_path, "analysis_outputs/cell_metadata.csv")}"
OUT_VZG="{pathlib.Path(res_path, "visualize.vzg")}"
TMP_PATH="{pathlib.Path(res_path, "tmp")}"

KEY="{key}"
STAINING="{args.staining}"
INPUT_PATH="${{INPUT_IMAGES}}"
RESULT_DIR="${{RES_PATH}}"
CMD="vpt run-segmentation + partition-transcripts + derive-entity-metadata + update-vzg"
start_run_log

mamba activate vpt

mkdir -p "${{RES_PATH}}"

vpt --verbose --processes 40 run-segmentation \\
  --segmentation-algorithm "${{EXPERIMENT_JSON}}" \\
  --input-images "${{INPUT_IMAGES}}" \\
  --input-micron-to-mosaic "${{INPUT_TRANSFORM}}" \\
  --output-path "${{OUT_ANALYSIS}}" \\
  --tile-size 2400 \\
  --tile-overlap 200

vpt --verbose partition-transcripts \\
  --input-boundaries "${{BOUNDARIES}}" \\
  --input-transcripts "${{INPUT_TRANSCRIPTS}}" \\
  --output-entity-by-gene "${{OUT_CBG}}"

vpt --verbose derive-entity-metadata \\
  --input-boundaries "${{BOUNDARIES}}" \\
  --output-metadata "${{OUT_META}}"

vpt --verbose --processes 10 update-vzg \\
  --input-vzg "${{VZG_PATH}}" \\
  --input-boundaries "${{BOUNDARIES}}" \\
  --input-entity-by-gene "${{OUT_CBG}}" \\
  --output-vzg "${{OUT_VZG}}" \\
  --input-metadata "${{OUT_META}}" \\
  --temp-path "${{TMP_PATH}}"

# remove leftover temp dir if vpt emptied it
rmdir --ignore-fail-on-non-empty "${{TMP_PATH}}" 2>/dev/null || true
""")
    f.close()

print(f"Wrote {len(data)} sbatch scripts to {SBATCH_DIR}")
