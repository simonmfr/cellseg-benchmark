import argparse
import json
from os import listdir
from os.path import join
from pathlib import Path

parser = argparse.ArgumentParser(
    description="Prepare scripts for vpt pipeline. Only cell-boundary staining."
)
parser.add_argument("staining", help="Name of cell-boundary staining.")
args = parser.parse_args()

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/sample_paths.json"
) as f:
    data = json.load(f)

experiment_json_path = f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/vpt_experiment_jsons/{staining}.json"

Path(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_vpt_2D_simple"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    res_path = f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/vpt_2D_DAPI_{args.staining}"
    for dire in listdir(value):
        if dire.endswith(".vzg"):
            vzg_path = join(value, dire)
    f = open(
        f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_vpt_2D_simple/{key}_{args.staining}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash

#SBATCH -p lrz-hgx-a100-80x4
#SBATCH -t 12:00:00
#SBATCH --mem=600G
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=40
#SBATCH -J vtp2D_{key}_{args.staining}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/vpt2D_{key}_{args.staining}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/vpt2D_{key}_{args.staining}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/vpt.sqsh"

mamba activate vpt
mkdir -p {res_path}

vpt --verbose --processes 40 run-segmentation \
--segmentation-algorithm {experiment_json_path} \
--input-images {join(value, "images")} \
--input-micron-to-mosaic {join(value, "images/micron_to_mosaic_pixel_transform.csv")} \
--output-path {join(res_path, "analysis_outputs")} \
--tile-size 2400 \
--tile-overlap 200

vpt --verbose partition-transcripts \
--input-boundaries {join(res_path, "analysis_outputs/cellpose2_micron_space.parquet")} \
--input-transcripts {join(value, "detected_transcripts.csv")} \
--output-entity-by-gene {join(res_path, "analysis_outputs/cell_by_gene.csv")}

vpt --verbose derive-entity-metadata \
--input-boundaries {join(res_path, "analysis_outputs/cellpose2_micron_space.parquet")} \
--output-metadata {join(res_path, "analysis_outputs/cell_metadata.csv")}

vpt --verbose --processes 10 update-vzg \
--input-vzg {vzg_path} \
--input-boundaries {join(res_path, "analysis_outputs/cellpose2_micron_space.parquet")} \
--input-entity-by-gene {join(res_path, "analysis_outputs/cell_by_gene.csv")} \
--output-vzg {join(res_path, "visualize.vzg")} \
--input-metadata {join(res_path, "analysis_outputs/cell_metadata.csv")} \
--temp-path {join(res_path, "tmp")}
            """)
    f.close()
