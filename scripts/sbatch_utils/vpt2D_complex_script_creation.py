import argparse
from os import listdir
from os.path import join
from pathlib import Path

import yaml

parser = argparse.ArgumentParser(description="Prepare scripts for vpt pipeline.")
parser.add_argument("staining1", help="Name of cell-boundary staining.")
parser.add_argument("staining2", help="Name of nucleus staining.")
args = parser.parse_args()

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sample_metadata.yaml"
) as f:
    data = yaml.save_load(f)

experiment_json_path = f"/dss/dsshome1/00/ra87rib/cellseg-benchmark/misc/vpt_experiment_jsons/{args.staining1}_{args.staining2}.json"

Path(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_vpt_2D_complex"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    res_path = f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/vpt_2D_DAPI_{args.staining1}_{args.staining2}"
    for dire in listdir(value["path"]):
        if dire.endswith(".vzg"):
            vzg_path = join(value["path"], dire)
    f = open(
        f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_vpt_2D_complex/{key}_{args.staining1}_{args.staining2}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash

#SBATCH -p lrz-hgx-a100-80x4
#SBATCH -t 12:00:00
#SBATCH --mem=600G
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=40
#SBATCH -J vtp2D_{key}_{args.staining1}_{args.staining2}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/vpt2D_{key}_{args.staining1}_{args.staining2}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/vpt2D_{key}_{args.staining1}_{args.staining2}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/vpt.sqsh"

mamba activate vpt
mkdir -p {res_path}

vpt --verbose --processes 40 run-segmentation \
--segmentation-algorithm {experiment_json_path} \
--input-images {join(value["path"], "images")} \
--input-micron-to-mosaic {join(value["path"], "images/micron_to_mosaic_pixel_transform.csv")} \
--output-path {join(res_path, "analysis_outputs")} \
--tile-size 2400 \
--tile-overlap 200

vpt --verbose partition-transcripts \
--input-boundaries {join(res_path, "analysis_outputs/cellpose2_micron_space.parquet")} \
--input-transcripts {join(value["path"], "detected_transcripts.csv")} \
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
