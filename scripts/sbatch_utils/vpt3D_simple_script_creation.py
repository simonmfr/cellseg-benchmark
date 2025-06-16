import json
import sys
from os import listdir
from os.path import join
from pathlib import Path

staining = sys.argv[1]

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/sample_paths.json"
) as f:
    data = json.load(f)

experiment_json_path = f"/dss/dsshome1/00/ra87rib/cellseg-benchmark/misc/vpt_experiment_jsons/{staining}_3D.json"

Path(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_vpt_3D_simple"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    res_path = f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/vpt_3D_DAPI_{staining}"
    vzg_path = None
    for dire in listdir(value):
        if dire.endswith(".vzg") or dire.endswith(".vzg2"):
            vzg_path = join(value, dire)
    assert vzg_path is not None, f"{key} provides not valid vzg file"
    f = open(
        f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_vpt_3D_simple/{key}_{staining}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash

#SBATCH -p lrz-hgx-a100-80x4
#SBATCH -t 1-12:00:00
#SBATCH --mem=600G
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=40
#SBATCH -J vtp3D_{key}_{staining}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/vpt3D_{key}_{staining}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/vpt3D_{key}_{staining}.err
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
