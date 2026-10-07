#!/usr/bin/env python
import argparse
import pathlib

parser = argparse.ArgumentParser(description="Scripts for vascular subtyping.")
parser.add_argument("cohort", help="Cohort name, e.g., 'foxf2'")

args = parser.parse_args()
BASE_PATH = (pathlib.Path(__file__).parents[2] / "data").resolve()
SBATCH_DIR = BASE_PATH / "misc/sbatches/sbatch_vascular_subtyping"
condition_col = "genotype" if args.cohort == "foxf2" else "condition"
methods = [
    "Baysor_2D_Cellpose_1_DAPI_PolyT_0.2",
    "Baysor_2D_Cellpose_1_DAPI_PolyT_0.8",
    "Baysor_2D_Cellpose_1_DAPI_Transcripts_0.2",
    "Baysor_2D_Cellpose_1_DAPI_Transcripts_0.8",
    "Baysor_2D_Cellpose_1_nuclei_model_1.0",
    "Baysor_2D_Cellpose_2_DAPI_PolyT_0.2",
    "Baysor_2D_Cellpose_2_DAPI_PolyT_0.8",
    "Baysor_2D_Cellpose_2_DAPI_Transcripts_0.2",
    "Baysor_2D_Cellpose_2_DAPI_Transcripts_0.8",
    "Cellpose_1_DAPI_PolyT",
    "Cellpose_1_DAPI_Transcripts",
    "Cellpose_1_nuclei_model",
    "Cellpose_1_Merlin",
    "Cellpose_2_DAPI_PolyT",
    "Cellpose_2_DAPI_Transcripts",
    "Ficture_segments",
    "Ficture_segments_dapi",
    "Negative_Control_Rastered_5",
    "Negative_Control_Rastered_10",
    "Negative_Control_Rastered_25",
    "Negative_Control_Visium",
    "Negative_Control_Voronoi",
    "Proseg_3D_Cellpose_1_DAPI_Transcripts",
    "Proseg_3D_Cellpose_1_DAPI_PolyT",
    "Proseg_3D_Cellpose_1_nuclei_model",
    "Proseg_3D_Cellpose_2_DAPI_PolyT",
    "Proseg_3D_Cellpose_2_DAPI_Transcripts",
    "Proseg_3D_vpt3D_DAPI_nuclei",
    "Proseg_3D_vpt3D_DAPI_PolyT",
    "Proseg_3D_vpt3D_DAPI_PolyT_nuclei",
    "SIS_DAPI_total_mrna",
    "vpt_2D_DAPI_PolyT",
    "vpt_2D_DAPI_nuclei",
    "vpt_2D_DAPI_PolyT_nuclei",
    "vpt_3D_DAPI_PolyT",
    "vpt_3D_DAPI_nuclei",
    "vpt_3D_DAPI_PolyT_nuclei",
    "Watershed_Merlin"
]

SBATCH_DIR.mkdir(parents=False, exist_ok=True)

for seg_method in methods:
    if seg_method == "Negative_Control_Rastered_5":
        time_limit = "10:00:00"
    elif any(
        keyword in seg_method for keyword in ["Baysor", "Cellpose"]
    ) or seg_method in [
        "Negative_Control_Rastered_10",
        "Negative_Control_Voronoi",
    ]:
        time_limit = "05:00:00"
    else:
        time_limit = "03:00:00"

    memory = "270G" if "Negative_Control" in seg_method else "65G"

    job_name = f"vasc_subty_{args.cohort}_{seg_method}"
    sbatch_file = SBATCH_DIR / f"{job_name}.sbatch"

    with open(sbatch_file, "w") as f:
        f.write(f"""#!/bin/bash
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t {time_limit}
#SBATCH --mem={memory}
#SBATCH -J {job_name}
#SBATCH -o {BASE_PATH}/misc/logs/merged/%x.log
#SBATCH --container-image="{BASE_PATH}/misc/enroot_images/benchmark_new.sqsh"

set -eu
cd $HOME/gitrepos/cellseg-benchmark

mamba activate seg_postprocessing
pip install trimap
python scripts/seg_postprocessing/vascular_subtyping.py {args.cohort} {seg_method} --condition-col {condition_col}
""")
