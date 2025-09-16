from pathlib import Path

cohorts = ["foxf2"]  # "aging"
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
    "Proseg_pure",
    "Proseg_Cellpose_1_DAPI_Transcripts",
    "Proseg_Cellpose_1_DAPI_PolyT",
    "Proseg_Cellpose_1_nuclei_model",
    "Proseg_Cellpose_2_DAPI_PolyT",
    "Proseg_Cellpose_2_DAPI_Transcripts",
    "Negative_Control_Rastered_5",
    "Negative_Control_Rastered_10",
    "Negative_Control_Rastered_25",
    "Negative_Control_Voronoi",
    "vpt_2D_DAPI_PolyT",
    "vpt_2D_DAPI_nuclei",
    "vpt_2D_DAPI_PolyT_nuclei",
    "vpt_3D_DAPI_PolyT",
    "vpt_3D_DAPI_nuclei",
    "vpt_3D_DAPI_PolyT_nuclei",
    "ComSeg",
]

base_path = "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark"
sbatch_dir = Path(f"{base_path}/misc/sbatches/sbatch_vasc_subty")
log_path = f"{base_path}/misc/logs/merged"
container_image = f"{base_path}/misc/cellseg_benchmark_2.sqsh"

sbatch_dir.mkdir(parents=True, exist_ok=True)

for cohort in cohorts:
    for method in methods:
        job_name = f"vasc_subty_{cohort}_{method}"
        sbatch_file = sbatch_dir / f"{job_name}.sbatch"

        with open(sbatch_file, "w") as f:
            f.write(f"""#!/bin/bash
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 2:00:00
#SBATCH --mem=32G
#SBATCH -J {job_name}
#SBATCH -o {log_path}/%x.log
#SBATCH --container-image="{container_image}"

set -eu
cd ~/gitrepos/cellseg-benchmark
git pull
mamba activate cellseg_benchmark
python scripts/vascular_subtyping.py {cohort} {method}
""")
