from pathlib import Path

base_path = "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark"
sbatch_path = f"{base_path}/misc/sbatches/sbatch_merge_adata"
container_image = f"{base_path}/misc/cellseg_benchmark.sqsh"
log_path = f"{base_path}/misc/logs/outputs"

methods = [
    "Baysor_2D_Cellpose_2_DAPI_Transcripts_0.8",
    "Cellpose_1_nuclei_model",
    "Baysor_2D_Cellpose_2_DAPI_PolyT_0.2",
    "Baysor_2D_Cellpose_1_DAPI_Transcripts_0.8",
    "Negative_Control_Rastered_5",
    "Cellpose_2_DAPI_PolyT",
    "Proseg_pure",
    "Proseg_Cellpose_1_DAPI_Transcripts",
    "Baysor_2D_Cellpose_1_DAPI_Transcripts_0.2",
    "Baysor_2D_Cellpose_1_DAPI_PolyT_0.2",
    "Baysor_2D_Cellpose_2_DAPI_Transcripts_0.2",
    "Cellpose_1_DAPI_PolyT",
    "Negative_Control_Rastered_10",
    "Negative_Control_Voronoi",
    "Proseg_Cellpose_2_DAPI_Transcripts",
    "Proseg_Cellpose_1_nuclei_model",
    "Cellpose_2_DAPI_Transcripts",
    "Cellpose_1_Merlin",
    "Negative_Control_Rastered_25",
    "Cellpose_1_DAPI_Transcripts",
    "Baysor_2D_Cellpose_2_DAPI_PolyT_0.8",
    "Baysor_2D_Cellpose_1_nuclei_model_1.0",
    "Proseg_Cellpose_2_DAPI_PolyT",
    "Baysor_2D_Cellpose_1_DAPI_PolyT_0.8",
    "Proseg_Cellpose_1_DAPI_PolyT",
]

Path(sbatch_path).mkdir(parents=False, exist_ok=True)

for method in methods:
    if method == "Negative_Control_Rastered_5":
        time_limit = "2-00:00:00"
    elif any(keyword in method for keyword in ["Baysor", "Cellpose"]) or \
         method in ["Negative_Control_Rastered_10", "Negative_Control_Voronoi"]:
        time_limit = "12:00:00"
    else:
        time_limit = "04:00:00"
    
    with open(f"{sbatch_path}/{method}.sbatch", "w") as f:
        f.write(f"""#!/bin/bash
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t {time_limit}
#SBATCH --mem=200G
#SBATCH -J merge_adata_{method}
#SBATCH -o {log_path}/%x.log
#SBATCH --container-image="{container_image}"
cd ~/gitrepos/cellseg-benchmark
git pull
mamba activate cellseg_benchmark
python scripts/merge_adata.py {method}
""")