import argparse
from pathlib import Path

parser = argparse.ArgumentParser(
    description="scripts for merging one method from different samples."
)
parser.add_argument("cohort", help="Cohort name, e.g., 'foxf2'")
parser.add_argument("--age", action="store_true", help="Consider age differentiation")
parser.add_argument(
    "--genotype", action="store_true", help="Consider genotype differentiation"
)
args = parser.parse_args()

base_path = "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark"
sbatch_path = f"{base_path}/misc/sbatches/sbatch_merge_adata"
container_image = f"{base_path}/misc/cellseg_benchmark_2.sqsh"
log_path = f"{base_path}/misc/logs/merged"

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
]

Path(sbatch_path).mkdir(parents=False, exist_ok=True)

for method in methods:
    if method == "Negative_Control_Rastered_5":
        time_limit = "2-00:00:00"
    elif any(keyword in method for keyword in ["Baysor", "Cellpose"]) or method in [
        "Negative_Control_Rastered_10",
        "Negative_Control_Voronoi",
    ]:
        time_limit = "12:00:00"
    else:
        time_limit = "06:00:00"

    memory = "700G" if "Negative_Control" in method else "500G"

    with open(f"{sbatch_path}/{args.cohort}_{method}.sbatch", "w") as f:
        f.write(f"""#!/bin/bash
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t {time_limit}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem={memory}
#SBATCH -J merge_adata_{args.cohort}_{method}
#SBATCH -o {log_path}/%x.log
#SBATCH --container-image="{container_image}"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

cd ~/gitrepos/cellseg-benchmark
git pull
mamba activate cellseg_benchmark
python scripts/merge_adata.py {args.cohort} {method} {"--age" if args.cohort == "aging" else "--genotype"}
""")
