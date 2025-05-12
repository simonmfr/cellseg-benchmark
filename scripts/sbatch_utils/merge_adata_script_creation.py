from pathlib import Path

Path(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_merge_adata"
).mkdir(parents=False, exist_ok=True)

methods = ["Baysor_2D_Cellpose_2_DAPI_Transcripts_0.8", "Cellpose_1_nuclei_model", "Baysor_2D_Cellpose_2_DAPI_PolyT_0.2",
           "Baysor_2D_Cellpose_1_DAPI_Transcripts_0.8", "Negative_Control_Rastered_5",
           "Cellpose_2_DAPI_PolyT", "Proseg_pure", "Proseg_Cellpose_1_DAPI_Transcripts",
           "Baysor_2D_Cellpose_1_DAPI_Transcripts_0.2", "Baysor_2D_Cellpose_1_DAPI_PolyT_0.2",
           "Baysor_2D_Cellpose_2_DAPI_Transcripts_0.2", "Cellpose_1_DAPI_PolyT",
           "Negative_Control_Rastered_10", "Negative_Control_Voronoi", "Proseg_Cellpose_2_DAPI_Transcripts",
           "Proseg_Cellpose_1_nuclei_model", "Cellpose_2_DAPI_Transcripts", "Cellpose_1_Merlin",
           "Negative_Control_Rastered_25", "Cellpose_1_DAPI_Transcripts", "Baysor_2D_Cellpose_2_DAPI_PolyT_0.8",
           "Baysor_2D_Cellpose_1_nuclei_model_1.0", "Proseg_Cellpose_2_DAPI_PolyT",
           "Baysor_2D_Cellpose_1_DAPI_PolyT_0.8", "Proseg_Cellpose_1_DAPI_PolyT"]

times = {}
for method in methods:
    if "Baysor" in method:
        times[method] = "07:00:00"
    elif "Negative_Control_Rastered_5" == method:
        times[method] = "1-00:00:00"
    elif method in ["Negative_Control_Rastered_10", "Negative_Control_Voronoi"]:
        times[method] = "12:00:00"
    else:
        times[method] = "04:00:00"

for method in methods:
    f = open(
        f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_merge_adata/{method}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t {times[method]}
#SBATCH --mem=200G
#SBATCH -J merge_adata_{method}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/merge_adata_{method}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/merge_adata_{method}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/cellseg_benchmark.sqsh"


cd ~/gitrepos/cellseg-benchmark
git pull
mamba activate cellseg_benchmark
python scripts/merge_adata.py {method}
""")
    f.close()