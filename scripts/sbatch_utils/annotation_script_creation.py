import sys
from os import listdir
from os.path import isdir, join
from pathlib import Path

sample = sys.argv[1]

path = join(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples", sample, "results"
)

Path(
    f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_annotation/{sample}"
).mkdir(parents=True, exist_ok=True)

times = {}
for f in listdir(path):
    if f == "Negative_Control_Rastered_5":
        times[f] = "1-12:00:00"
    elif "Baysor" in f:
        times[f] = "02:00:00"
    elif f == "Negative_Control_Rastered_10" or f == "Negative_Control_Voronoi":
        times[f] = "04:00:00"
    else:
        times[f] = "01:00:00"

for method in listdir(path):
    if isdir(join(path, method, "sdata.zarr")):
        f = open(
            f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_annotation/{sample}/{method}.sbatch",
            "w",
        )
        f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t {times[method]}
#SBATCH --mem=25G
#SBATCH --cpus-per-task=5
#SBATCH -J annotation_{sample}_{method}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/annotation_{sample}_{method}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/annotation_{sample}_{method}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/annotation.sqsh"

mamba activate annotation
cd ~/gitrepos/cellseg-benchmark
git pull
python scripts/run_mapmycells.py {sample} {method} /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark
                    """)
        f.close()
