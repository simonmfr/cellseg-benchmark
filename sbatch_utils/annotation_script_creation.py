import sys
from pathlib import Path
from os.path import join, isdir
from os import listdir

sample = sys.argv[1]

path = join("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples", sample, "results")

Path(
    f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatch_annotation/{sample}"
).mkdir(parents=True, exist_ok=True)

for method in listdir(path):
    if isdir(join(path, method, "sdata.zarr")) and isdir(join(path, method, "sdata.explorer")):
        f = open(
            f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatch_annotation/{sample}/{method}.sbatch",
            "w",
        )
        f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 00:20:00
#SBATCH --mem=64G
#SBATCH -J annotation_{sample}_{method}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/annotation_{sample}_{method}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/annotation_{sample}_{method}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sopa_tmp.sqsh"

mamba activate sopa
python /dss/dssfs03/pn52re/pn52re-dss-0001/Git/cellseg-benchmark/scripts/cell_annotation.py \
 {sample} {method} /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark 0
                    """)
        f.close()
