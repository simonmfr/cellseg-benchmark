from pathlib import Path

import yaml

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sample_metadata.yaml"
) as f:
    data = yaml.safe_load(f)

Path(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_voronoi"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(
        f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_voronoi/{key}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash
            
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 10:00:00
#SBATCH --mem=128G
#SBATCH -J voronoi_{key}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/voronoi_{key}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/voronoi_{key}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/enroot_images/benchmark_py3_12.sqsh"

cd ~/gitrepos/spatialdata
git pull -q
cd ~/gitrepos/cellseg-benchmark
git pull -q

mamba activate sopa
mkdir -p /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/Negative_Control_Voronoi
python scripts/segmentation/voronoi_segmentation.py {value["path"]} \
 /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/Negative_Control_Voronoi
""")
    f.close()
