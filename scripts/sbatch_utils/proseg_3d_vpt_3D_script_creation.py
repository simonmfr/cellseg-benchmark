import argparse
from pathlib import Path

import yaml

parser = argparse.ArgumentParser(
    description="Prepare scripts for ProSeg with prior segmentation."
)
parser.add_argument("vpt_flavor", choices=["nuclei", "PolyT", "PolyT_nuclei"], help="vpt flavor.")
parser.add_argument("vpt_dim", choices=["2D", "3D"], help="vpt dimension.")
parser.add_argument("--voxel", default=1, type=int, help="number of z-layers.")
args = parser.parse_args()

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sample_metadata.yaml"
) as f:
    data = yaml.safe_load(f)

Path(
    f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_Proseg_3D_vpt{args.vpt_dim}_{args.vpt_flavor}"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(
        f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_Proseg_3D_vpt{args.vpt_dim}_{args.vpt_flavor}/{key}_{args.voxel}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 18:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=30
#SBATCH -J Proseg_3D_{key}_vpt{args.vpt_dim}_{args.vpt_flavor}_vxl_{args.voxel}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/Proseg_3D_{key}_vpt{args.vpt_dim}_{args.vpt_flavor}_vxl_{args.voxel}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/Proseg_3D_{key}_vpt{args.vpt_dim}_{args.vpt_flavor}_vxl_{args.voxel}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/enroot_images/benchmark_py3_12.sqsh"

cd ~/gitrepos/spatialdata
git pull -q
cd ~/gitrepos/cellseg-benchmark
git pull -q

mamba activate sopa
mkdir -p /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/Proseg_3D_vpt{args.vpt_dim}_{args.vpt_flavor}
python scripts/segmentation/proseg_3D.py {value["path"]} {key} vpt_{args.vpt_dim}_{args.vpt_flavor} \
--voxel-layers {args.voxel} --output-cell-polygon-layers cell-polygons-layers.geojson.gz
""")
    f.close()
