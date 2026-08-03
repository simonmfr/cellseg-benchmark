#!/usr/bin/env python
import pathlib
import yaml

BASE_PATH = pathlib.Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")

with open(
    f"{BASE_PATH}/misc/sample_metadata.yaml"
) as f:
    data = yaml.safe_load(f)

pathlib.Path(
    f"{BASE_PATH}/misc/sbatches/sbatch_Intensities_3D"
).mkdir(parents=False, exist_ok=True)

for key, value in data.items():
    f = open(
        f"{BASE_PATH}/misc/sbatches/sbatch_Intensities_3D/{key}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 12:00:00
#SBATCH --mem=150G
#SBATCH -J intensities_3D_{key}
#SBATCH -o {BASE_PATH}/misc/logs/outputs/intensities_3D_{key}.out
#SBATCH -e {BASE_PATH}/misc/logs/errors/intensities_3D_{key}.err
#SBATCH --container-image="{BASE_PATH}/misc/enroot_images/benchmark_new.sqsh"
            
mamba activate segmentation
python $HOME/gitrepos/cellseg-benchmark/scripts/seg_postprocessing/intensities_3D_wrapper.py {key} {value["path"]}
""")
    f.close()
