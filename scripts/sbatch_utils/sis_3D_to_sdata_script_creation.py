#!/usr/bin/env python
"""Submit SLURM array job to convert all available SIS outputs to SpatialData zarr stores."""

import subprocess
import pathlib
from cellseg_benchmark import BASE_PATH

BASE = pathlib.Path(BASE_PATH)
SCRIPT = "~/gitrepos/cellseg-benchmark/scripts/seg_postprocessing/sis_to_sdata.py"

samples = sorted([
    p.parent
    for p in (BASE / "samples").glob("*/results/*/sis_out")
    if (p / "cell_by_gene.h5ad").exists() and (p / "cell_polygons.geojson").exists()
])

if not samples:
    print("No samples with sis_out found.")
    exit(0)

print(f"Found {len(samples)} samples.")

paths = " ".join(f'"{s}"' for s in samples)
sbatch = f"""#!/bin/bash
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 02:00:00
#SBATCH --mem=64G
#SBATCH -J SIS_to_sdata
#SBATCH --array=0-{len(samples)-1}
#SBATCH -o {BASE}/misc/logs/outputs/SIS_to_sdata_%a.out
#SBATCH -e {BASE}/misc/logs/errors/SIS_to_sdata_%a.err
#SBATCH --container-image="{BASE}/misc/enroot_images/benchmark.sqsh"

PATHS=({paths})
mamba activate segmentation
python {SCRIPT} ${{PATHS[$SLURM_ARRAY_TASK_ID]}}
"""

sbatch_file = BASE / "misc/sbatches/sbatch_SIS_to_sdata/SIS_to_sdata_array.sbatch"
sbatch_file.parent.mkdir(parents=True, exist_ok=True)
sbatch_file.write_text(sbatch)
print(f"To call: sbatch {sbatch_file}")