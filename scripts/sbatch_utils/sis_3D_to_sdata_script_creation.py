#!/usr/bin/env python
"""Create SLURM array job to convert all available SIS outputs to sdatas."""
import pathlib
import cellseg_benchmark as cb

samples = sorted([
    p.parent
    for p in (pathlib.Path(cb.BASE_PATH) / "samples").glob("*/results/*/sis_out")
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
#SBATCH -o {pathlib.Path(cb.BASE_PATH)}/misc/logs/outputs/SIS_to_sdata_%a.out
#SBATCH -e {pathlib.Path(cb.BASE_PATH)}/misc/logs/errors/SIS_to_sdata_%a.err
#SBATCH --container-image="{pathlib.Path(cb.BASE_PATH)}/misc/enroot_images/benchmark.sqsh"

PATHS=({paths})
mamba activate segmentation
python "~/gitrepos/cellseg-benchmark/scripts/seg_postprocessing/sis_to_sdata.py" ${{PATHS[$SLURM_ARRAY_TASK_ID]}}
"""

sbatch_file = pathlib.Path(cb.BASE_PATH) / "misc/sbatches/sbatch_SIS_to_sdata/SIS_to_sdata_array.sbatch"
sbatch_file.parent.mkdir(parents=True, exist_ok=True)
sbatch_file.write_text(sbatch)
print(f"To call: sbatch {sbatch_file}")