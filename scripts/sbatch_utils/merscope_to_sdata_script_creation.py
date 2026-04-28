#!/usr/bin/env python
import argparse
import pathlib
import yaml

BASE_PATH = "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark"

parser = argparse.ArgumentParser(description="Create a SLURM array sbatch for a given cohort to convert Merscope outputs to sdata.")
parser.add_argument("cohort", help="Cohort name (e.g. aging).")
parser.add_argument("--segmentation", choices=["cellpose", "watershed", "both"], default="both")
parser.add_argument("--explorer", action="store_true", help="Write 10X Xenium explorer files (cellpose only).")
args = parser.parse_args()

with open(pathlib.Path(BASE_PATH) / "misc/sample_metadata.yaml") as f:
    metadata = yaml.safe_load(f)

seg_configs = [
    ("cellpose",  "path",           "Cellpose_1_Merlin"),
    ("watershed", "path_watershed", "Watershed_Merlin"),
]
if args.segmentation != "both":
    seg_configs = [c for c in seg_configs if c[0] == args.segmentation]

jobs = []
for key, val in metadata.items():
    if key.split("_")[0] != args.cohort:
        continue
    for seg, path_key, result_dir in seg_configs:
        if path_key not in val:
            continue
        data_path = pathlib.Path(val[path_key])
        if (data_path / "cell_by_gene.csv").exists():
            save_path = pathlib.Path(BASE_PATH) / "samples" / key / "results" / result_dir
            jobs.append((str(data_path), str(save_path), seg))

if not jobs:
    print(f"No valid samples found for cohort '{args.cohort}'.")
    exit(0)

data_paths = " ".join(f'"{dp}"' for dp, _, __ in jobs)
save_paths = " ".join(f'"{sp}"' for _, sp, __ in jobs)
seg_flags  = " ".join(f'"{sf}"' for _, __, sf in jobs)

sbatch = f"""#!/bin/bash
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 12:00:00
#SBATCH --mem=250G
#SBATCH -J merscope_{args.cohort}
#SBATCH --array=0-{len(jobs)-1}
#SBATCH -o {pathlib.Path(BASE_PATH)}/misc/logs/outputs/merscope_to_sdata_{args.cohort}_%a.out
#SBATCH -e {pathlib.Path(BASE_PATH)}/misc/logs/errors/merscope_to_sdata_{args.cohort}_%a.err
#SBATCH --container-image="{pathlib.Path(BASE_PATH)}/misc/enroot_images/benchmark.sqsh"

DATA_PATHS=({data_paths})
SAVE_PATHS=({save_paths})
SEG_FLAGS=({seg_flags})

mamba activate segmentation
mkdir -p "${{SAVE_PATHS[$SLURM_ARRAY_TASK_ID]}}"

EXPLORER_FLAG=""
if [ "${{SEG_FLAGS[$SLURM_ARRAY_TASK_ID]}}" = "cellpose" ] && [ "{str(args.explorer).lower()}" = "true" ]; then
    EXPLORER_FLAG="--explorer"
fi

python ~/gitrepos/cellseg-benchmark/scripts/seg_postprocessing/merscope_to_sdata.py \\
  "${{DATA_PATHS[$SLURM_ARRAY_TASK_ID]}}" \\
  "${{SAVE_PATHS[$SLURM_ARRAY_TASK_ID]}}" \\
  --segmentation "${{SEG_FLAGS[$SLURM_ARRAY_TASK_ID]}}" \\
  $EXPLORER_FLAG
"""

sbatch_dir = pathlib.Path(BASE_PATH) / "misc/sbatches/sbatch_merscope_to_sdata"
sbatch_dir.mkdir(parents=True, exist_ok=True)
sbatch_file = sbatch_dir / f"{args.cohort}_{args.segmentation}.sbatch"
sbatch_file.write_text(sbatch)
print(f"[{args.cohort}] {len(jobs)} jobs. Call: sbatch {sbatch_file}")