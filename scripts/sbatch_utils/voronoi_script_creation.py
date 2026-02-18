from pathlib import Path
import yaml

YAML = "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sample_metadata.yaml"
SBATCH_DIR = "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_visium"
OUT_DIR = "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{k}/results/Negative_Control_Visium"

with open(YAML) as f:
    data = yaml.safe_load(f)

Path(SBATCH_DIR).mkdir(exist_ok=True)

for k, v in data.items():
    out = OUT_DIR.format(k=k)

    text = f"""#!/bin/bash
#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 10:00:00
#SBATCH --mem=128G
#SBATCH -J voronoi_{k}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/voronoi_{k}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/voronoi_{k}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/enroot_images/benchmark.sqsh"

mamba activate segmentation
mkdir -p {out}

"$CONDA_PREFIX/bin/time" -v \
  -o "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/voronoi_{k}_$(date +%Y%m%d_%H%M%S).time" \
python ~/gitrepos/cellseg-benchmark/scripts/segmentation/voronoi_segmentation.py {v["path"]} {out}
"""

    Path(SBATCH_DIR, f"{k}.sbatch").write_text(text)