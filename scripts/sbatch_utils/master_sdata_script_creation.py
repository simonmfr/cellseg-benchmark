from pathlib import Path

import yaml

with open(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sample_metadata.yaml"
) as f:
    data = yaml.save_load(f)

Path(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_master_sdata"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(
        f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_master_sdata/{key}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=150G
#SBATCH -J master_sdata_{key}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/master_sdata_{key}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/master_sdata_{key}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/cellseg_benchmark_2.sqsh"


cd ~/gitrepos/cellseg-benchmark
git pull
mamba activate cellseg_benchmark
python scripts/master_sdata.py {key} {value["path"]} z3 /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark \
--genotype {value["genotype"]} \
--age_months {value["age_months"]} \
--run_date {value["run_date"]} \
--animal_id {value["animal_id"]} \
--organism {value["organism"]} \
--slide {value["slide"]} \
--region {value["region"]} \
--cohort {value["cohort"]} \
""")
    f.close()
