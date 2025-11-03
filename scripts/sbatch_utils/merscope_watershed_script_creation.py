from pathlib import Path

data = {}
data["foxf2_s1_r0"] = (
    "/dss/dssfs03/pn52re/pn52re-dss-0000/202402-Foxf2/merfish_output/watershed/20240229_mousebrain-Slide01-ws-WT-ECKO/region_0-ECKO000"
)
data["foxf2_s1_r1"] = (
    "/dss/dssfs03/pn52re/pn52re-dss-0000/202402-Foxf2/merfish_output/watershed/20240229_mousebrain-Slide01-ws-WT-ECKO/region_1-WT000"
)
data["foxf2_s4_r0"] = (
    "/dss/dssfs03/pn52re/pn52re-dss-0000/202402-Foxf2/merfish_output/watershed/20240304_mousebrain-Slide04-ws-ECKO-PCKO/region_0-PCKO421"
)
data["foxf2_s4_r1"] = (
    "/dss/dssfs03/pn52re/pn52re-dss-0000/202402-Foxf2/merfish_output/watershed/20240304_mousebrain-Slide04-ws-ECKO-PCKO/region_1-ECKO557"
)

Path(
    "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_merscope_watershed"
).mkdir(parents=False, exist_ok=True)
for key, value in data.items():
    f = open(
        f"/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/sbatches/sbatch_merscope_watershed/{key}.sbatch",
        "w",
    )
    f.write(f"""#!/bin/bash

#SBATCH -p lrz-cpu
#SBATCH --qos=cpu
#SBATCH -t 12:00:00
#SBATCH --mem=300G
#SBATCH -J merscope_watershed_{key}
#SBATCH -o /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/outputs/merscope_watershed_{key}.out
#SBATCH -e /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/logs/errors/merscope_watershed_{key}.err
#SBATCH --container-image="/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/enroot_images/benchmark_py3_12.sqsh"

cd ~/gitrepos/spatialdata
git pull -q
cd ~/gitrepos/cellseg-benchmark
git pull -q

mamba activate sopa
mkdir -p /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/Watershed_Merlin
python scripts/segmentation/merscope_watershed_sdata.py {value} \
 /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/samples/{key}/results/Watershed_Merlin
""")
    f.close()
