import argparse
import logging

# TODO: ADD TO CONTAINER
import subprocess
import warnings
from pathlib import Path

import pandas as pd
from spatialdata import read_zarr

from cellseg_benchmark.metrics.boundary_based import (
    count_assigned_transcripts,
)

subprocess.run(
    ["conda", "install", "-n", "cellseg_benchmark", "-c", "conda-forge", "rtree", "-y"]
)

warnings.filterwarnings("ignore")

# Logger setup
logger = logging.getLogger("compute_assigned_transcripts")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

# CLI args
parser = argparse.ArgumentParser(
    description="Compute fraction of cell-assigned transcripts for all segmentaiton methods from a master-sdata."
)
parser.add_argument("sample", help="Sample, e.g., 'foxf2_s2_r1'")
args = parser.parse_args()

# Paths
base_path = Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")
sample_path = base_path / "samples" / args.sample
save_path = base_path / "samples" / args.sample / "misc"
save_path.mkdir(parents=True, exist_ok=True)

# Load sdata
logger.info("Loading sdata...")
sdata = read_zarr(sample_path / "sdata_z3.zarr")

logger.info("Computing assigned transcripts...")
assigned_transcripts_results = count_assigned_transcripts(
    sdata, sdata_transcripts_key=next(iter(sdata.points.keys()))
)

# Flatten to df and save
df = pd.DataFrame.from_dict(assigned_transcripts_results, orient="index")
df.reset_index(inplace=True)
df = df.rename(columns={"index": "seg_method"})
logger.info("Saving results...")
df.to_csv(os.path.join(save_path, "assigned_transcripts_results.csv"), index=False)
logger.info("Done.")
