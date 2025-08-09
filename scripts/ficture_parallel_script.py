import argparse
import logging
import warnings
from os.path import join
from pathlib import Path
import socket

from dask.distributed import LocalCluster, Client
import os

# ---- import and run your pipeline (no scheduler= in dask.compute) ----
from cellseg_benchmark.ficture_parallel import run_pipeline_parallel
# 1 worker per core (good for NumPy/Shapely). Tune if needed.

warnings.filterwarnings("ignore")

logger = logging.getLogger("shape_mapping")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

ncores = int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))
cluster = LocalCluster(
    n_workers=ncores, threads_per_worker=1,
    dashboard_address=":8787",         # exposes the dashboard
    processes=True,                     # spawn processes, not threads
    # memory_limit="auto",             # or "3GiB" etc.
    local_directory=os.environ.get("SLURM_TMPDIR","/tmp")+"/dask-worker-space",
)
client = Client(cluster)
logger.info(f"DASHBOARD_NODE: {socket.gethostname()}")
logger.info(f"DASHBOARD: {cluster.dashboard_link}")

parser = argparse.ArgumentParser(
    description="Creates the master sdata, collecting all segmentations from subfolders."
)
parser.add_argument("sample", help="sample name.")
parser.add_argument("data_path", help="Path to data folder with merscope data.")
parser.add_argument(
    "zmode", choices=["z3", "3d"], help="mode of master sdata. Either 'z3' or '3d'."
)
parser.add_argument("data_dir", help="output data folder.")
parser.add_argument(
    "--n_ficture",
    default=21,
    type=int,
    help="Consider Ficture model with n_ficture factors.",
)
args = parser.parse_args()

sdata_path = join(args.data_dir, "samples", args.sample)
seg_methods = [
    method
    for method in os.listdir(join(sdata_path, "results"))
    if os.path.isdir(join(sdata_path, "results", method, "sdata.zarr"))
]

run_pipeline_parallel(
    args=args,
    results_path=Path(sdata_path),
    compute_ficture_methods=seg_methods,
    logger=logger,
    scheduler=None,   # let distributed scheduler handle it
)
client.close()