import argparse
import json
import logging
import re
from collections import defaultdict
from os.path import join
from pathlib import Path

import numpy as np
import pandas as pd
from spatialdata import read_zarr, transform
from tqdm import tqdm

from cellseg_benchmark.metrics.ovrl import compute_ovrl, compute_mean_vsi_per_polygon, plot_vsi_overview

logger = logging.getLogger("integrate_adatas")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

parser = argparse.ArgumentParser(
    description="Run ovrlpy über eine cohorte."
)
parser.add_argument("cohort", help="Cohort name, e.g., 'foxf2'")
parser.add_argument(
    "--age",
    default=False,
    action="store_true",
    help="Age information is available. Otherwise, assume age: 6m.",
)
parser.add_argument(
    "--genotype",
    default=False,
    action="store_true",
    help="Genotype information is available. Otherwise, assume WT.",
)
args = parser.parse_args()

base_path = Path("/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark")
samples_path = base_path / "samples"
save_path = base_path / "metrics" / args.cohort
save_path.mkdir(parents=True, exist_ok=True)

with open(base_path / "sample_paths.json") as f:
    sample_paths_file = json.load(f)

excluded_samples = {
    "foxf2": ["foxf2_s2_r0", "foxf2_s3_r0", "foxf2_s3_r1"],
    "aging": [],
}.get(args.cohort, [])

logger.info("Loading data...")
sdata_list = []
available_names = set()

for sample_dir in samples_path.glob(f"{args.cohort}*"):  # restriction to cohort folders
    if sample_dir.name in excluded_samples:
        continue
    sdata = read_zarr(sample_dir / "sdata_z3.zarr", selection=("shapes", "points"))
    current_names = ["_".join(k.split("_")[1:]) for k in sdata.shapes.keys()]
    available_names.update(current_names)
    sdata_list.append((sample_dir.name, sdata))

logger.info("Compute ovrlpy outputs...")
compute_ovrl(sdata_list, str(base_path), logger=logger)

metrics = defaultdict(list)
for sample, sdata in tqdm(sdata_list):
    logger.info(f"Compute mean vsi per polygon for all methods in {sample}...")
    transform_matrix = np.loadtxt(
        join(
            sample_paths_file[sample],
            "images",
            "micron_to_mosaic_pixel_transform.csv",
        )
    ).reshape(3, 3)
    integrity_matrix = np.load(join(str(base_path), "samples", sample, 'vertical_doublets_ovrlpy_output.npz'))['integrity']
    signal_matrix = np.load(join(str(base_path), "samples", sample, 'vertical_doublets_ovrlpy_output.npz'))['signal']
    for boundary_name in sdata.shapes.keys():
        logger.info(f"Computation for {boundary_name}...")
        png_path = samples_path / sample / "results" / "_".join(boundary_name.split("_")[1:]) / "ovrlpy_overview_plot.png"
        boundary = transform(sdata[boundary_name], to_coordinate_system="micron")
        tmp = compute_mean_vsi_per_polygon(integrity_matrix, boundary, transform_matrix)
        samples = sample.split("_")
        cohort, slide, region = samples[0], samples[1], samples[2]
        tmp["cohort"] = cohort
        tmp["slide"] = slide
        tmp["region"] = region
        tmp["sample"] = sample
        if args.age:
            # Assume, that all region names contain months in the end (e.g. region_0-WT279_12m)
            tmp["age"] = int(
                re.split(r"[_\-]", sample_paths_file[sample].split("/")[-1])[-1].split(
                    "m"
                )[0]
            )
        else:
            tmp["age"] = 6
            # Requires region names to follow example: region_1-KO885
        tmp["genotype"] = re.search(
            r"region_\d+-([A-Za-z_]*?)(?=\d)", sample_paths_file[sample]
        ).group(1)
        tmp["condition"] = tmp["genotype"] + "_" + tmp["age"].astype(str)
        metrics["_".join(boundary_name.split("_"))].append(tmp)
        plot_vsi_overview(integrity_map=integrity_matrix, signal_map = signal_matrix, boundaries_aligned = boundary, vsi_mean=tmp["mean_integrity"].values, sample_name = sample, png_path=png_path)

logger.info(f"Writing metrics to {save_path}")
for method, metric in tqdm(metrics.items()):
    pd.concat(metric).to_csv(save_path / f"mean_vsi_{method}.csv", index=False)
