#!/usr/bin/env python
"""Compute the ovrlpy VSI and signal maps of one sample (input for the Z integrity metric)."""

import argparse
import logging
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parents[2]))

from cellseg_benchmark import BASE_PATH
from cellseg_benchmark.metrics.ovrl import compute_ovrl

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]: %(message)s")

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("sample", help="Sample name.")
parser.add_argument("data_path", help="Path to the MERSCOPE output folder.")
args = parser.parse_args()

compute_ovrl(args.sample, pathlib.Path(BASE_PATH, "samples", args.sample), args.data_path, logger=logging.getLogger("ovrlpy"))
