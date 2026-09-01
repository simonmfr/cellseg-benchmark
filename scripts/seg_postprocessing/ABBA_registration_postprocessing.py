import argparse
import logging
import warnings
from pathlib import Path

import cellseg_benchmark
import pandas as pd


warnings.filterwarnings("ignore")

logger = logging.getLogger("ABBA_registration_postprocessing")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

def main():
    parser = argparse.ArgumentParser(
        description=("Create points table for registration in QuPath.")
    )
    parser.add_argument("points_file", help="Path to file with mapped points.")
    parser.add_argument("out_dir", help="Path to output directory.")
    parser.add_argument("depth", default=3, help="base path of data structure.")
    args = parser.parse_args()

    mapped = pd.read_csv(args.points_file, index_col=0, sep="\t")
    mapped['atlas_region_id'] = mapped['atlas_region_id'].fillna(-1).astype(int).astype(str)
    ontology_lookup = cellseg_benchmark.ABBA_registration.load_allen_ontology("/Volumes/T7 Shield/registration_qupaths/aging_5_1/Adult Mouse Brain - Allen Brain Atlas V3p1-Ontology.json")
    mapped_coarse = cellseg_benchmark.ABBA_registration.add_depth_restricted_annotation(
        df=mapped,
        ontology_df=ontology_lookup,
        target_depth=args.depth,
        use_acronym=False,
    )
    mapped_coarse.to_csv(Path(args.out_dir) / f"{Path(args.points_file).name}_coarse{args.depth}", sep="\t")