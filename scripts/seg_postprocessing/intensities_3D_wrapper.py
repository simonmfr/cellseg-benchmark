import argparse
import logging
import pathlib
import os
import subprocess
import sys

logger = logging.getLogger("intensities_3D_wrapper")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

def main():
    parser = argparse.ArgumentParser(
        description=(
            "wrapper for intensities_3D to apply sample-wise."
        )
    )
    parser.add_argument("sample", type=str, help="Sample name.")
    parser.add_argument("data_path", help="Path to merfish output folder.")
    parser.add_argument("--method_names", type=bool, default=False, help="Method names.")
    parser.add_argument("--method", nargs=argparse.REMAINDER, type=str, help="If --method_names then provide a list of methods with flavors, e.g. Proseg_3D_Cellpose_1_nuclei_model. Otherwise provide list of segmentation algorithms, e.g. Proseg_3D. If not provided all intensities of available 3D segmentation methods will be computed.")
    args = parser.parse_args()

    BASE_PATH = "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark"
    methods_3D = ["Proseg_3D", "vpt_3D", "Watershed_Merlin", "SIS"]
    sample_path = pathlib.Path(BASE_PATH) / "samples" / args.sample

    methods = []
    if args.method_names:
        logger.info("Check existence of provided methods")
        for method_name in args.method:
            if not os.path.exists(sample_path / "results" / method_name / "sdata.zarr"):
                logger.error(f"Method {method_name} does not exist or doesn't contain a sdata.zarr file.")
            else:
                methods.append(method_name)
    else:
        for method_name in os.listdir(sample_path / "results"):
            if os.path.exists(sample_path / "results" / method_name / "sdata.zarr") and any([method_name.startswith(x) for x in methods_3D]):
                logger.info(f"Method {method_name} identified for computation.")
                methods.append(method_name)

    if len(methods) == 0:
        logger.info("No methods were provided or identified.")
        return

    intensities_3D_script = pathlib.Path(__file__).parent / "intensities_3D.py"
    for method in methods:
        cmd = [sys.executable, str(intensities_3D_script), str(sample_path / "results" / method), args.data_path, "--method", method]

        exit_code = subprocess.run(cmd, stdout=sys.stdout, stderr=sys.stderr)
        if exit_code.returncode != 0:
            logger.error(f"Intensity computation in 3D failed for method {method} with exit code {exit_code.returncode}")
        else:
            logger.info(f"Intensity computation in 3D succeeded for method {method}")
    logger.info("Done.")

if __name__ == "__main__":
    main()
