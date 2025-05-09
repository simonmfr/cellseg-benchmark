import pandas as pd
import logging
from cellseg_benchmark.adata_utils import *

if __name__ == "__main__":
    logger = logging.getLogger("integration_harmony")
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
    logger.addHandler(handler)

    adatas = [("foxf2_s1_r0", read_h5ad("/Users/jonasflor/Downloads/foxf2_s1_r0.h5ad")),
              ("foxf2_s1_r1", read_h5ad("/Users/jonasflor/Downloads/foxf2_s1_r1.h5ad"))]
    adata = adatas[0][1]
    obs = pd.read_csv("/Users/jonasflor/Downloads/adata_obs_annotated_1_0.csv")[["cell_type_mmc_incl_low_quality_revised",
           "cell_type_mmc_incl_low_quality_clusters",
           "cell_type_mmc_incl_low_quality",
           "cell_type_mmc_incl_mixed_revised",
           "cell_type_mmc_incl_mixed_clusters",
           "cell_type_mmc_incl_mixed",
           "cell_type_mmc_raw_revised",
           "cell_type_mmc_raw_clusters",
           "cell_type_mmc_raw",
           "cell_id"]]
    new_obs = adata.obs.merge(
        obs, how="left", left_index=True, right_on="cell_id"
    )
    for col in new_obs.columns:
        if isinstance(new_obs[col].dtype, pd.CategoricalDtype):
            new_obs[col] = new_obs[col].cat.add_categories("Low-Read-Cells")
        new_obs[col] = new_obs[col].fillna("Low-Read-Cells")
    new_obs.index = adata.obs.index
    adatas[0][1].obs = new_obs

    adata = adatas[1][1]
    obs = pd.read_csv("/Users/jonasflor/Downloads/adata_obs_annotated_1_1.csv")[
        ["cell_type_mmc_incl_low_quality_revised",
         "cell_type_mmc_incl_low_quality_clusters",
         "cell_type_mmc_incl_low_quality",
         "cell_type_mmc_incl_mixed_revised",
         "cell_type_mmc_incl_mixed_clusters",
         "cell_type_mmc_incl_mixed",
         "cell_type_mmc_raw_revised",
         "cell_type_mmc_raw_clusters",
         "cell_type_mmc_raw",
         "cell_id"]]
    new_obs = adata.obs.merge(
        obs, how="left", left_index=True, right_on="cell_id"
    )
    for col in new_obs.columns:
        if isinstance(new_obs[col].dtype, pd.CategoricalDtype):
            new_obs[col] = new_obs[col].cat.add_categories("Low-Read-Cells")
        new_obs[col].fillna("Low-Read-Cells", inplace=True)
    new_obs.index = adata.obs.index
    adatas[1][1].obs = new_obs

    tmp = merge_adatas_deb(adatas, logger=logger, do_qc=False, save_path="/Users/jonasflor/Desktop/debug_pics")
    tmp = filter_cells(tmp, save_path="/Users/jonasflor/Desktop/debug_pics", logger=logger)
    tmp = filter_genes(tmp, save_path="/Users/jonasflor/Desktop/debug_pics", logger=logger)
    tmp = normalize(tmp, "/Users/jonasflor/Desktop/debug_pics", logger=logger)
    dimensionality_reduction(tmp, "/Users/jonasflor/Desktop/debug_pics", logger=logger)
    tmp = integration_harmony(tmp, "full_name", "/Users/jonasflor/Desktop/debug_pics", logger=logger)