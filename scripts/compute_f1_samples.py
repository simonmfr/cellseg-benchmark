import argparse
import logging
import warnings
from os import listdir
from os.path import exists, join
from re import split

import pandas as pd
from matplotlib import pyplot as plt
from scanpy import read_h5ad
import seaborn as sns

from cellseg_benchmark._constants import factor_to_celltype, true_cluster, column_order, index_order
from cellseg_benchmark.metrics import compute_f1

warnings.filterwarnings("ignore")

logger = logging.getLogger("compute_f1")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
logger.addHandler(handler)

parser = argparse.ArgumentParser(description="Compute F1 statistics based on Ficture.")
parser.add_argument("method", help="method name.")
parser.add_argument("cohort", help="Cohort name.")
parser.add_argument("data", choices=["area", "variance", "means"], help="Cohort name.")
parser.add_argument(
    "--correct_celltypes", action="store_true", help="Compute F1 for correct cell type mapping."
)
parser.add_argument(
    "--weighted", action="store_true", help="If data is weighted."
)
parser.add_argument(
    "--celltype_name", default="cell_type_mmc_raw_revised", help="Compute F1 for correct cell type mapping."
)
parser.add_argument(
    "--flavor", default="f1", choices=["f1", "macro", "micro", "all"], help="Which flavor to compute F1 for."
)
parser.add_argument(
    "subset", nargs=argparse.REMAINDER, help="List of celltypes to compute F1 for."
)
args = parser.parse_args()

base_path = "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/"
data_path = join(base_path, "analysis", args.cohort, args.method)

adata = read_h5ad(join(data_path, "adatas", "adata_integrated.h5ad.gz"))

def rename(colnames):
    mapper = {}
    for col in colnames:
        mapper[col] = factor_to_celltype[split("_", col)[1]]
    return mapper

if args.correct_celltypes:
    correct_celltypes = {}
    for i in factor_to_celltype.keys():
        correct_celltypes[factor_to_celltype[i]] = true_cluster[factor_to_celltype[i]]
else:
    correct_celltypes = None

if not args.subset:
    subset = None

obsm_key = f"ficture_{args.data}{'_weight' if args.weighted else ''}"

general_stats_dic = {}
for sample in listdir(join(base_path, "samples")):
    if sample.startswith(args.cohort) and exists(join(base_path, "samples", sample, "results", "Ficture", "general_stats.csv")):
        tmp = pd.read_csv(join(base_path, "samples", sample, "results/Ficture/general_stats.csv"), index_col=0)
        tmp.rename(columns=factor_to_celltype, inplace=True)
        general_stats_dic[sample] = tmp

data = {}
for key in general_stats_dic.keys():
    data[key] = adata[adata.obs['sample'] == key].obsm[obsm_key].copy()
    data[key].rename(columns=rename(data[key].columns), inplace=True)
    data[key]['celltype'] = adata[adata.obs['sample'] == key].obs[args.celltype_name].values

f1 = compute_f1(data, general_stats=general_stats_dic, flavor=args.flavor, correct_celltypes=correct_celltypes, subset=subset, weighted=args.weighted)
f1.to_csv(join(base_path, "metrics", args.cohort, f"{args.method}_f1{'_weighted' if args.weighted else ''}_{'celltypes' if args.correct_celltypes else 'matrix'}.csv"))
if args.correct_celltypes:
    sns.set_theme(rc={"figure.figsize": (20, 16)})
    sns.barplot(
        data=f1
    )
    plt.xticks(rotation=90)
    plt.ylim(0, 1)
    plt.savefig(join(base_path, "metrics", args.cohort, "ficture", f"{args.method}_f1_{args.data}{'_weighted' if args.weighted else ''}_barplot.png"))
else:
    data = f1.astype(float)
    data.index = pd.CategoricalIndex(data.index, categories=index_order)
    data.sort_index(level=0, inplace=True)
    data = data[column_order]
    if args.weighted:
        sns.set_theme(rc={"figure.figsize": (20, 16)})
        sns.heatmap(data, cmap="YlOrRd", annot=True)
        plt.savefig(join(base_path, "metrics", args.cohort, "ficture"
                         f"{args.method}_f1_{args.data}{'_weighted' if args.weighted else ''}_heatmap.png"))
    else:
        sns.set_theme(rc={"figure.figsize": (20, 16)})
        sns.heatmap(data, fmt=".3f", cmap="YlOrRd", vmin=0, vmax=1, annot=True)
        plt.savefig(join(base_path, "metrics", args.cohort, "ficture", f"{args.method}_f1_{args.data}{'_weighted' if args.weighted else ''}_heatmap.png"))