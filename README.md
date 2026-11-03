# Cell segmentation benchmark

## Repository structure
- `cellseg_benchmark/`: package
- `scripts/`: segmentation and metric scripts; `scripts/sbatch_utils/` generates per-method sbatch scripts
- `notebooks/`: development and analysis
- `configs/`: configuration files, e.g. VPT
- `archive/`: symlink to raw MERSCOPE data on DSS
- `data/`: symlink to processed data on DSS

## Development
[ruff](https://docs.astral.sh/ruff/) formats and lints via a pre-commit hook (`.pre-commit-config.yaml`). Rules are set in `pyproject.toml` ([rule list](https://docs.astral.sh/ruff/rules/)).
```
pip install ruff pre-commit
pre-commit install
pre-commit run --all-files  # optional: check whole codebase
```
The hook checks changed files on every commit. Fix reported errors, re-add, commit again. Skip once with `git commit --no-verify`. Manual run: `ruff format && ruff check --fix`.

## Container setup (enroot, LRZ)
Run enroot on a compute node. Images in `misc/enroot_images/testing/`: `segmentation.sqsh` (envs `segmentation`, `SIS_segmentation`), `analysis.sqsh` (envs `seg_postprocessing`, `analysis`).

### Base image
```
enroot import docker://ubuntu
enroot create --name <container_name> ubuntu.sqsh

enroot start -r -w <container_name>
apt-get update && apt-get install -y vim curl git nano
echo "HOME=/home/ubuntu" >> /etc/profile
exit

enroot start -w <container_name>
cd /tmp
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
cat > /home/ubuntu/activate_env.sh <<'EOL'
eval "$(/home/ubuntu/miniforge3/bin/conda shell.bash hook)"
[ -f /home/ubuntu/miniforge3/etc/profile.d/mamba.sh ] && . /home/ubuntu/miniforge3/etc/profile.d/mamba.sh
EOL
echo "BASH_ENV=/home/ubuntu/activate_env.sh" >> /etc/environment
source ~/.bashrc  # should activate conda base
exit
```
Make project folders visible by adding `<outside_path> <inside_path> none x-create=dir,bind` to the container's `/etc/fstab` ([format](https://github.com/NVIDIA/enroot/blob/main/doc/image-format.md)).

### Repository clone and auto-pull
Jobs run scripts from `/home/ubuntu/gitrepos/cellseg-benchmark`, so each image needs this clone. Scripts import `cellseg_benchmark` from it; segmentation envs don't pip-install the package. Clone as user so pulls can write:
```
enroot start -w <container_name>
mkdir -p ~/gitrepos && git clone https://github.com/simonmfr/cellseg-benchmark.git ~/gitrepos/cellseg-benchmark
exit
```
As root, add an auto-pull to `/etc/rc`. It runs at every job start and isn't saved to the image.
```
enroot start -r -w <container_name>
cat > /tmp/rc_insert.sh <<'EOF'
# BEGIN cellseg-benchmark auto-update
REPO="/home/ubuntu/gitrepos/cellseg-benchmark"
if command -v git >/dev/null 2>&1 && [ -d "$REPO/.git" ]; then
  export GIT_TERMINAL_PROMPT=0
  git -C "$REPO" pull --ff-only -q || true
fi
# END cellseg-benchmark auto-update
EOF
grep -q "BEGIN cellseg-benchmark auto-update" /etc/rc || \
  { awk 'NR==1 {print; system("cat /tmp/rc_insert.sh"); next} {print}' /etc/rc > /tmp/rc.new && cat /tmp/rc.new > /etc/rc; }
chmod +x /etc/rc
exit
```

### Analysis env
As user, with DEA (R) and Jupyter support:
```
mamba create -n cellseg_benchmark python=3.12 -c conda-forge -c bioconda r-base rpy2 anndata2ri zlib cmake -y
mamba activate cellseg_benchmark
pip install -e ~/gitrepos/cellseg-benchmark
pip install setuptools==80.0.0
cd ~/gitrepos
git clone https://github.com/AllenInstitute/cell_type_mapper.git
git clone https://github.com/jonas2612/spatialdata.git
pip install ./cell_type_mapper ./spatialdata
```
R packages for DEA, in `R` from the same env:
```R
install.packages(c("dplyr", "data.table", "Matrix", "BiocManager"))
BiocManager::install(c("edgeR", "variancePartition", "limma", "BiocParallel", "MAST", "SummarizedExperiment"))
```
Jupyter kernel: run `python -m ipykernel install --user --name cellseg_benchmark`, then add `"env": {"R_HOME": "/home/ubuntu/miniforge3/envs/cellseg_benchmark/lib/R"}` to its `kernel.json` (path from `jupyter kernelspec list`) so R is found.

BANKSY (`scripts/seg_postprocessing/banksy_*.py`):
```
mamba install -n analysis -c conda-forge -c bioconda bioconductor-banksy=1.6.0 bioconductor-spatialexperiment
```

### Segmentation env
Install into `segmentation`. Cellpose, Cellpose-SAM, CellSAM, Baysor, ComSeg and Proseg run through [Sopa](https://github.com/prism-oncology/sopa) ([setup](https://prism-oncology.github.io/sopa/getting_started/)).
+ **GPU:** expose NVIDIA devices.
  ```bash
  ROOTFS=/raid/enroot/data/user-$(id -u)/<container_name>
  printf 'NVIDIA_VISIBLE_DEVICES=all\nNVIDIA_DRIVER_CAPABILITIES=compute,utility\n' >> $ROOTFS/etc/environment
  ```
+ **Cellpose, ComSeg:** `pip install cellpose==3.1.1.1 comseg==1.8.5`
+ **Proseg:** install cargo (`curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh`), then `cargo install proseg`.
+ **Visium:** `pip install rasterio`
+ **CellSAM:** `--no-deps` keeps the CUDA PyTorch build.
  ```bash
  mamba activate segmentation
  python -c "import torch; assert torch.__version__=='2.10.0+cu128', torch.__version__"
  pip install --no-deps --index-url https://download.pytorch.org/whl/cu128 torchvision==0.25.0+cu128
  pip install --no-deps git+https://github.com/facebookresearch/segment-anything.git
  pip install --no-deps git+https://github.com/vanvalenlab/cellSAM.git
  pip install kornia
  python -c "import torchvision, kornia, segment_anything, cellSAM; print('ok')"
  ```
  Download weights before running on compute nodes, using a DeepCell token ([how to get one](https://github.com/vanvalenlab/cellSAM)) saved at the path below:
  ```bash
  export DEEPCELL_ACCESS_TOKEN=$(tr -d '[:space:]' < /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/deepcell_token.txt)
  export HOME=/home/ubuntu
  python -c "from cellSAM import get_model; get_model()"
  ```
  Don't pin `tifffile`; use `tifffile.memmap` instead of `aszarr()` if needed.
+ **Cellpose-SAM:** venv on top of `segmentation`, so `segmentation` keeps Cellpose 3.
  ```bash
  mamba activate segmentation
  python -m venv --system-site-packages /opt/cellposesam
  source /opt/cellposesam/bin/activate
  pip install cellpose==4.2.1.1
  python -c "import torch, sopa, cellpose; print(torch.__version__, sopa.__version__, cellpose.version, cellpose.__file__)"  # 2.10.0+cu128 2.1.10 4.2.1.1 /opt/cellposesam/...
  mkdir -p /home/ubuntu/.cellpose/models
  CELLPOSE_LOCAL_MODELS_PATH=/home/ubuntu/.cellpose/models python -c "from cellpose import models; models.CellposeModel(pretrained_model='cpsam_v2')"
  chmod -R a+rX /opt/cellposesam /home/ubuntu/.cellpose
  deactivate && python -c "import cellpose; print(cellpose.version)"  # 3.1.1.1
  ```
+ **Spots-in-Space:** separate env `SIS_segmentation`; it needs Cellpose 3.0.7 and other scipy/pyarrow versions.

### Export
Jobs use the `.sqsh`, so export after changes. Clean caches inside the container first (`mamba clean -a -y; pip cache purge`). To edit an existing image, start with `enroot create --name <container_name> <image>.sqsh`.
```
enroot export -f -o <image>.sqsh <container_name>
```
