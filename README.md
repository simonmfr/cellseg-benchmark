# Cell segmentation benchmark

## Structure of this repository
- `cellseg_benchmark`
  Function definitions (the package)
- `scripts`
  Scripts for metric calculations and segmentation algorithms, incl. `scripts/sbatch_utils` for generating per-method sbatch scripts
- `notebooks`
  Jupyter notebooks for development and analysis
- `configs`
  Configuration files (e.g. VPT segmentation configs)
- `archive`
  Symlink to raw MERSCOPE data on DSS
- `data`
  Symlink to processed data on DSS

## Development
We are using `ruff` and the [ruff pre-commit hook](https://github.com/astral-sh/ruff-pre-commit) to check and format the code and docstrings

**Installation**: Install ruff and pre-commit in your environment and install the pre-commit hooks for ruff defined in `.pre-commit-config.yml`
```
pip install ruff
pip install pre-commit
pre-commit install
# (optional: run against all files & fix any errors that are in your current codebase)
pre-commit run --all-files
```

**Basic usage**: The ruff config is located in `pyproject.toml`. See the [ruff documentation of rules](https://docs.astral.sh/ruff/rules/) for all possible rules that we can enable / disable.
As we have installed the pre-commit hook, ruff formatting and liniting will run automatically for all changed files whenever you do git commit. If there are errors, you will get a detailled messaged of the offending code and the error. Fix the errors, add the changed file and try to commit again.

You can also manually run the ruff formatter and checker on all files with: 
```
ruff format
ruff check --fix
```
or
```
pre-commit run --all-files
```

If necessary, you can also temporarily disable all pre-commit hooks when committing by using the `--no-verify` flag with `git commit`.

## Container setup (enroot LRZ)

First, import the ubuntu docker container and setup python and mamba:

```
enroot import docker://ubuntu
enroot create --name <container_name> ubuntu.sqsh
enroot start -r <container_name>

# --- inside container ---
apt-get update
apt-get install vim
apt-get install curl
apt-get install git
apt-get install nano

echo "HOME=/home/ubuntu" >> /etc/profile
exit

# --- restart container without root access ---
enroot start <container_name>

cd tmp
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh

cat > /home/ubuntu/activate_env.sh <<'EOL'
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/ubuntu/miniforge3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/ubuntu/miniforge3/etc/profile.d/conda.sh" ]; then
        . "/home/ubuntu/miniforge3/etc/profile.d/conda.sh"
    else
        export PATH="/home/ubuntu/miniforge3/bin:$PATH"
    fi
fi
unset __conda_setup

if [ -f "/home/ubuntu/miniforge3/etc/profile.d/mamba.sh" ]; then
    . "/home/ubuntu/miniforge3/etc/profile.d/mamba.sh"
fi
# <<< conda initialize <<<
EOL

# set $BASH_ENV activate_env.sh
echo "BASH_ENV=/home/ubuntu/activate_env.sh" >> /etc/environment

# re-source .bashrc -> this should activate conda base environment
source ~/.bashrc
exit
```

**Container configuration:**
First, check that your project folders, e.g. data and output folders, are available to the container. To do this, locate `/etc/fstab`, open the file and add `<outside_path> <inside_path> none x-create=dir,bind`. For further information for this consult https://github.com/NVIDIA/enroot/blob/main/doc/image-format.md.
As we want to use the cellseg_benchmark code, create a gitrepos directory and clone the cellseg_benchmark repository
```
mkdir ~/gitrepos
cd ~/gitrepos
git clone  https://github.com/simonmfr/cellseg-benchmark.git
```
To auto-pull the repository to the lastest version each time the container is started, exit now the container and restart it with root access. Inside the container create a temporary file with the automatic update code:
```
cat > /tmp/rc_insert.sh <<'EOF'
# BEGIN cellseg-benchmark auto-update
REPO="/home/ubuntu/gitrepos/cellseg-benchmark" #if you have stored the cellseg_benchmark repository at a different places, change this to the updated path
[ -d "$REPO/.git" ] || REPO="${HOME}/gitrepos/cellseg-benchmark"

if command -v git >/dev/null 2>&1 && [ -d "$REPO/.git" ]; then
  export GIT_TERMINAL_PROMPT=0
  git -C "$REPO" pull --ff-only -q || true
fi
# END cellseg-benchmark auto-update
EOF
```
Insert this snippet into `/etc/rc` and ensure that `/etc/rc` remains executable:
```
if ! grep -q "BEGIN cellseg-benchmark auto-update" /etc/rc; then
  awk 'NR==1 {print; system("cat /tmp/rc_insert.sh"); next} {print}' /etc/rc > /tmp/rc.new \
    && cat /tmp/rc.new > /etc/rc
fi
chmod +x /etc/rc
```

**Environment setup:**
Now exit and restart the container again, this time without root access. To install cellseg_benchmark with dea possibilities and jupyter notebook support, do the following:
```
mamba create -n cellseg_benchmark python=3.12 -c conda-forge -c bioconda r-base rpy2 anndata2ri zlib cmake -y
mamba activate cellseg_benchmark
cd ~/gitrepos/cellseg_benchmark
pip install -e .
pip install setuptools==80.0.0
cd ..
git clone https://github.com/AllenInstitute/cell_type_mapper.git
git clone https://github.com/jonas2612/spatialdata.git
cd cell_type_mapper
pip install .
cd ../spatialdata
pip install .
```
If you want to register the kernel to work with juypter notebooks, do the following:
```
python -m ipykernel install --user --name cellseg_benchmark
nano /home/ubuntu/miniforge3/jupyter/kernels/cellseg_benchmark/kernel.json
```
The last command opens the kernel specs of the jupyter kernel. Add the following:
```
"env": {
  "R_HOME": "/home/ubuntu/miniforge3/envs/cellseg_benchmark/lib/R #or wherever the R version you want to use is located
}
```
This ensures the recognition of the R environment within jupyter.

We're still missing the R packages needed for differential expression testing. For this execute the following:
```R
install.packages(c("dplyr", "data.table", "Matrix"))
install.packages("BiocManager")
BiocManager::install(c(
  "edgeR",
  "variancePartition",
  "limma",
  "BiocParallel",
  "MAST",
  "SummarizedExperiment"
))
```

For BANKSY brain-region mapping (`scripts/seg_postprocessing/banksy_*.py`), install into the `analysis` env:
```
mamba install -n analysis -c conda-forge -c bioconda bioconductor-banksy=1.6.0 bioconductor-spatialexperiment
```

**Segmentation setup:**
We're assuming that the container is properly setup with python.
+ **GPU setup:** Expose NVIDIA devices in the enroot container.
  ```bash
  ROOTFS=/raid/enroot/data/user-$(id -u)/<container_name>
  printf 'NVIDIA_VISIBLE_DEVICES=all\nNVIDIA_DRIVER_CAPABILITIES=compute,utility\n' \
    >> $ROOTFS/etc/environment
  ```
+ **Cellpose and ComSeg installation:** `pip install cellpose==3.1.1.1 comseg==1.8.5` in the `segmentation` env.
+ **Proseg installation:** install cargo (`curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh`), then `cargo install proseg`.
+ **CellSAM installation:**
  ```bash
  mamba activate segmentation

  python -c "import torch; assert torch.__version__=='2.10.0+cu128', torch.__version__"
  pip install --no-deps --index-url https://download.pytorch.org/whl/cu128 torchvision==0.25.0+cu128
  pip install --no-deps git+https://github.com/facebookresearch/segment-anything.git
  pip install --no-deps git+https://github.com/vanvalenlab/cellSAM.git
  pip install kornia

  python -c "import torchvision, kornia, segment_anything, cellSAM; print('ok')"
  ```

  `--no-deps` prevents replacement of the CUDA-enabled PyTorch installation.

  Create a DeepCell token as described in the [CellSAM repository](https://github.com/vanvalenlab/cellSAM), save it at the path below, and download the weights before running on compute nodes:

  ```bash
  export DEEPCELL_ACCESS_TOKEN=$(tr -d '[:space:]' < \
    /dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark/misc/deepcell_token.txt)
  export HOME=/home/ubuntu

  python -c "from cellSAM import get_model; get_model()"
  ```

  Do not pin `tifffile`; use `tifffile.memmap` instead of `aszarr()` if needed.
+ **Cellpose-SAM installation:** venv on top of `segmentation`, sharing its packages while `segmentation` keeps Cellpose 3.
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

Cellpose, Cellpose-SAM, CellSAM, Baysor, ComSeg and Proseg were run through [Sopa](https://github.com/prism-oncology/sopa). See the [Sopa documentation](https://prism-oncology.github.io/sopa/getting_started/) for installation and method setup.

