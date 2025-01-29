#! /bin/bash
# Helper script to start jupyter notebook on CPU
#SBATCH -o start_notebook_%j.out
#SBATCH -e start_notebook_%j.out
#SBATCH -J jupyter
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100Gb
#SBATCH -t 12:00:00

# activate conda
# NOTE: conda activate svd does not work. Other workaround would be to source $CONDA_ROOT/etc/profile.d/conda.sh
source activate cellseg

# run script
BASE=/home/hspitzer/projects/cellseg_benchmark
# prevent other people from reading the log file
chmod og-r start_notebook_${SLURM_JOB_ID}.out
jupyter lab --no-browser --ip=0.0.0.0 --port=8081 --notebook-dir $BASE
