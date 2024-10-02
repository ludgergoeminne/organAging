#!/bin/bash

#SBATCH --job-name=aft
#SBATCH --partition=bigmem
#SBATCH --array=1-5
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --output=output/aft_kfold_%a.out
#SBATCH --error=output/aft_kfold_%a.err


# Check what versions of anaconda are available with ‘module avail anaconda’

module load miniconda3
conda activate
conda activate allpip

python aft_clocks.py $SLURM_ARRAY_TASK_ID 
