#!/bin/bash

#SBATCH --job-name=mort
#SBATCH --partition=long
#SBATCH --array=1-5
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --output=output/mortality_kfold_%a.out
#SBATCH --error=output/mortality_kfold_%a.err


# Check what versions of anaconda are available with ‘module avail anaconda’

module load miniconda3
conda activate
conda activate allpip

python full_first_visit_mortality.py $SLURM_ARRAY_TASK_ID 