#!/bin/bash

#SBATCH --job-name=mort
#SBATCH --partition=long
#SBATCH --array=1-5
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --output=output/reduced_chronological_kfold_%a.out
#SBATCH --error=output/reduced_chronological_kfold_%a.err


# Check what versions of anaconda are available with ‘module avail anaconda’

module load miniconda3
conda activate
conda activate allpip

python reduced_dataset_chronological.py $SLURM_ARRAY_TASK_ID 
