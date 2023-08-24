#!/bin/bash

#SBATCH --job-name latentspace_projections
#SBATCH --array=1944,2034,2136,5454,5662
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=amilan
#SBATCH --time=12:00:00
#SBATCH -o /home/imtaylor@colostate.edu/logs/latentspace_%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ian.taylor@colostate.edu

module purge
module load anaconda
conda activate netreg

Rscript simulation-run.R $SLURM_ARRAY_TASK_ID
