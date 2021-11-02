#!/bin/bash

#SBATCH --job-name latentspace_projections
#SBATCH --array=0-3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=shas
#SBATCH --time=12:00:00
#SBATCH -o /home/imtaylor@colostate.edu/logs/latentspace_%A_%a.out
#SBATCH -e /home/imtaylor@colostate.edu/logs/latentspace_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ian.taylor@colostate.edu

module purge
module load R

Rscript simulation-run.R $SLURM_ARRAY_TASK_ID
