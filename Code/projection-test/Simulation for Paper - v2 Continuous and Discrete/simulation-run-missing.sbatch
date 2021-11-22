#!/bin/bash

#SBATCH --job-name latentspace_projections
#SBATCH --array=1100-1110,1112-1114,1116-1119,1340-1354,1356-1359,1573,1580-1599
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --partition=shas
#SBATCH --time=24:00:00
#SBATCH -o /home/imtaylor@colostate.edu/logs/latentspace_%A_%a.out
#SBATCH -e /home/imtaylor@colostate.edu/logs/latentspace_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ian.taylor@colostate.edu

module purge
module load R

Rscript simulation-run.R $SLURM_ARRAY_TASK_ID
