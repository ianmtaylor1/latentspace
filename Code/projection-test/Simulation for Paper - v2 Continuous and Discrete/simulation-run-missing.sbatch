#!/bin/bash

#SBATCH --job-name latentspace_projections
#SBATCH --array=123,140-159,380-399,408,600,603-606,609-610,613,617,620-639,860-879,1100-1102,1104-1119,1328,1340-1344,1346-1359,1580-1581,1583-1599
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=shas
#SBATCH --time=24:00:00
#SBATCH -o /home/imtaylor@colostate.edu/logs/latentspace_%A_%a.out
#SBATCH -e /home/imtaylor@colostate.edu/logs/latentspace_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ian.taylor@colostate.edu

module purge
module load R

Rscript simulation-run.R $SLURM_ARRAY_TASK_ID
