#!/bin/bash

#SBATCH --job-name latentspace_projections
#SBATCH --array=1000-1679,147,152,157,382,396-398,631,871-879
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=shas
#SBATCH --time=24:00:00
#SBATCH -o /home/imtaylor@colostate.edu/logs/latentspace_%A_%a.out
#SBATCH -e /home/imtaylor@colostate.edu/logs/latentspace_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ian.taylor@colostate.edu

module purge
module load R

Rscript simulation-run.R $SLURM_ARRAY_TASK_ID
