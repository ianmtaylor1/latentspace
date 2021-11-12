#!/bin/bash

#SBATCH --job-name latentspace_projections
#SBATCH --array=100,120-159,163-167,238,241,252,278,360-399,600-639,840-879,1000-1679
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
