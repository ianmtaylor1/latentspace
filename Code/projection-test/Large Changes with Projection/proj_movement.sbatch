#!/bin/bash

#SBATCH --job-name latentspace
#SBATCH --array=1-30
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=shas
#SBATCH --time=16:00:00
#SBATCH --chdir=/projects/imtaylor@colostate.edu/latentspace/proj_movement/
#SBATCH -o /home/imtaylor@colostate.edu/logs/latentspace_%A_%a.out
#SBATCH -e /home/imtaylor@colostate.edu/logs/latentspace_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ian.taylor@colostate.edu

module purge
module load R

Rscript /projects/imtaylor@colostate.edu/latentspace/proj_movement/measure-change.R \
	run $SLURM_ARRAY_TASK_ID \
	outfilebase proj_movement

