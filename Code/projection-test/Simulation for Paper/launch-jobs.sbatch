#!/bin/bash

#SBATCH --job-name latentspace_projections
#SBATCH --array=1-320
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=shas
#SBATCH --time=12:00:00
#SBATCH --chdir=/projects/imtaylor@colostate.edu/latentspace/proj_for_paper/
#SBATCH -o /home/imtaylor@colostate.edu/logs/latentspace_%A_%a.out
#SBATCH -e /home/imtaylor@colostate.edu/logs/latentspace_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ian.taylor@colostate.edu

module purge
module load R

Rscript /projects/imtaylor@colostate.edu/latentspace/proj_for_paper/capturing-delta.R \
	run $SLURM_ARRAY_TASK_ID \
	outfilebase projections_for_paper

