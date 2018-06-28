#!/bin/bash

# run this script as follows:
#sbatch --array=1-20 latentspace.sbatch <args-for-R>

#SBATCH --job-name latentspace
#SBATCH --array=1-20
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=shas
#SBATCH --time=16:00:00
#SBATCH --workdir=/projects/imtaylor@colostate.edu/latentspace/
#SBATCH -o /home/imtaylor@colostate.edu/logs/latentspace_%A_%a.out
#SBATCH -e /home/imtaylor@colostate.edu/logs/latentspace_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ian.taylor@colostate.edu

module purge
module load R

Rscript /home/imtaylor@colostate.edu/independent_covariate.R \
    jobname nodcor$SLURM_ARRAY_TASK_ID \
    n 20 \
    run $SLURM_ARRAY_TASK_ID \
    workers 1 \
    reps 200
