#!/bin/bash

#SBATCH --nodes=2
#SBATCH --ntasks-per-node=12
#SBATCH --partition=shas
#SBATCH --time=24:00:00
#SBATCH --workdir=~/projects

module purge
module load R
module load intel
module load impi

srun -n 1 R --slave -f ~/independent_covariate.R
