#!/bin/bash
#SBATCH --partition=shared-cpu
#SBATCH --mail-type=ALL
#SBATCH --time=03:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000

module load foss/2019a R/3.6.0

srun R CMD BATCH scenarioDB25.r