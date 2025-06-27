#!/bin/bash
#SBATCH --account=MATH030465
#SBATCH --partition=short
#SBATCH --mail-type=ALL
#SBATCH --time=00:06:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=3000MB

module load languages/R/4.3.3

R CMD BATCH mf.r