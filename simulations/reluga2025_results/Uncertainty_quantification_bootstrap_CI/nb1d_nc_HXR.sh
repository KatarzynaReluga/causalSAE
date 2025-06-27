#!/bin/bash
#SBATCH --account=MATH030465
#SBATCH --partition=short
#SBATCH --mail-type=ALL
#SBATCH --time=07:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=3000MB

module load languages/R/4.4.2

R CMD BATCH nb1d_nc_HXR.r