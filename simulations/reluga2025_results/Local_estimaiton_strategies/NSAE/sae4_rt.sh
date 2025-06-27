#!/bin/bash
#SBATCH --account=MATH030465
#SBATCH --partition=short
#SBATCH --mail-type=ALL
#SBATCH --time=02:30:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1500MB

module load lang/r/4.3.0-bioconductor-gcc

srun R CMD BATCH sae4_rt.r