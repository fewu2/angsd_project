#!/bin/bash 

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=multiqc
#SBATCH --time=36:00:00   # HH/MM/SS
#SBATCH --mem=64G         # memory requested, units available : K,M,G,T

cd /athena/angsd/scratch/few4007/angsd-files/project_2

# multiqc trim_out -o multiqc_out
multiqc . -o . --interactive
