#!/bin/bash 

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=checkStranded
#SBATCH --time=12:00:00   # HH/MM/SS
#SBATCH --mem=12G         # memory requested, units available : K,M,G,T

# This script requires the rqsec conda env
# This is not part of the script, but as a note, I used the following code to create the .bed version of the gtf
# conda activate angsd
# gtf2bed < "gencode.v47.annotation.gtf" > "gencode.v47.annotation.bed"   


cd "/athena/angsd/scratch/few4007/angsd-files/project_2"

bedfile="gencode.v47.annotation.bed"

for bam in align_out_star/SRR*Aligned.sortedByCoord.out.bam; do
  echo "Running infer_experiment.py on $bam..."
  infer_experiment.py -r "$bedfile" -i "$bam"
done
