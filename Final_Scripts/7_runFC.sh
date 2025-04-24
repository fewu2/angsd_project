#!/bin/bash 

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=FC
#SBATCH --time=36:00:00   # HH/MM/SS
#SBATCH --mem=32G         # memory requested, units available : K,M,G,T

cd "/athena/angsd/scratch/few4007/angsd-files/project_2" 

# Path to the annotation GTF file (adjust as necessary)
inGTF="gencode.v47.annotation.gtf"
outfile="v4AllCountsRNA.txt"
indir="align_out_star"

# Run featureCounts with the appropriate settings
# featureCounts -T 16 -a $inGTF -o $outfile -t exon -g gene_id -p --countReadPairs -s 2 -Q 10 ${indir}/*Aligned.sortedByCoord.out.bam
# featureCounts -T 16 -a $inGTF -o $outfile -t exon -g gene_id -p -s 2 -Q 10 ${indir}/*Aligned.sortedByCoord.out.bam
featureCounts -T 16 -a $inGTF -o $outfile -t gene -g gene_id -p --countReadPairs -s 2 -Q 10 ${indir}/*Aligned.sortedByCoord.out.bam
