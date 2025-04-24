#!/bin/bash 

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=runQoRTS
#SBATCH --time=44:00:00   # HH/MM/SS
#SBATCH --mem=48G         # memory requested, units available : K,M,G,T

# Note we need to switch into the qorts conda environment
cd "/athena/angsd/scratch/few4007/angsd-files/project_2" 

indir="align_out_star"  
# original run for qorts_out did not use --stranded parameter
outdir="qorts_out_2"  
inGTF="gencode.v47.annotation.gtf"

# Loop over all BAM files in the directory and run QoRTs
for bam in ${indir}/*Aligned.sortedByCoord.out.bam; do
    # Extract the sample name by removing the suffix from the filename
    sample=$(basename "$bam" Aligned.sortedByCoord.out.bam)

    # Define output directory for this sample
    sample_outdir="$outdir/$sample"

    # Only run if the output directory does not already exist
    if [ ! -d "$sample_outdir" ]; then
        echo "Running QoRTs for $sample..."
        java -jar /athena/angsd/scratch/mef3005/share/envs/qorts/share/qorts-1.3.6-1/QoRTs.jar QC \
            --stranded \
            --numThreads 24 \
            --maxReadLength 76 \
            "$bam" \
            "$inGTF" \
            "$sample_outdir"
    else
        echo "Skipping $sample (output already exists)."
    fi
done
