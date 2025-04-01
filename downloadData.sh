#!/bin/bash 

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=download-Data
#SBATCH --time=05:00:00   # HH/MM/SS
#SBATCH --mem=64G         # memory requested, units available : K,M,G,T

echo "Starting at:" `date`
echo "This is job #:" $SLURM_JOB_ID
echo "Running on node:" `hostname`
echo "Running on cluster:" $SLURM_CLUSTER_NAME

csvDir=".."
cd /athena/angsd/scratch/few4007/angsd-files/project_2/RawFastqZipped

for SRR in $(grep -o 'SRR[0-9]\+' $csvDir/SraRunTable.csv); do
    echo "Processing $SRR"
    prefetch $SRR
    fastq-dump --split-files "$SRR/$SRR.sra"
done