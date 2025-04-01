#!/bin/bash 

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=make_index
#SBATCH --time=05:00:00   # HH/MM/SS
#SBATCH --mem=64G         # memory requested, units available : K,M,G,T

echo "Starting at: $(date)"
echo "This is job #: $SLURM_JOB_ID"
echo "Running on node: $(hostname)"
echo "Running on cluster: $SLURM_CLUSTER_NAME"

assembly="GRCh38.primary_assembly.genome.fa"
gtf="gencode.v47.annotation.gtf"
index_out="HS_index"

# Set working directory 
cd /athena/angsd/scratch/few4007/angsd-files/project_2/

# Check if genome and annotation files exist
if [ ! -f "$assembly" ]; then
    echo "Specified assembly input is missing. Aborted."
    exit 1
fi

if [ ! -f "$gtf" ]; then
    echo "Specified .GTF annotation file is missing. Aborted."
    exit 1
fi

# Check if output directory exists, if not, create it
if [ ! -d "$index_out" ]; then
    mkdir -p "$index_out"
    STAR --runMode genomeGenerate \
     --genomeDir "$index_out" \
     --genomeFastaFiles "$assembly" \
     --sjdbGTFfile "$gtf" \
     --runThreadN 16
else
    echo "Input index already exists! Please delete it if you want to make a new one."
fi
