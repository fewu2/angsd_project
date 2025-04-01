#!/bin/bash 

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=unzip
#SBATCH --time=05:00:00   # HH/MM/SS
#SBATCH --mem=30G         # memory requested, units available : K,M,G,T

input_dir="/athena/angsd/scratch/few4007/angsd-files/project_2/RawFastqZipped"
cd $input_dir

for file in $input_dir/*_1.fastq; do
  base=$(basename "$file" _1.fastq)
  
  # Check if a compressed version already exists for both fastqs
  # If not, compress.
  if [ ! -f "$input_dir/${base}_1.fastq.gz" ]; then
    echo "Compressing $input_dir/${base}_1.fastq to $input_dir/${base}_1.fastq.gz"
    gzip "$input_dir/${base}_1.fastq"  
  fi

  if [ ! -f "$input_dir/${base}_2.fastq.gz" ]; then
    echo "Compressing $input_dir/${base}_2.fastq to $input_dir/${base}_2.fastq.gz"
    gzip "$input_dir/${base}_2.fastq"  
  fi
done