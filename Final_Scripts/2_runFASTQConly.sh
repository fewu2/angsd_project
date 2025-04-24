#!/bin/bash 

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=fastqc
#SBATCH --time=24:00:00   # HH/MM/SS
#SBATCH --mem=32G         # memory requested, units available : K,M,G,T

cd "/athena/angsd/scratch/few4007/angsd-files/project_2"
input_dir="RawFastqZipped"
output_dir="fastqcOnly_out"

# Loop through all .fastq.gz files in the directory
# Data is paired-end, and both strands are needed as inputs
# To avoid parsing each pair twice, we loop only through the first fastq file of each sample
# Then use basename to pull the second fastq for that sample
# This is a naive approach. Please be aware that for a given sample, if the FIRST strand is missing but not the second...
# ... this script will miss the sample entirely and fail to throw a warning. If the second strand is missing,
# ... this script will throw a warning and skip the sample.
for file in "$input_dir"/*_1.fastq.gz; do
  base=$(basename "$file" _1.fastq.gz)
  
  # Check if both compressed forward and reverse fastq files exist, as well as output directory
  if [ -f "$input_dir/${base}_1.fastq.gz" ] && [ -f "$input_dir/${base}_2.fastq.gz" ] && [ ! -d "$output_dir/$base" ]; then
    # Create output directories for trimmed and fastqc results
    mkdir -p "$output_dir/$base"
    
    echo "Running fastqc on $base..."
    
    # Run fastqc only
    fastqc  --extract \
        -o "${output_dir}/${base}" \
        "$input_dir/${base}_1.fastq.gz" \
        "$input_dir/${base}_2.fastq.gz"
    
    echo "FastQC completed for $base."

  # This script is desisnged to be run multiple times if full job is not completed in one run
  # Thus, check if an output directory already exists for the given sample
  # If so, do not overwrite the output directory, but throw a warning
  elif [ ! -d "${output_dir}/${base}/fastqc_out" ]; then

    # If output directory does not exist, the issue will be that the second strand is missing
    # Throw a warning and skip the sample
    if [ ! -f "${input_dir}/${base}_2.fastq.gz" ]; then
      echo "${base} is missing the second strand!"
    fi
  else
    echo "${base} already has an output directory. Please delete it if you want new results! Skipping..."
  fi
done
