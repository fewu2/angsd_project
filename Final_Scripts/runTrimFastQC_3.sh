#!/bin/bash 

#SBATCH --partition=angsd_class
#SBATCH --nodes=4
#SBATCH --ntasks=1
#SBATCH --job-name=fastqc-trim
#SBATCH --time=36:00:00   # HH/MM/SS
#SBATCH --mem=64G         # memory requested, units available : K,M,G,T

cd "/athena/angsd/scratch/few4007/angsd-files/project_2"
input_dir="RawFastqZipped"
output_dir="trim_out"

# Loop through all .fastq files in the directory
# To avoid running trim_galore twice on each sample, we loop only through the first fastq file of each sample
# Then use basename to pull the second fastq for that sample
# Before finally submitting a paired-end input to trim-galore. 
for file in "$input_dir"/*_1.fastq.gz; do
  base=$(basename "$file" _1.fastq.gz)
  
  # Check if both compressed forward and reverse fastq files exist
  if [ -f "$input_dir/${base}_1.fastq.gz" ] && [ -f "$input_dir/${base}_2.fastq.gz" ] && [ ! -d "$output_dir/$base" ]; then
    # Create output directories for trimmed and fastqc results
    mkdir -p "$output_dir/$base/fastqc_out"
    
    echo "Running trim_galore and fastqc on $base..."
    
    # Run trim_galore with --fastqc option, this will run fastqc alongside trimming
    trim_galore --paired \
      --fastqc \
      --illumina \
      -j 4 \
      --fastqc_args "--extract --outdir=${output_dir}/${base}/fastqc_out" \
      --output_dir "$output_dir/$base" \
      "$input_dir/${base}_1.fastq.gz" \
      "$input_dir/${base}_2.fastq.gz"
    
    echo "Trimming and FastQC completed for $base."
  elif [ ! -d "$output_dir/$base/fastqc_out" ]; then
    if [ ! -f "$input_dir/${base}_2.fastq.gz" ]; then
      echo "${base} is missing the second strand!"
    fi
  fi


done
