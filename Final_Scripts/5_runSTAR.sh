#!/bin/bash 

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=runSTAR
#SBATCH --time=36:00:00   # HH/MM/SS
#SBATCH --mem=64G         # memory requested, units available : K,M,G,T

# Set inputs
threads=16
indir="../HS_index"
readsdir="../trim_out"
outdir="align_out_star"

cd /athena/angsd/scratch/few4007/angsd-files/project_2/$outdir
# Loop through all sample subdirectories in expected path format of trim_out
for sampdir in ${readsdir}/*/; do
    # Get files in expected naming format/convention
    fwd=$(ls "${sampdir}"*_1_val_1.fq.gz 2>/dev/null)
    bkwd=$(ls "${sampdir}"*_2_val_2.fq.gz 2>/dev/null)

    # Make sure both files are present
    if [[ -f "$fwd" && -f "$bkwd" ]]; then
        # Extract the base name of sample from the filename
        sample=$(basename "$fwd" | sed 's/_1_val_1\.fq\.gz//')
        echo "Running STAR on ${sample} reads..."

        # Run STAR
        STAR --runMode alignReads \
            --genomeDir "$indir" \
            --readFilesIn "$fwd" "$bkwd" \
            --readFilesCommand zcat \
            --outFileNamePrefix "$sample" \
            --outSAMtype BAM SortedByCoordinate \
            --runThreadN "$threads"
    else
        echo "Missing files in $sampdir, skipping."
    fi
done
