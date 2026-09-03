#!/bin/bash

#SBATCH --job-name=fastp_take1          # job name
#SBATCH --partition=capitulab           # partition
#SBATCH --ntasks=1                      # number of tasks across all nodes
#SBATCH --cpus-per-task=4               # number of cpus per task
#SBATCH --mem=12GB                      # total memory requested
#SBATCH --time=72:00:00                 # Run time (D-HH:MM:SS)
#SBATCH --output=job-%j.out             # Output file. %j is replaced with job ID
#SBATCH --error=job-%j.err              # Error file. %j is replaced with job ID
#SBATCH --mail-type=ALL                 # will send email for begin,end,fail
#SBATCH --mail-user=magwin@clemson.edu


# Set input and output base directories
INPUT_DIR="/scratch/magwin/BCY_floret_dev/raw_reads/"
OUTPUT_DIR="/scratch/magwin/BCY_floret_dev/pre-processing/after_filtering/"
THREADS=4

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Find all R1 files recursively
find "$INPUT_DIR" -name "*_1.fq.gz" | while read R1; do
    # Derive corresponding R2
    R2="${R1/_1.fq.gz/_2.fq.gz}"

    # Extract sample name (removes path and suffix)
    SAMPLE=$(basename "$R1" _1.fq.gz)

    # Create output filenames
    OUT_R1="$OUTPUT_DIR/${SAMPLE}_1_trimmed.fq.gz"
    OUT_R2="$OUTPUT_DIR/${SAMPLE}_2_trimmed.fq.gz"
    HTML="$OUTPUT_DIR/${SAMPLE}_fastp.html"
    JSON="$OUTPUT_DIR/${SAMPLE}_fastp.json"

    echo "Processing sample: $SAMPLE"
    fastp \
        -i "$R1" -I "$R2" \
        -o "$OUT_R1" -O "$OUT_R2" \
        --detect_adapter_for_pe \
        -h "$HTML" -j "$JSON" \
        -w "$THREADS"
done
