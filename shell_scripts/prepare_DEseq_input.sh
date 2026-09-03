#!/bin/bash

#SBATCH --job-name=prep_DEseq_Hap2      # job name
#SBATCH --partition=capitulab
#SBATCH --ntasks=1                      # number of tasks across all nodes
#SBATCH --cpus-per-task=1               # number of cpus per task
#SBATCH --mem=5GB                       # total memory requested
#SBATCH --time=72:00:00                 # Run time (D-HH:MM:SS)
#SBATCH --output=job-%j.out             # Output file. %j is replaced with job ID
#SBATCH --error=job-%j.err              # Error file. %j is replaced with job ID
#SBATCH --mail-type=ALL                 # will send email for begin,end,fail
#SBATCH --mail-user=magwin@clemson.edu

# This script prepares your input files so that they can be easily read into DESeq
# make sure to look at which columns are being extracted!
INPUT_DIR="/scratch/magwin/BCY_floret_dev_analysis/genomes/Hap1/read_mapping/"

# Find all R1 files recursively
find "$INPUT_DIR" -name "*ReadsPerGene.out.tab" | while read R1; do

       # Extract sample name (removes path and suffix)
       SAMPLE=$(basename "$R1" ReadsPerGene.out.tab)

       # Create output filenames
       OUT_FILE="${INPUT_DIR}/${SAMPLE}ReadsPerGene_trimmed.out.tab"

       #Extract columns 1 and 2, skipping the 4 line STAR header
       ## Use column 2 for unstranded reads, and column 4 for stranded reads!!##
       awk -F '\t' 'BEGIN {OFS="\t"} NR>4 { print $1, $2 }' "$R1" > "$OUT_FILE"

done
