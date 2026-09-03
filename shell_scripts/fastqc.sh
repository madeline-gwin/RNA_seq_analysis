#!/bin/bash

#SBATCH --job-name=fastqc_take2         # job name
#SBATCH --partition=capitulab           # partition
#SBATCH --ntasks=1                      # number of tasks across all nodes
#SBATCH --cpus-per-task=1               # number of cpus per task
#SBATCH --mem=5GB                       # total memory requested
#SBATCH --time=72:00:00                 # Run time (D-HH:MM:SS)
#SBATCH --output=job-%j.out             # Output file. %j is replaced with job ID
#SBATCH --error=job-%j.err              # Error file. %j is replaced with job ID
#SBATCH --mail-type=ALL                 # will send email for begin,end,fail
#SBATCH --mail-user=magwin@clemson.edu

module load biocontainers
module load fastqc
module load multiqc

# Take 1 (before filtering/trimming)
#Input_dir="/scratch/magwin/BCY_floret_dev/raw_reads/"
#Output_dir="/scratch/magwin/BCY_floret_dev/pre-processing/before_filtering/"

# Take 2 (after filtering/trimming)
Input_dir="/scratch/magwin/BCY_floret_dev/pre-processing/after_filtering/"
Output_dir="/scratch/magwin/BCY_floret_dev/pre-processing/after_filtering/QC/"


if [[ -d "$Input_dir" ]]; then
    for f in $(find "$Input_dir" -name "*.fq.gz"); do
        if [[ -f "$f" ]]; then
            echo "Checking quality of $f"
            fastqc -o "$Output_dir" -t 6 "$f"
        else
            echo "$f is not a valid file"
        fi
    done
    echo "Running MultiQC on results..."
    multiqc "$Output_dir" -o "$Output_dir"
else
    echo "$Input_dir is not a valid directory"
fi
