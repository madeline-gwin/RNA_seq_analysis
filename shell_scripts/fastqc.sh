#!/bin/bash

module load fastqc
module load multiqc

Input_dir="/scratch/magwin/callus_RNA/pre-processing/after_filtering/"
Output_dir="/scratch/magwin/callus_RNA/pre-processing/after_filtering/QC"

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

