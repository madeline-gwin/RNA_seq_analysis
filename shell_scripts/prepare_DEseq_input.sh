#!/bin/bash

#This script prepares your input files so that they can be easily read into R

SPECIES="BCY_hap2"

for fn in /scratch/magwin/genome/read_mapping/$SPECIES/*.tab; do
    awk -F '\t' 'NR>4 { print $1, $4 }' "$fn" >tmp_file 
    mv tmp_file "$fn"
done    
