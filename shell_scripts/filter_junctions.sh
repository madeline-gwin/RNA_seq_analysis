#!/bin/bash

#SBATCH --job-name=filter_junctions_Hap2_GTF    # job name
#SBATCH --partition=capitulab                   # partition
#SBATCH --ntasks=1                              # number of tasks across all nodes
#SBATCH --cpus-per-task=1                       # number of cpus per task
#SBATCH --mem=10GB                              # total memory requested
#SBATCH --time=72:00:00                         # Run time (D-HH:MM:SS)
#SBATCH --output=job-%j.out                     # Output file. %j is replaced with job ID
#SBATCH --error=job-%j.err                      # Error file. %j is replaced with job ID
#SBATCH --mail-type=ALL                         # will send email for begin,end,fail
#SBATCH --mail-user=magwin@clemson.edu

set -o pipefail

##############################
# HARD-CODED VARIABLES
##############################

# Path to directory containing SJ.out.tab files and where outputs will be saved
#JUNCTIONDIR="/scratch/magwin/BCY_floret_dev/genomes/Hap1/genome_indices/Junctions"
JUNCTIONDIR="/scratch/magwin/BCY_floret_dev/genomes/Hap2/genome_indices/Junctions"

# Prefix for the final filtered junction list filename
#SJ_LISTNAME="BCY_Hap1"
SJ_LISTNAME="BCY_Hap2"

# Filter settings:
SCAFFOLD_STRING="NA"         # String indicating scaffold sequences, "NA" disables scaffold filtering
REMOVE_NC_JUNC="yes"         # "yes" to remove non-canonical junctions, "no" to keep them
UNIQUE_NUM=2                 # Minimum number of uniquely mapping reads to support junction

##############################
# SCRIPT LOGIC STARTS HERE
##############################

# Count total junctions (unique by first 6 columns)
TOTJUNC=$(cat ${JUNCTIONDIR}/*SJ.out.tab | cut -f1-6 | sort | uniq | wc -l)

# Count annotated junctions (column 6 == 1)
ANNJUNC=$(cat ${JUNCTIONDIR}/*SJ.out.tab | awk '($6==1)' | cut -f1-6 | sort | uniq | wc -l)

# Count non-canonical junctions (column 5 == 0)
NCJUNC=$(cat ${JUNCTIONDIR}/*SJ.out.tab | awk '($5==0)' | cut -f1-6 | sort | uniq | wc -l)

echo "The total number of junctions was ${TOTJUNC}, ${ANNJUNC} were already annotated"

if [[ "$SCAFFOLD_STRING" == "NA" ]]; then
    echo "No scaffold sequence specified"
    if [[ "$REMOVE_NC_JUNC" == "yes" ]]; then
        echo "Removing ${NCJUNC} non-canonical junctions"
        cat ${JUNCTIONDIR}/*SJ.out.tab | \
        awk '($5 > 0 && $6==0)' | \
        awk -v num="$UNIQUE_NUM" '($7>=num || ++a[$1,$2,$3,$4,$5,$6]==num)' | \
        cut -f1-6 | sort | uniq > ${JUNCTIONDIR}/${SJ_LISTNAME}_SJ.filtered.tab
    else
        echo "Not removing ${NCJUNC} non-canonical junctions"
        cat ${JUNCTIONDIR}/*SJ.out.tab | \
        awk '($6==0)' | \
        awk -v num="$UNIQUE_NUM" '($7>=num || ++a[$1,$2,$3,$4,$5,$6]==num)' | \
        cut -f1-6 | sort | uniq > ${JUNCTIONDIR}/${SJ_LISTNAME}_SJ.filtered.tab
    fi
else
    SCAFFOLDJUNC=$(cat ${JUNCTIONDIR}/*SJ.out.tab | awk -v var="$SCAFFOLD_STRING" '($1 ~ var)' | cut -f1-6 | sort | uniq | wc -l)
    echo "Removing ${SCAFFOLDJUNC} junctions from scaffold sequence"
    if [[ "$REMOVE_NC_JUNC" == "yes" ]]; then
        echo "Removing ${NCJUNC} non-canonical junctions"
        cat ${JUNCTIONDIR}/*SJ.out.tab | \
        awk -v var="$SCAFFOLD_STRING" '($1 !~ var && $5 > 0 && $6==0)' | \
        awk -v num="$UNIQUE_NUM" '($7>=num || ++a[$1,$2,$3,$4,$5,$6]==num)' | \
        cut -f1-6 | sort | uniq > ${JUNCTIONDIR}/${SJ_LISTNAME}_SJ.filtered.tab
    else
        echo "Not removing ${NCJUNC} non-canonical junctions"
        cat ${JUNCTIONDIR}/*SJ.out.tab | \
        awk -v var="$SCAFFOLD_STRING" '($1 !~ var && $6==0)' | \
        awk -v num="$UNIQUE_NUM" '($7>=num || ++a[$1,$2,$3,$4,$5,$6]==num)' | \
        cut -f1-6 | sort | uniq > ${JUNCTIONDIR}/${SJ_LISTNAME}_SJ.filtered.tab
    fi
fi

FINALNUM=$(wc -l < ${JUNCTIONDIR}/${SJ_LISTNAME}_SJ.filtered.tab)
echo "After removing junctions supported by less than ${UNIQUE_NUM} uniquely mapped reads,"
echo "the final filtered list contains ${FINALNUM} junctions and can be used for 2nd-pass read mapping"
echo "This list can be found at: ${JUNCTIONDIR}/${SJ_LISTNAME}_SJ.filtered.tab"
