#!/bin/bash
set -o pipefail

##############################
# HARD-CODED VARIABLES
##############################

# Define your species
SPECIES="Marigold"

# Path to directory containing SJ.out.tab files and where outputs will be saved
JUNCTIONDIR="/scratch/magwin/genome/genome_indices/$SPECIES/Junctions"

# Prefix for the final filtered junction list filename
SJ_LISTNAME="$SPECIES"

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

