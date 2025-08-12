#!/bin/bash

# Load STAR module (adjust to your environment if needed)
module load biocontainers
module load star

set -o pipefail

######################################
### HARD-CODED VARIABLE DEFINITIONS ###
######################################

# Queuing system: "Slurm" or "PBS"
QUEUE="Slurm"

# What is your species?
SPECIES="BCY"

# Input directory (contains trimmed reads)
RM_INPUT="/scratch/magwin/callus_RNA/pre-processing/after_filtering/$SPECIES"

# File suffixes for forward and reverse reads
FORWARD="_1_trimmed.fq.gz"
REVERSE="_2_trimmed.fq.gz"

# Output directories
#CJ_OUTPUTDIR="/scratch/magwin/genome/genome_indices/$SPECIES/Junctions"
CJ_OUTPUTDIR="/scratch/magwin/genome/genome_indices/BCY_hap1/Junctions"

#RM_OUTPUTDIR="/scratch/magwin/genome/read_mapping/$SPECIES/"
RM_OUTPUTDIR="/scratch/magwin/genome/read_mapping/BCY_hap1/"

mkdir -p "$CJ_OUTPUTDIR"
mkdir -p "$RM_OUTPUTDIR"

# Genome index directory
#GEN_DIR="/scratch/magwin/genome/genome_indices/$SPECIES/"
GEN_DIR="/scratch/magwin/genome/genome_indices/BCY_hap1/"

# Job log (used in SLURM mode)
##Pass1
#RM_JOB_LOG="/scratch/magwin/genome/genome_indices/$SPECIES/Junctions/job_log.txt"
#RM_JOB_LOG="/scratch/magwin/genome/genome_indices/BCY_hap1/Junctions/job_log.txt"
##Pass2
#RM_JOB_LOG="/scratch/magwin/genome/read_mapping/$SPECIES/job_log.txt"
RM_JOB_LOG="/scratch/magwin/genome/read_mapping/BCY_hap1/job_log.txt"

# Junctions file (if using second pass)
#JUNCTIONS="/scratch/magwin/genome/genome_indices/$SPECIES/Junctions/${SPECIES}_SJ.filtered.tab"
JUNCTIONS="/scratch/magwin/genome/genome_indices/BCY_hap1/Junctions/BCY_hap1_SJ.filtered.tab"

# STAR binary (only used for PBS)
STAR_FILE="/path/to/star/executable/if/using/PBS"

# STAR parameters
RM_PASS="second"              # "first" or "second"
RM_NTHREAD=8
SEEDSEARCH=50			#Emily = 50
MAX_MIS=10
MAX_N=10			#Emily = 10
MINSCORE_READL=0.66		
MINMATCH_READL=0.66
UNMAP_F="Fastx"              # Options: None, Within, Fastx
GENOMIC_COORDINATE_BAMSORTED="yes"
QUANT="TranscriptomeSAM"     # or "-" to disable
PLATFORM="ILLUMINA"


##################################
### JOB ARRAY ID & FILE HANDLING #
##################################

if [[ "${QUEUE}" == "Slurm" ]]; then
    PBS_ARRAYID=${SLURM_ARRAY_TASK_ID}
    echo "${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}" >> ${RM_JOB_LOG}
    echo "Processing array index ${PBS_ARRAYID}"
else
    echo "Processing array ${PBS_ARRAYID} through PBS queuing system"
fi

# Get list of forward reads
if [[ ! -d "$RM_INPUT" ]]; then
    echo "ERROR: RM_INPUT ($RM_INPUT) is not a directory."
    exit 1
fi

readarray -t all_forward_reads < <(find "$RM_INPUT" -name "*$FORWARD" | sort)

f1="${all_forward_reads[$PBS_ARRAYID]}"

if [[ -z "$f1" || ! -f "$f1" ]]; then
    echo "ERROR: Forward read file not found for task ID $PBS_ARRAYID"
    exit 1
fi

name=$(basename "${f1%%$FORWARD}")

#########################################
### EXTRACT FLOWCELL NAME AFTER f1 SET ###
#########################################

FLOWCELL_NAME=$(zcat "$f1" | head -n 1 | cut -d ':' -f 3)

###################################
### HANDLE PAIRED-END OR SINGLE ###
###################################

PE="True"  # Set to "False" for single-end reads

if [[ "$PE" == "True" ]]; then
    f2="${f1%%$FORWARD}$REVERSE"
    if [[ -f "$f2" ]]; then
        echo "Mapping PE reads for sample $name"
    else
        echo "ERROR: Paired read $f2 does not exist"
        exit 1
    fi
else
    f2=""
    echo "Mapping SE reads for sample $name"
fi

#########################################
### CREATE READ GROUP ID FOR SAM HEADER #
#########################################

LANE_NUM=$(grep -o "L00[1-4]" <<< "$name" | cut -c 4)
SAMPLE_NAME=${name%%[!0-9]*}
ID="${SAMPLE_NAME}:${FLOWCELL_NAME}.${LANE_NUM}"

echo "File name indicates the sample name is ${SAMPLE_NAME} and the lane number is ${LANE_NUM}"
echo "The read group ID field will be ${ID}"

##########################################
### CHOOSE OUTPUT FORMAT FOR ALIGNMENT ###
##########################################

if [[ "$GENOMIC_COORDINATE_BAMSORTED" == "yes" ]]; then
    FORMAT="BAM SortedByCoordinate"
    echo "Output Genomic Alignments will be sorted BAM files"
else
    FORMAT="SAM"
    echo "Output Genomic Alignments will be unsorted SAM files"
fi

# Choose STAR command
if [[ "${QUEUE}" == "Slurm" ]]; then
    STAR="STAR"
else
    STAR="${STAR_FILE}"
fi

##############################
### STAR ALIGNMENT COMMAND ###
##############################

if [[ "$RM_PASS" == "first" ]]; then
    echo "In first pass Mode"
    $STAR \
        --runThreadN $RM_NTHREAD \
        --genomeDir $GEN_DIR \
        --readFilesIn $f1 $f2 \
        --readFilesCommand gunzip -c \
        --seedSearchStartLmax $SEEDSEARCH \
        --outFileNamePrefix $CJ_OUTPUTDIR/"$name" \
        --outFilterMismatchNmax $MAX_MIS \
        --outFilterMultimapNmax $MAX_N \
        --outFilterScoreMinOverLread $MINSCORE_READL \
        --outFilterMatchNminOverLread $MINMATCH_READL \
        --outReadsUnmapped $UNMAP_F \
        --outSAMtype SAM \
        --quantMode - \
        --outSAMattrRGline ID:${ID} LB:${SAMPLE_NAME} PL:${PLATFORM} SM:${SAMPLE_NAME} PU:${ID} \
        --outFilterType BySJout \
        --outSJfilterReads Unique
elif [[ "$RM_PASS" == "second" ]]; then
    if [[ ! -z "$JUNCTIONS" ]]; then
        echo "In second pass mode using $NUM_JUNCTIONS junction files"
        echo "Junctions are as follows: $JUNCTIONS"
        $STAR \
            --runThreadN $RM_NTHREAD \
            --genomeDir $GEN_DIR \
            --readFilesIn $f1 $f2 \
            --readFilesCommand gunzip -c \
            --seedSearchStartLmax $SEEDSEARCH \
            --outFileNamePrefix $RM_OUTPUTDIR/"$name" \
            --outFilterMismatchNmax $MAX_MIS \
            --outFilterMultimapNmax $MAX_N \
            --outFilterScoreMinOverLread $MINSCORE_READL \
            --outFilterMatchNminOverLread $MINMATCH_READL \
            --outReadsUnmapped $UNMAP_F \
            --outSAMtype $FORMAT \
            --outSAMattributes NH HI AS nM NM MD \
            --quantMode TranscriptomeSAM GeneCounts\
            --outSAMattrRGline ID:${ID} LB:${SAMPLE_NAME} PL:${PLATFORM} SM:${SAMPLE_NAME} PU:${ID} \
            --outFilterType BySJout \
            --sjdbFileChrStartEnd $JUNCTIONS
    else
        echo "Mapping without incorporating un-annotated junctions"
        $STAR \
            --runThreadN $RM_NTHREAD \
            --genomeDir $GEN_DIR \
            --readFilesIn $f1 $f2 \
            --readFilesCommand gunzip -c \
            --seedSearchStartLmax $SEEDSEARCH \
            --outFileNamePrefix $RM_OUTPUTDIR/"$name" \
            --outFilterMismatchNmax $MAX_MIS \
            --outFilterMultimapNmax $MAX_N \
            --outFilterScoreMinOverLread $MINSCORE_READL \
            --outFilterMatchNminOverLread $MINMATCH_READL \
            --outReadsUnmapped $UNMAP_F \
            --outSAMtype $FORMAT \
            --outSAMattributes NH HI AS nM NM MD \
            --quantMode $QUANT \
            --outSAMattrRGline ID:${ID} LB:${SAMPLE_NAME} PL:${PLATFORM} SM:${SAMPLE_NAME} PU:${ID} \
            --outFilterType BySJout
    fi
else
    echo "Error: Unsure of whether first or second pass mode, exiting..."
    exit 1
fi

