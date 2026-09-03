#!/bin/bash

#SBATCH --job-name=STAR_Hap2_GTF        # job name
#SBATCH --partition=capitulab           # partition
#SBATCH --ntasks=1                      # number of tasks across all nodes
#SBATCH --cpus-per-task=1               # number of cpus per task
#SBATCH --mem=35GB                      # total memory requested
#SBATCH --time=72:00:00                 # Run time (D-HH:MM:SS)
#SBATCH --output=job-%j.out             # Output file. %j is replaced with job ID
#SBATCH --error=job-%j.err              # Error file. %j is replaced with job ID
#SBATCH --mail-type=ALL                 # will send email for begin,end,fail
#SBATCH --mail-user=magwin@clemson.edu

#What file format is your genome annotation file? (GTF or GFF3)
ANNOTATION_FORMAT="GTF"

#Where is the Genome Annotation file (GFF3 or GTF)?
GEN_ANN="/scratch/magwin/BCY_floret_dev_analysis/genomes/Hap2_GTF/BidensCY_H2_genes.gtf"

#Where is the Genome FASTA file?
GEN_FASTA="/scratch/magwin/BCY_floret_dev_analysis/genomes/Hap2_GTF/BidensCY_H2.fasta"

#Where do you want the files for your genome index?
GEN_DIR="/scratch/magwin/BCY_floret_dev_analysis/genomes/Hap2_GTF/genome_indices/"

# mkake output directory, if it doesnt exist
mkdir -p $GEN_DIR

#Specify the length of genomic sequence to be used in constructing the splice junctions database.
#This length should be equal to ReadLength-1, where ReadLength is the length of reads
#(ex. for 2x100bp paired-end reads, the ideal value is 99)
SPLICE_JUN=149

set -o pipefail

# read in the star module
module load biocontainers
module load star

# Define Exon Parent Gene and Transcript Tags
if [ "$ANNOTATION_FORMAT" == "GFF3" ]; then
    echo "Using GFF3 annotation to generate a genome index"
    TRANSCRIPT_TAG="Parent"
    GENE_TAG="Parent"  # Or change based on your GFF3's attribute for genes
elif [ "$ANNOTATION_FORMAT" == "GTF" ]; then
    echo "Using GTF annotation to generate a genome index"
    TRANSCRIPT_TAG="transcript_id"
    GENE_TAG="gene_id"
else
    echo "Please specify whether annotation file is in GTF or GFF3 format"
    exit 1
fi

cd $GEN_DIR

echo "Beginning indexing"
    STAR \
--genomeSAindexNbases 13 \
--runMode genomeGenerate \
--genomeDir $GEN_DIR \
--genomeFastaFiles $GEN_FASTA \
--sjdbGTFtagExonParentTranscript $TRANSCRIPT_TAG \
--sjdbGTFtagExonParentGene $GENE_TAG \
--sjdbGTFfile $GEN_ANN \
--sjdbOverhang $SPLICE_JUN \
