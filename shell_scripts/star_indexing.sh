#!/bin/bash

#What file format is your genome annotation file? (GTF or GFF3)
ANNOTATION_FORMAT="GFF3"

#Where is the Genome Annotation file (GFF3 or GTF)?
GEN_ANN="/scratch/magwin/genome/Annotation_files/Tagetes_erecta.scaffold.genes.v2.gff3"

#Where is the Genome FASTA file?
GEN_FASTA="/scratch/magwin/genome/FASTA_files/Tagetes_erecta.scaffold.genome.v2.fa"

#Where do you want the files for your genome index?
GEN_DIR="/scratch/magwin/genome/genome_indices/Marigold"

#Specify the length of genomic sequence to be used in constructing the splice junctions database.
#This length should be equal to ReadLength-1, where ReadLength is the length of reads
#(ex. for 2x100bp paired-end reads, the ideal value is 99)
SPLICE_JUN=149

set -o pipefail

# read in the star module
ml STAR

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

