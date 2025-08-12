# RNA_seq_analysis
Built from [Emily Yaklich's dev_RNAseq pipeline](https://github.com/emilyyaklich/dev_RNAseq). 
## Program Information
FastQC
MultiQC
FastP
Star
DESeq2
# Part 1: Read Mapping and Quantification
## Checking Downloaded Raw Reads for Corruption
```
# compile md5sum into a file
## print MD5sum for files from sequencing company into check.txt 
cat MD5.txt >> check.txt

# get md5sum for every raw read file downloaded from sequencing company and add to check.txt
## this should include forward and reverse reads
md5sum 01.RawData/D0R1/* >> check.txt
md5sum 01.RawData/D0R2/* >> check.txt
md5sum 01.RawData/D0R3/* >> check.txt
md5sum 02.Report_sampleID.zip >> check.txt

#check that only 1 unique md5sum exists per file
## sorts lines in check.txt, removes all duplicate lines, and sorts by the second column
## if corruption had occured, there would be two entries for one file with two different md5sums
sort check.txt | uniq | sort -k2
```
## Pre-processing
#### Check Quality of Raw Reads Before Filtering
Run [fastqc.sh](https://github.com/madeline-gwin/RNA_seq_analysis/blob/main/shell_scripts/fastqc.sh) to generate FastQC files and a combined MuitiQC file. This information can help guide you through filtering and adapter trimming processes. 
#### Trimming
Run [fastp.sh](https://github.com/madeline-gwin/RNA_seq_analysis/blob/main/shell_scripts/fastp.sh) to trim adapter sequences and low quality reads. 
#### Check Quality of Reads After Filtering
Run [fastqc.sh](https://github.com/madeline-gwin/RNA_seq_analysis/blob/main/shell_scripts/fastqc.sh) to generate FastQC files and a combined MuitiQC file to verify that your trimming did what it needed to. 

## Genome Indexing
Run [star_indexing.sh](https://github.com/madeline-gwin/RNA_seq_analysis/blob/main/shell_scripts/star_indexing.sh) to generate a genome index. 

## Two Pass Mapping
Two pass mapping is the gold standard in RNA sequencing analysis because it increases confidence of reads being mapped across exon-exon junctions.
#### Collect Junctions
Run [star_mapping.sh](https://github.com/madeline-gwin/RNA_seq_analysis/blob/main/shell_scripts/star_mapping.sh) on the **first pass** option. 
#### Filter Junctions
Run [filter_junctions.sh](https://github.com/madeline-gwin/RNA_seq_analysis/blob/main/shell_scripts/filter_junctions.sh).
#### Map Reads
Run [star_mapping.sh](https://github.com/madeline-gwin/RNA_seq_analysis/blob/main/shell_scripts/star_mapping.sh) on the **second pass** option.

## Prepare DESeq2 Input Files
The output from read mapping is too complex for DESeq2, so this script extracts only the gene ID's and the read counts.
Run [prepare_DEseq_input.sh](https://github.com/madeline-gwin/RNA_seq_analysis/blob/main/shell_scripts/prepare_DEseq_input.sh)

# Part 2: Expression Analysis
Expression analysis occurs exclusively using R. 

## Load Gene Count Data into R
Run [load_GC_data.R](https://github.com/madeline-gwin/RNA_seq_analysis/blob/main/R_scripts/load_GC_data.R) to read your gene count dataset into R, specifically as a DESeqDataSet. This script also generates a PCA plot to visualize variation between samples before any filtering.

## Perform Statistical Analysis using DESeq2
Run [run_DGE.R](https://github.com/madeline-gwin/RNA_seq_analysis/blob/main/R_scripts/run_DGE.R) to filter reads so that only reads with 1 or more occurences across at least 2 samples are included in the DESeqDataSet to be analyzed (summed_counts_filt). Another PCA plot is generated to show variation between samples from these filtered counts, allowing for comparison with the unfiltered gene counts. Then, DESeq will assign statistical values to the expression levels between the treatment conditions, creating DESeqDataSets that assign significance of differential gene expression between these conditions.

## Analyze Differential Gene Expression 
Run [analyze_DGE_deseq.R](https://github.com/madeline-gwin/RNA_seq_analysis/blob/main/R_scripts/analyze_DGE_deseq.R) to extract only significant differentially expressed genes and create an upset plot that shows the number of shared DEG's across sample conditions. This script requires [Functions.R](https://github.com/madeline-gwin/RNA_seq_analysis/blob/main/R_scripts/Functions.R) to run. 

## Generate a Volcano Plot
Run [volcano_plot.R](https://github.com/madeline-gwin/RNA_seq_analysis/blob/main/R_scripts/volcano_plot.R) to generate a volcano plot that shows the upregulated, downregulated, and non-significant genes. 
