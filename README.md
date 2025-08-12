# RNA_seq_analysis
Built from [Emily Yaklich's dev_RNAseq pipeline](https://github.com/emilyyaklich/dev_RNAseq). 
## Program Information
FastQC
MultiQC
FastP
Star
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
### Check Quality of Raw Reads Before Filtering
Run [fastqc.sh](https://github.com/madeline-gwin/RNA_seq_analysis/blob/main/shell_scripts/fastqc.sh) to generate FastQC files and a combined MuitiQC file. This information can help guide you through filtering and adapter trimming processes. 
### Trimming
Run [fastp.sh](https://github.com/madeline-gwin/RNA_seq_analysis/blob/main/shell_scripts/fastp.sh) to trim adapter sequences and low quality reads. 
### Check Quality of Reads After Filtering
Run [fastqc.sh](https://github.com/madeline-gwin/RNA_seq_analysis/blob/main/shell_scripts/fastqc.sh) to generate FastQC files and a combined MuitiQC file to verify that your trimming did what it needed to. 

## Genome Indexing
Run [star_indexing.sh](https://github.com/madeline-gwin/RNA_seq_analysis/blob/main/shell_scripts/star_indexing.sh) to generate a genome index. 

## 2 Pass Mapping
Two pass mapping is the gold standard in RNA sequencing analysis because it increases confidence of reads being mapped across exon-exon junctions.
### Collect Junctions
Run [star_mapping.sh](https://github.com/madeline-gwin/RNA_seq_analysis/blob/main/shell_scripts/star_mapping.sh) on the **first pass** option. 
### Filter Junctions
Run [filter_junctions.sh](https://github.com/madeline-gwin/RNA_seq_analysis/blob/main/shell_scripts/filter_junctions.sh).
### Map Reads
Run [star_mapping.sh](https://github.com/madeline-gwin/RNA_seq_analysis/blob/main/shell_scripts/star_mapping.sh) on the **second pass** option.
