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

## get md5sum for every raw read file downloaded from sequencing company
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
