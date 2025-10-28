# Name: plot_MDS
# Author: MG
# Date: 8/1/2025
# Version:4.2.1
# Description: will read the data matrix from load_GC_data_and_sum_reps into
# DESeq2, read in metadata,pre-filter and plotMDS

#BiocManager::install("readx1")
library(readxl)
#BiocManager::install("DESeq2")
library(DESeq2)
#BiocManager::install("dplyr")
library(dplyr)
#BiocManager::install("Glimma")
library(Glimma)
#BiocManager::install("sva")
library(sva)

species<-"Lettuce"
DE_dir<-(file.path("/scratch/magwin/genome/read_mapping",species,"DEseq"))
setwd(DE_dir)

# read in the data matrix
summed_counts<-readRDS(file.path(DE_dir,"dds_set.Rdata"))
dim(summed_counts)

#Define sample names, stage, and tissue type
samples=c(summed_counts$condition)
print(samples)
#dev_stage=c("T1","T1","T1","T2","T2","T2","T3","T3","T3")
dev_stage=c("T1","T1","T1","T2","T2","T2","T3","T3")
#tissue_type=c("callus","callus","callus","callus","callus","callus","callus","callus","callus")
tissue_type=c("callus","callus","callus","callus","callus","callus","callus","callus")

# create data frame
metadata<-data.frame(samples, tissue_type,dev_stage)

# create the factors of interest
metadata$dev_stage<-factor(metadata$dev_stage)
metadata$tissue_type<-factor(metadata$tissue_type)

# create the model
summed_counts <- DESeqDataSetFromMatrix(counts(summed_counts), colData = metadata, design = ~0 + dev_stage)

# pre-filter for reads where at least 3 samples have a count of 1 or higher
keep<-rowSums(counts(summed_counts)>=1)>=3
length(which(keep==1))
summed_counts_filt<-summed_counts[keep,]

# will load the plot...need to save within the html
glimmaMDS(summed_counts_filt, groups=metadata)
glimmaMDS(summed_counts, groups=metadata)
