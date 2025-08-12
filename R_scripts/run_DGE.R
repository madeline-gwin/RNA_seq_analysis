# Name: run_DGE.R
# Author: MG
# Date: 08/06/2025
# Version:4.3.1
# Description: Will run the differential gene expression on samples

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("DESeq2")
library(DESeq2)
#BiocManager::install("Glimma")
library(Glimma)
#BiocManager::install("sva")
library(sva)
library(dplyr)
library(edgeR)
library(ggplot2)

# read in and process data
species<-"Marigold"
DE_dir<-(file.path("C:/Users/madel/OneDrive - Clemson University/CapituLab/callus_RNA/read_counts",species,"DESeq_output"))
setwd(DE_dir)

# read in the data matrix
summed_counts<-readRDS(file.path(DE_dir,"dds_set.Rdata"))

samples=c(summed_counts$condition)
print(samples)

dev_stage=c("T1","T1","T1","T2","T2","T2","T3","T3","T3")
#dev_stage=c("T1","T1","T1","T2","T2","T2","T3","T3")
dev_stage <- as.factor(dev_stage)

metadata<-data.frame(samples, dev_stage)

# create the factors of interest
metadata$dev_stage<-factor(metadata$dev_stage)

# create the model (wrt dev_stage)
summed_counts<-DESeqDataSetFromMatrix(counts(summed_counts),colData = metadata, design=~0+dev_stage)

# pre-filter for reads where at least 3 samples have a count of 1 or higher
keep<-rowSums(counts(summed_counts)>=1)>=2
length(which(keep==1))
summed_counts_filt<-summed_counts[keep,]

# plot MDS of adjusted counts...can compare this with previous plot pre- filtering for 0's
glimmaMDS(summed_counts_filt, group=metadata, col=c("#E97F89", "#6DA45F", "#E1CD86", "#E0A3D1"))

#make a PCA plot of filtered data
# Assume dds is your DESeqDataSet object
# Transform counts (rlog or vst)
rld <- rlog(summed_counts_filt, blind = TRUE)  # or use vst(dds_set)

# Perform PCA
pca_data <- plotPCA(rld, intgroup = "dev_stage", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Plot using ggplot2
#dev.new()
#print(
ggplot(pca_data, aes(PC1, PC2, color = dev_stage)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  labs(title="PCA of Gene Expression Data After Filtering")
#)

ggsave("pca_plot_summed_counts_filt.png", width = 6, height = 5, dpi = 300)

###################################################################
###################################################################
# run the DGE analysis
DESeq_dataset_results<-DESeq(summed_counts_filt,parallel=FALSE)


#if having problems with 0s
#res <- results(DESeq_dataset_results)
#res$log10pvalue <- -log10(res$pvalue)

# save DGE results to the output pathway
saveRDS(DESeq_dataset_results, file=file.path(DE_dir,"deseq_dataset_results.RData"))

# set up the pairwise contrasts (IvII, IIvIII,IvIII) play w alpha p values and see how data changes##
result_T2_vs_T1 <- results(DESeq_dataset_results, name = "dev_stageT2", contrast = c("dev_stage", "T2", "T1"))
result_T3_vs_T1 <- results(DESeq_dataset_results, name = "dev_stageT3", contrast = c("dev_stage", "T3", "T1"))
result_T3_vs_T2 <- results(DESeq_dataset_results, contrast = c("dev_stage", "T3", "T2"))

# write to CSV file
write.csv(as.data.frame(result_T2_vs_T1), file=file.path(DE_dir,"result_T2_vs_T1.csv"))
write.csv(as.data.frame(result_T3_vs_T1), file=file.path(DE_dir,"result_T3_vs_T1.csv"))
write.csv(as.data.frame(result_T3_vs_T2), file=file.path(DE_dir,"result_T3_vs_T2.csv"))
