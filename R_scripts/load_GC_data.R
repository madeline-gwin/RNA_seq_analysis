# Name: load GC data and sum reps
# Author: MG
# Date: 07/29/2025
# Description: will load in the parsed gene count data from hpc

#load DESeq2
#install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)
R.version

#This is a lot of data that needs to be read in, so 100GB allocated 

#define species
species<-"Lettuce"

# Define your directory containing the count files
count_dir<-(file.path("C:/Users/madel/OneDrive - Clemson University/CapituLab/callus_RNA/read_counts",species))
setwd(count_dir)

# List all files ending with ReadsPerGene.out.tab
GC_files <- list.files(path = count_dir, pattern = "ReadsPerGene.out.tab$", full.names = FALSE)
print (GC_files)

# Extract sample names by removing the suffix
sample_names <- sub("ReadsPerGene.out.tab", "", GC_files)
print(sample_names)

# Optional: extract condition from sample name (e.g., LSD0 from LSD0R1)
conditions <- sub("R[0-9]+$", "", sample_names)
print(conditions)

# Extract replicate number from sampleName by capturing digits after 'R'
replicates <- sub(".*R", "", sample_names)
print(replicates)

# Build the sample table
sampleTable <- data.frame(
  sampleName = sample_names,
  fileName = GC_files,
  condition = conditions,
  replicate = replicates
)
print(sampleTable)

# Create DESeq dataset
dds_set <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable,
  directory = count_dir,
  design = ~ condition
)

# Define path using species variable
output_dir <- file.path("C:/Users/madel/OneDrive - Clemson University/CapituLab/callus_RNA/read_counts",species,"DESeq_output")

# Save the data set
saveRDS(dds_set, file = file.path(output_dir, "dds_set.Rdata"))
dim(dds_set)

## Make a PCA plot of your un-filtered data
  # Assume dds is your DESeqDataSet object
  # Transform counts (rlog or vst)
rld <- rlog(dds_set, blind = TRUE)  # or use vst(dds_set)

# Perform PCA
pca_data <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Plot using ggplot2
dev.new()
print(
  ggplot(pca_data, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme_minimal() +
    labs(title="PCA of Gene Expression Data No Filtering")
)

ggsave(file.path(output_dir,"pca_plot_dds_set.png"), width = 6, height = 5, dpi = 300)
