#get TPM from gene count data

#load packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GenomicFeatures", "GenomicRanges", "rtracklayer", "DESeq2", "edgeR"))

library(rtracklayer)

#Define Variables
#define species
species<-"Marigold"
abbreviation<- "MG"
day<-"T3R3"

# Define your directory containing the count files
count_dir<-(file.path("C:/Users/madel/OneDrive - Clemson University/CapituLab/callus_RNA/read_counts",species))
setwd(count_dir)

file_to_read <- list.files(pattern = "D4R3_read_counts.csv$", full.names = TRUE)

output_dir <- "C:/Users/madel/OneDrive - Clemson University/CapituLab/callus_RNA/output_tables/TPM/"

#annotation file load 
annotation_file="Tagetes_erecta.scaffold.genes.v2.gff3"
#annotation_file<-"Hildegardhap1v1_2 (1).gff"
#annotation_file<-"bidensv5hap1softreorderbig12 (1).gff"
annotation <- import(file.path("C:/Users/madel/OneDrive - Clemson University/CapituLab/callus_RNA/genome/Annotation_files/",annotation_file))

#####################################################################
#counts <- read.delim(file_to_read, row.names = 1, header = TRUE, sep=",")
counts <- read.csv(file_to_read, header = TRUE)
head(counts)

rownames(counts) <- counts$V1
counts <- counts["V2"]
colnames(counts) <- "Counts"
head(counts)


# Extract gene lengths (adjust based on your annotation structure)
gene_lengths <- width(annotation)
cds_regions <- annotation[annotation$type == "CDS"]
head(mcols(cds_regions))
# Use the appropriate attribute (gene_id, Parent, or ID)
# Many GFF files use "Parent" to refer to mRNA, which links to a gene
if ("gene_id" %in% colnames(mcols(cds_regions))) {
  gene_ids <- cds_regions$gene_id
} else if ("Parent" %in% colnames(mcols(cds_regions))) {
  gene_ids <- cds_regions$Parent  # This usually refers to mRNA, which links to a gene
} else {
  stop("No gene ID or Parent attribute found in the GFF file.")
}

gene_ids <- as.character(gene_ids)
valid_idx <- !is.na(gene_ids) & gene_ids != ""
cds_widths <- as.numeric(width(cds_regions))
cds_lengths <- tapply(cds_widths[valid_idx], gene_ids[valid_idx], sum, na.rm = TRUE)

cds_lengths <- unlist(cds_lengths)
head(cds_lengths)

# Subset CDS lengths to match genes in adjusted counts
common_genes <- intersect(rownames(counts), names(cds_lengths))

# Keep only matching genes
#counts <- counts[common_genes, ]
# Force counts to be a matrix even if there’s only one column
counts <- as.matrix(counts[common_genes, , drop = FALSE])

#cds_lengths <- cds_lengths[common_genes]

# Ensure cds_lengths is in the same order as counts
cds_lengths <- cds_lengths[rownames(counts)]

# Convert to a column matrix
cds_lengths_mat <- matrix(cds_lengths, nrow = length(cds_lengths), ncol = ncol(counts))

# Now divide
counts_per_kb <- counts / (cds_lengths_mat / 1000)

# Step 1: Calculate Counts Per Kilobase (CPK)
#counts_per_kb <- counts / (cds_lengths / 1000)

# Step 2: Calculate Scaling Factor (sum of CPK per sample)
scaling_factors <- colSums(counts_per_kb)

# Step 3: Normalize to TPM
TPM <- t(t(counts_per_kb) / scaling_factors * 1e6)

# Check output
head(TPM)

filename <- file.path(output_dir, paste0(abbreviation, "_", day, "_TPM_results.csv"))
write.csv(TPM, file = filename, row.names = TRUE)
