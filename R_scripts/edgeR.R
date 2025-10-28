if (!requireNamespace("BiocManager", quietly = TRUE))
+     install.packages("BiocManager")

BiocManager::install("edgeR")

library(edgeR)

setwd("/Users/rns0027/Library/CloudStorage/Box-Box/CapituLab general/Reid/projects/dev RNA Ls capitula/total RNA/ShortStack/run 1")

counts <- read.delim("Counts.txt", comment.char = "#")

count_matrix <- counts[, grep("S[1-5]_", colnames(counts))]

# If the known_miRNA column has a label, I want to name it after that, if it says "NA", I want to to name it after the Name column which is all Cluster IDs #

results <- read.delim("Results.txt", comment.char = "#")

rownames(count_matrix) <- make.unique(ifelse(!is.na(results$known_miRNAs), results$known_miRNAs, results$Name))

library(edgeR)
dge <- DGEList(counts = count_matrix)

# Filtering out lowly expressed loci using edgeR's smart filtration #

group <- factor(c(
  rep("S1", 3),  # S1_1, S1_2, S1_3
  rep("S2", 3),  # S2_1, S2_2, S2_3
  rep("S3", 3),  # S3_1, S3_2, S3_3
  rep("S4", 2),  # S4_2, S4_3 (S4_1 failed)
  rep("S5", 3)   # S5_1, S5_2, S5_3
))

dge$samples$group <- group

keep <- filterByExpr(dge, group = group)

dge <- dge[keep, , keep.lib.sizes = FALSE]

cat("Kept", sum(keep), "of", length(keep), "loci.\n")
# Kept 81429 of 89956 loci. #

# Now I want to normalize with TMM (Trimmed Mean of M-values), edgeR's default #

dge <- calcNormFactors(dge)

# Now I need to estimate dispersion (biological variability ) #

# Now I want to fit a model #

# define design matrix #
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# re-estimate dispersion (biological variability) with this design in mind instead of blindly across samples #
dge <- estimateDisp(dge, design)

# fitting a generalized linear model #
fit <- glmFit(dge, design)

### CIRCLE BACK TO THIS ###


# creating a PCA plot #

install.packages("ggplot2")

library(ggplot2)

logCPM <- cpm(dge, log = TRUE, prior.count = 1)

pca <- prcomp(t(logCPM), scale. = TRUE)






pca_df <- data.frame(PC1 = pca$x[,1],
                     PC2 = pca$x[,2],
                     group = group)

# Defining colors that match my mRNA PCA
custom_colors <- c(
  "S1" = "#9fd653",  # light green
  "S2" = "#e282be",  # light pink
  "S3" = "#1620f3",  # blue
  "S4" = "#f98c6a",  # salmon
  "S5" = "#e9f200"   # yellow
)

# Plot with manual colors
ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Small RNA Hilde_ShortStack_v1",
       x = "PC1", y = "PC2") +
  theme_classic()


# Now I want to find DEGs #

library(limma)

dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

contrast_matrix <- makeContrasts(
  S2_vs_S1 = S2 - S1,
  S3_vs_S1 = S3 - S1,
  S4_vs_S1 = S4 - S1,
  S5_vs_S1 = S5 - S1,
  S3_vs_S2 = S3 - S2,
  S4_vs_S2 = S4 - S2,
  S5_vs_S2 = S5 - S2,
  S4_vs_S3 = S4 - S3,
  S5_vs_S3 = S5 - S3,
  S5_vs_S4 = S5 - S4,
  levels = design
)

de_results_list <- list()

for (contrast_name in colnames(contrast_matrix)) {
  qlf <- glmQLFTest(fit, contrast = contrast_matrix[, contrast_name])
  de_results_list[[contrast_name]] <- topTags(qlf, n = Inf)$table
}

for (name in names(de_results_list)) {
  write.csv(de_results_list[[name]],
            paste0("DE_results_", name, ".csv"),
            row.names = TRUE)
}

setwd("/Users/rns0027/Library/CloudStorage/Box-Box/CapituLab general/Reid/projects/dev RNA Ls capitula/total RNA/ShortStack/run 1/edgeR")

# volcano plots #
library(ggplot2)

make_volcano_plot <- function(results_table, comparison_name, lfc_threshold = 1, fdr_threshold = 0.05) {
  results_table$significance <- with(results_table, ifelse(FDR < fdr_threshold & abs(logFC) > lfc_threshold,
                                                           "Significant", "Not Significant"))

  ggplot(results_table, aes(x = logFC, y = -log10(FDR), color = significance)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
    labs(title = paste("Volcano Plot:", comparison_name),
         x = "log2 Fold Change", y = "-log10 FDR") +
    theme_classic()
}

for (name in names(de_results_list)) {
  print(make_volcano_plot(de_results_list[[name]], name))
}

# Now I want a list of all the SDEG #

lfc_threshold <- 1
fdr_threshold <- 0.05

sig_genes_list <- list()

for (name in names(de_results_list)) {
  result <- de_results_list[[name]]
  sig <- result[result$FDR < fdr_threshold & abs(result$logFC) > lfc_threshold, ]
  sig_genes_list[[name]] <- rownames(sig)
  
  write.table(sig_genes_list[[name]],
              file = paste0("significant_genes_", name, ".txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# Now I'm going to take the lists of DELs and combine them into one .xlsx #

install.packages("openxlsx")
library(openxlsx)

file_info <- file.info(file_list)
empty_files <- rownames(file_info[file_info$size == 0, ])
print(empty_files)

non_empty_files <- file_list[file.info(file_list)$size > 0]

wb <- createWorkbook()

for (file in non_empty_files) {
  sheet_name <- gsub("^significant_genes_|\\.txt$", "", file)
  df <- read.delim(file, stringsAsFactors = FALSE)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, x = df)
}

saveWorkbook(wb, "All_Significant_Genes.xlsx", overwrite = TRUE)

# Now I want to generate upset plots #

install.packages("UpSetR")

library(UpSetR)

de_lists <- lapply(de_results_list, function(df) {
  rownames(df[df$FDR < 0.05, ])
})

input_upset <- fromList(de_lists)

upset(input_upset, order.by = "freq", nintersects = 30)






# I made a list of any locus that is differentially expressed between any stage, then created label_map.csv with this list and created a condensed version of the name ()(e.g. Ath-miR156a;Fan-miR156j;lsa-miR156i;Lsa-miR156h > miR156) #
# Not I want to create expression plots of all of them #

### Something went wrong here with naming I think, I'm starting over ###
label_map <- read.csv("label_map.csv", stringsAsFactors = FALSE)

label_map$known_miRNAs <- trimws(label_map$known_miRNAs)
results$known_miRNAs <- trimws(results$known_miRNAs)

merged <- merge(label_map, results, by = "known_miRNAs")

plot_lookup <- merged[, c("Name", "condensed_name")]

clusters_to_plot <- plot_lookup$Name[plot_lookup$Name %in% rownames(logCPM)]

filtered_lookup <- plot_lookup[plot_lookup$Name %in% clusters_to_plot, ]

plot_expr <- logCPM[filtered_lookup$Name, , drop = FALSE]

rownames(plot_expr) <- filtered_lookup$condensed_name[match(rownames(plot_expr), filtered_lookup$Name)]

install.packages("tidyr")

library(tidyr)
library(dplyr)

plot_df <- as.data.frame(t(plot_expr))
plot_df$sample <- rownames(plot_df)
plot_df$stage <- group  # assuming 'group' is already defined


plot_long <- pivot_longer(plot_df, cols = -c(sample, stage), names_to = "locus", values_to = "logCPM")

library(ggplot2)

# Summary stats: mean ± standard error
plot_summary <- plot_long %>%
  group_by(stage, gene) %>%
  summarize(mean_logCPM = mean(logCPM),
            se_logCPM = sd(logCPM) / sqrt(n()),
            .groups = "drop")

# Plot
ggplot(plot_summary, aes(x = stage, y = mean_logCPM, color = gene, group = gene)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_logCPM - se_logCPM,
                    ymax = mean_logCPM + se_logCPM),
                width = 0.2) +
  labs(title = "Expression Patterns Across Developmental Stages",
       x = "Stage",
       y = "Mean log2 CPM ± SE") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
### Something went wrong here with naming I think, I'm starting over ###


library(edgeR)
library(ggplot2)
library(tidyr)
library(dplyr)

setwd("/Users/rns0027/Library/CloudStorage/Box-Box/CapituLab general/Reid/projects/dev RNA Ls capitula/total RNA/ShortStack/run 1/edgeR")

results <- read.delim("Results.txt", comment.char = "#")
counts <- read.delim("Counts.txt", row.names = 2)
label_map <- read.csv("label_map.csv", stringsAsFactors = FALSE)

label_map$known_miRNAs <- trimws(label_map$known_miRNAs)
results$known_miRNAs <- trimws(results$known_miRNAs)

merged <- merge(label_map, results, by = "known_miRNAs")
plot_lookup <- merged[, c("Name", "condensed_name")]

dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log = TRUE)


custom_colors <- c(
  "S1" = "#9fd653",  # light green
  "S2" = "#e282be",  # light pink
  "S3" = "#1620f3",  # blue
  "S4" = "#f98c6a",  # salmon
  "S5" = "#e9f200"   # yellow
)

