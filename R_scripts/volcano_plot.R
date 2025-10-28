#Script derived from Biostatsquid
#https://biostatsquid.com/volcano-plots-r-tutorial/

# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
#install.packages("ggrepel")
library(ggrepel) # for nice annotations

# Set input path
species<-"BCY_hap1"
comparison <- "T3_vs_T1"
#comparison <- "T2_vs_T1"
#comparison <- "T3_vs_T2"
file_name <- paste0("result_", comparison, ".csv")

DE_dir<-(file.path("C:/Users/madel/OneDrive - Clemson University/CapituLab/callus_RNA/read_counts",species,"DESeq_output"))
setwd(DE_dir)

# Import DGE results
df <- read.csv(file.path(DE_dir,file_name), row.names = 1)


# Biostatsquid theme
theme_set(theme_classic(base_size = 20) +
              theme(
                axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
                axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
                plot.title = element_text(hjust = 0.5)
              ))


# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
df$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
df$diffexpressed[df$log2FoldChange > 0.6 & df$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df$diffexpressed[df$log2FoldChange < -0.6 & df$pvalue < 0.05] <- "DOWN"
# Create a new column "delabel" to de, that will contain the name of the top 30 differentially expressed genes (NA in case they are not)
#df$gene_symbol <- rownames(df)
#df$delabel <- ifelse(df$gene_symbol %in% head(df[order(df$padj), "gene_symbol"], 30), df$gene_symbol, NA)



ggplot(data = df, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "grey", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = 'dashed') + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("#088C0A", "grey", "#FFDB6D"), # to set the colours of our variable  
      labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 200), xlim = c(-15, 15)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Gene Expression', #legend_title, 
       x = expression("log"[2]*"FoldChange"), 
       y = expression("-log"[10]*"p-value"), 
       title = paste(comparison, species, "Callus")) + 
  scale_x_continuous(breaks = seq(-15, 15, 2))  # to customize the breaks in the x axis


#+ geom_text_repel(max.overlaps = Inf) # To show all labels 

#save plot
ggsave(paste0("volcano_plot_", comparison, ".png"), width = 10, height = 7, dpi = 300)

