# create tables to better analyze RNA seq data

###########################################################
#All genes#

species<-"BCY_hap1"
abbreviation<-"BCY"

DE_dir<-(file.path("C:/Users/madel/OneDrive - Clemson University/CapituLab/callus_RNA/read_counts",species,"DESeq_output"))
setwd(DE_dir)

output_dir<- "C:/Users/madel/OneDrive - Clemson University/CapituLab/callus_RNA/output_tables/All_genes/"

T3_vs_T1<- read.csv(file.path(DE_dir,"result_T3_vs_T1.csv"))
T2_vs_T1<- read.csv(file.path(DE_dir,"result_T2_vs_T1.csv"))
T3_vs_T2<- read.csv(file.path(DE_dir,"result_T3_vs_T2.csv"))

Gene_ID_3v1 <- c(T3_vs_T1$X)
Gene_ID_3v2 <- c(T3_vs_T2$X)
Gene_ID_2v1 <- c(T2_vs_T1$X)

# Gene|LFC
T3_vs_T1_LFC<- data.frame(Gene_ID_3v1, T3_vs_T1$log2FoldChange)
T2_vs_T1_LFC<- data.frame(Gene_ID_2v1, T2_vs_T1$log2FoldChange)
T3_vs_T2_LFC<- data.frame(Gene_ID_3v2, T3_vs_T2$log2FoldChange)

write.csv(T3_vs_T1_LFC, file=file.path(output_dir, paste0(abbreviation,"_All_T3_vs_T1_LFC.csv")))
write.csv(T2_vs_T1_LFC, file=file.path(output_dir, paste0(abbreviation,"_All_T2_vs_T1_LFC.csv")))
write.csv(T3_vs_T2_LFC, file=file.path(output_dir, paste0(abbreviation,"_All_T3_vs_T2_LFC.csv")))

# Gene|padj
T3_vs_T1_padj<- data.frame(Gene_ID_3v1, T3_vs_T1$padj)
T2_vs_T1_padj<- data.frame(Gene_ID_2v1, T2_vs_T1$padj)
T3_vs_T2_padj<- data.frame(Gene_ID_3v2, T3_vs_T2$padj)

write.csv(T3_vs_T1_padj, file=file.path(output_dir, paste0(abbreviation,"_All_T3_vs_T1_padj.csv")))
write.csv(T2_vs_T1_padj, file=file.path(output_dir, paste0(abbreviation,"_All_T2_vs_T1_padj.csv")))
write.csv(T3_vs_T2_padj, file=file.path(output_dir, paste0(abbreviation,"_All_T3_vs_T2_padj.csv")))

#Combined LFC and padj
T3_vs_T1_combined<- data.frame(Gene_ID_3v1, T3_vs_T1$log2FoldChange, T3_vs_T1$padj)
T2_vs_T1_combined<- data.frame(Gene_ID_2v1, T2_vs_T1$log2FoldChange, T2_vs_T1$padj)
T3_vs_T2_combined<- data.frame(Gene_ID_3v2, T3_vs_T2$log2FoldChange, T3_vs_T2$padj)

write.csv(T3_vs_T1_combined, file=file.path(output_dir, paste0(abbreviation,"_All_T3_vs_T1_combined.csv")))
write.csv(T2_vs_T1_combined, file=file.path(output_dir, paste0(abbreviation,"_All_T2_vs_T1_combined.csv")))
write.csv(T3_vs_T2_combined, file=file.path(output_dir, paste0(abbreviation,"_All_T3_vs_T2_combined.csv")))


#############################################################
#DE genes only
species<-"Marigold"
abbreviation<-"MG"

DE_dir<-(file.path("C:/Users/madel/OneDrive - Clemson University/CapituLab/callus_RNA/read_counts",species,"DESeq_output"))
setwd(DE_dir)

output_dir<- "C:/Users/madel/OneDrive - Clemson University/CapituLab/callus_RNA/output_tables/DE_only"

T3_vs_T1<- read.csv(file.path(DE_dir,"pairwise_T3_vs_T1.csv"))
T2_vs_T1<- read.csv(file.path(DE_dir,"pairwise_T2_vs_T1.csv"))
T3_vs_T2<- read.csv(file.path(DE_dir,"pairwise_T3_vs_T2.csv"))

Gene_ID_3v1 <- c(T3_vs_T1$Gene)
Gene_ID_3v2 <- c(T3_vs_T2$Gene)
Gene_ID_2v1 <- c(T2_vs_T1$Gene)

# Gene|LFC
T3_vs_T1_LFC<- data.frame(Gene_ID_3v1, T3_vs_T1$log2FoldChange)
T2_vs_T1_LFC<- data.frame(Gene_ID_2v1, T2_vs_T1$log2FoldChange)
T3_vs_T2_LFC<- data.frame(Gene_ID_3v2, T3_vs_T2$log2FoldChange)

write.csv(T3_vs_T1_LFC, file=file.path(output_dir, paste0(abbreviation,"_DE_T3_vs_T1_LFC.csv")))
write.csv(T2_vs_T1_LFC, file=file.path(output_dir, paste0(abbreviation,"_DE_T2_vs_T1_LFC.csv")))
write.csv(T3_vs_T2_LFC, file=file.path(output_dir, paste0(abbreviation,"_DE_T3_vs_T2_LFC.csv")))

# Gene|padj
T3_vs_T1_padj<- data.frame(Gene_ID_3v1, T3_vs_T1$padj)
T2_vs_T1_padj<- data.frame(Gene_ID_2v1, T2_vs_T1$padj)
T3_vs_T2_padj<- data.frame(Gene_ID_3v2, T3_vs_T2$padj)

write.csv(T3_vs_T1_padj, file=file.path(output_dir, paste0(abbreviation,"_DE_T3_vs_T1_padj.csv")))
write.csv(T2_vs_T1_padj, file=file.path(output_dir, paste0(abbreviation,"_DE_T2_vs_T1_padj.csv")))
write.csv(T3_vs_T2_padj, file=file.path(output_dir, paste0(abbreviation,"_DE_T3_vs_T2_padj.csv")))

#Combined LFC and padj
T3_vs_T1_combined<- data.frame(Gene_ID_3v1, T3_vs_T1$log2FoldChange, T3_vs_T1$padj)
T2_vs_T1_combined<- data.frame(Gene_ID_2v1, T2_vs_T1$log2FoldChange, T2_vs_T1$padj)
T3_vs_T2_combined<- data.frame(Gene_ID_3v2, T3_vs_T2$log2FoldChange, T3_vs_T2$padj)

write.csv(T3_vs_T1_combined, file=file.path(output_dir, paste0(abbreviation,"_DE_T3_vs_T1_combined.csv")))
write.csv(T2_vs_T1_combined, file=file.path(output_dir, paste0(abbreviation,"_DE_T2_vs_T1_combined.csv")))
write.csv(T3_vs_T2_combined, file=file.path(output_dir, paste0(abbreviation,"_DE_T3_vs_T2_combined.csv")))

