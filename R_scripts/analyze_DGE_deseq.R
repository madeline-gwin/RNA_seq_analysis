# Name: analyze DGE deseq sunflower inflo ruvseq
# Author: ELL (based off EY code)
# Date: 08/30/2024
# Version:4.2.1
# Description: Will analyze the output from the DESeq DGE with combatseq for pairwise dev stages
# need Functions.R written by ED

#Download Packages
library(dplyr)
library(ggplot2)
library(UpSetR)
library(Glimma)
source("C:/Users/madel/OneDrive - Clemson University/CapituLab/callus_RNA/scripts/Functions.R")

#####################################################################
#Define Variables

#Define Species
species<-"BCY_hap1"

#Define Directory with DESeq Output
DE_dir<-(file.path("C:/Users/madel/OneDrive - Clemson University/CapituLab/callus_RNA/read_counts",species,"DESeq_output"))
setwd(DE_dir)
######################################################################
# now analyze 

##read##
DEData_pairwise_cs<-ImportCSVs(file.path(DE_dir),0.05)
# filter out significant results
mydataSig_pairwise_cs<-lapply(DEData_pairwise_cs,SigDEdf,PvaluesCol=7,CritP=0.05)

# see which genes overlap (input into upset plot)
SigOverlap_pairwise_cs <- GeneSets(
  mydataSig_pairwise_cs$result_T2_vs_T1[[1]],
  mydataSig_pairwise_cs$result_T3_vs_T1[[1]],
  mydataSig_pairwise_cs$result_T3_v_T2[[1]]
)
names(SigOverlap_pairwise_cs)
lapply(SigOverlap_pairwise_cs, length)
SigOverlapGraph_pairwise_cs<-lapply(mydataSig_pairwise_cs, function(x) {x$Gene})

# create an upset plot of DE expression by pairwise dev_stage
png("sequential_pairwise_upset.png", res=215, width =1500, height=800)
upset(fromList(SigOverlapGraph_pairwise_cs),order.by="freq",nsets=13,nintersects=20, text.scale = 1.5)
dev.off()

T2_vs_T1 <- as.data.frame(mydataSig_pairwise_cs$result_T2_vs_T1)
T3_vs_T1 <- as.data.frame(mydataSig_pairwise_cs$result_T3_vs_T1)
T3_vs_T2 <- as.data.frame(mydataSig_pairwise_cs$result_T3_vs_T2)

InCommon2<- intersect(T2_vs_T1$Gene, T3_vs_T1$Gene)
InCommonAll<- intersect(InCommon2, T3_vs_T2$Gene)                    

write.csv(x = InCommonAll, file = "InCommonAll.csv", row.names = FALSE)
write.csv(T2_vs_T1, file='pairwise_T2_vs_T1.csv')
write.csv(T3_vs_T1, file='pairwise_T3_vs_T1.csv')
write.csv(T3_vs_T2, file='pairwise_T3_vs_T2.csv')

