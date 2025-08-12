# Name: analyze DGE deseq sunflower inflo ruvseq
# Author: ELL (based off EY code)
# Date: 08/30/2024
# Version:4.2.1
# Description: Will analyze the output from the DESeq DGE with combatseq for pairwise dev stages
# need Functions.R written by ED

library(dplyr)
library(ggplot2)
library(UpSetR)
library(Glimma)
source("C:/Users/madel/OneDrive - Clemson University/CapituLab/callus_RNA/scripts/Functions.R")

species<-"Lettuce"
DE_dir<-(file.path("C:/Users/madel/OneDrive - Clemson University/CapituLab/callus_RNA/read_counts",species,"DESeq_output"))
setwd(DE_dir)

# now analyze 
# read in the data, CHANGE PATH ERIKA To where results are on computer, results in a folder on their own, all csv in folder will
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
png("sequential_pairwise_upset.png", res=215, width = 1800, height=1000)
upset(fromList(SigOverlapGraph_pairwise_cs),order.by="freq",nsets=13,nintersects=20, text.scale = 1.5)
dev.off()

T2_v_T1 <- as.data.frame(mydataSig_pairwise_cs$result_T2_vs_T1)
T3_v_T1 <- as.data.frame(mydataSig_pairwise_cs$result_T3_vs_T1)
T3_v_T2 <- as.data.frame(mydataSig_pairwise_cs$result_T3_vs_T2)
InCommonAll <- as.data.frame(SigOverlap_pairwise_cs$InCommonAll)
print(InCommonAll)

write.csv(T2_v_T1, file='pairwise_T2_vs_T1.csv')
write.csv(T3_v_T1, file='pairwise_T3_vs_T1.csv')
write.csv(T3_v_T2, file='pairwise_T3_vs_T2.csv')

write.csv(SigOverlap_pairwise_cs[["InCommonAll"]], file= 'InCommonAll.csv')
