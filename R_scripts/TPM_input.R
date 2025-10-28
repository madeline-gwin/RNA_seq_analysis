#Create TPM input

abbreviation<- "MG"
species<-"Marigold"
count_dir<-(file.path("C:/Users/madel/OneDrive - Clemson University/CapituLab/callus_RNA/read_counts",species))
setwd(count_dir)

##BCY##
D0_R1 <- read.delim(file.path(count_dir, "BCYD0R1ReadsPerGene.out.tab"), header=FALSE, sep="")
D0_R2 <- read.delim(file.path(count_dir, "BCYD0R2ReadsPerGene.out.tab"), header=FALSE, sep="")
D0_R3 <- read.delim(file.path(count_dir, "BCYD0R3ReadsPerGene.out.tab"), header=FALSE, sep="")

D2_R1 <- read.delim(file.path(count_dir, "BCYD2R1ReadsPerGene.out.tab"), header=FALSE, sep="")
D2_R2 <- read.delim(file.path(count_dir, "BCYD2R2ReadsPerGene.out.tab"), header=FALSE, sep="")
D2_R3 <- read.delim(file.path(count_dir, "BCYD2R3ReadsPerGene.out.tab"), header=FALSE, sep="")

D3_R1 <- read.delim(file.path(count_dir, "BCYD3R1ReadsPerGene.out.tab"), header=FALSE, sep="")
D3_R2 <- read.delim(file.path(count_dir, "BCYD3R2ReadsPerGene.out.tab"), header=FALSE, sep="")
D3_R3 <- read.delim(file.path(count_dir, "BCYD3R3ReadsPerGene.out.tab"), header=FALSE, sep="")

write.csv(D0_R1, file=file.path(count_dir, paste0(abbreviation,"D0R1_read_counts.csv")))
write.csv(D0_R2, file=file.path(count_dir, paste0(abbreviation,"D0R2_read_counts.csv")))
write.csv(D0_R3, file=file.path(count_dir, paste0(abbreviation,"D0R3_read_counts.csv")))
write.csv(D2_R1, file=file.path(count_dir, paste0(abbreviation,"D2R1_read_counts.csv")))
write.csv(D2_R2, file=file.path(count_dir, paste0(abbreviation,"D2R2_read_counts.csv")))
write.csv(D2_R3, file=file.path(count_dir, paste0(abbreviation,"D2R3_read_counts.csv")))
write.csv(D3_R1, file=file.path(count_dir, paste0(abbreviation,"D3R1_read_counts.csv")))
write.csv(D3_R2, file=file.path(count_dir, paste0(abbreviation,"D3R2_read_counts.csv")))
write.csv(D3_R3, file=file.path(count_dir, paste0(abbreviation,"D3R3_read_counts.csv")))


##Marigold##
D0_R1 <- read.delim(file.path(count_dir, "MGD0R1ReadsPerGene.out.tab"), header=FALSE, sep="")
D0_R2 <- read.delim(file.path(count_dir, "MGD0R2ReadsPerGene.out.tab"), header=FALSE, sep="")
D0_R3 <- read.delim(file.path(count_dir, "MGD0R3ReadsPerGene.out.tab"), header=FALSE, sep="")

D3_R1 <- read.delim(file.path(count_dir, "MGD3R1ReadsPerGene.out.tab"), header=FALSE, sep="")
D3_R2 <- read.delim(file.path(count_dir, "MGD3R2ReadsPerGene.out.tab"), header=FALSE, sep="")
D3_R3 <- read.delim(file.path(count_dir, "MGD3R3ReadsPerGene.out.tab"), header=FALSE, sep="")

D4_R1 <- read.delim(file.path(count_dir, "MGD4R1ReadsPerGene.out.tab"), header=FALSE, sep="")
D4_R2 <- read.delim(file.path(count_dir, "MGD4R2ReadsPerGene.out.tab"), header=FALSE, sep="")
D4_R3 <- read.delim(file.path(count_dir, "MGD4R3ReadsPerGene.out.tab"), header=FALSE, sep="")

write.csv(D0_R1, file=file.path(count_dir, paste0(abbreviation,"D0R1_read_counts.csv")))
write.csv(D0_R2, file=file.path(count_dir, paste0(abbreviation,"D0R2_read_counts.csv")))
write.csv(D0_R3, file=file.path(count_dir, paste0(abbreviation,"D0R3_read_counts.csv")))
write.csv(D3_R1, file=file.path(count_dir, paste0(abbreviation,"D3R1_read_counts.csv")))
write.csv(D3_R2, file=file.path(count_dir, paste0(abbreviation,"D3R2_read_counts.csv")))
write.csv(D3_R3, file=file.path(count_dir, paste0(abbreviation,"D3R3_read_counts.csv")))
write.csv(D4_R1, file=file.path(count_dir, paste0(abbreviation,"D4R1_read_counts.csv")))
write.csv(D4_R2, file=file.path(count_dir, paste0(abbreviation,"D4R2_read_counts.csv")))
write.csv(D4_R3, file=file.path(count_dir, paste0(abbreviation,"D4R3_read_counts.csv")))


##Lettuce##
D0_R1 <- read.delim(file.path(count_dir, "LSD0R1ReadsPerGene.out.tab"), header=FALSE, sep="")
D0_R2 <- read.delim(file.path(count_dir, "LSD0R2ReadsPerGene.out.tab"), header=FALSE, sep="")
D0_R3 <- read.delim(file.path(count_dir, "LSD0R3ReadsPerGene.out.tab"), header=FALSE, sep="")

D2_R1 <- read.delim(file.path(count_dir, "LsD2R1ReadsPerGene.out.tab"), header=FALSE, sep="")
D2_R2 <- read.delim(file.path(count_dir, "LsD2R2ReadsPerGene.out.tab"), header=FALSE, sep="")
D2_R3 <- read.delim(file.path(count_dir, "LsD2R3ReadsPerGene.out.tab"), header=FALSE, sep="")

D3_R1 <- read.delim(file.path(count_dir, "LsD3R1ReadsPerGene.out.tab"), header=FALSE, sep="")
D3_R2 <- read.delim(file.path(count_dir, "LsD3R2ReadsPerGene.out.tab"), header=FALSE, sep="")
D3_R3 <- read.delim(file.path(count_dir, "LsD3R3ReadsPerGene.out.tab"), header=FALSE, sep="")

write.csv(D0_R1, file=file.path(count_dir, paste0(abbreviation,"D0R1_read_counts.csv")))
write.csv(D0_R2, file=file.path(count_dir, paste0(abbreviation,"D0R2_read_counts.csv")))
write.csv(D0_R3, file=file.path(count_dir, paste0(abbreviation,"D0R3_read_counts.csv")))
write.csv(D2_R1, file=file.path(count_dir, paste0(abbreviation,"D2R1_read_counts.csv")))
write.csv(D2_R2, file=file.path(count_dir, paste0(abbreviation,"D2R2_read_counts.csv")))
write.csv(D2_R3, file=file.path(count_dir, paste0(abbreviation,"D2R3_read_counts.csv")))
write.csv(D3_R1, file=file.path(count_dir, paste0(abbreviation,"D3R1_read_counts.csv")))
write.csv(D3_R2, file=file.path(count_dir, paste0(abbreviation,"D3R2_read_counts.csv")))
write.csv(D3_R3, file=file.path(count_dir, paste0(abbreviation,"D3R3_read_counts.csv")))






