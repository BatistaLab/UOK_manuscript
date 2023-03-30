# Analysis of differential gene expression at the gene_input level for the m6A-IP data
# Christina Fitzsimmons
# Created 2021-03-09
# Updated 2021-04-23


# Set working directory 
setwd ("/Users/fitzsimmonscm/Documents/Projects_and_Data/Batista_Lab/FH_SDHB_Project/Rprojects/UOK_manuscript/m6A_IP/")
library(tidyverse)
library(QNB)
library(DESeq2)

#Telling deseq2 htseq-count script where to find the data and what it is called
sampleFiles <- list.files(pattern="*.count") #load files, make sure no other files include the word count
sampleFiles # check the order of the files

directory <- "./05_read_counting/gene_level_counts/input/"
sampleFiles <- list.files(path = "./05_read_counting/gene_level_counts/input/", pattern="*.count")       # load files
sampleFiles # this checks the order of the files

status <- factor(c(rep("WT",3), rep("MUT",3)))
sampleTable <- data.frame (sampleName = sampleFiles, fileName = sampleFiles, status=status)


head(sampleTable)
dds <- formula(~status)
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design= dds)

head(sampleFiles)  #check files
head(dds)     #check dds object

#filter the dataset - remove all zeros, influences padj
#In this case, we are removing things with less than 5 counts--but we're doing this after normalization

nrow(dds) #this is at the gene level. At the peak level analysis, we should have approx 25,000
dds <- estimateSizeFactors(dds) 
idx <- rowSums(counts(dds, normalized=FALSE) >= 5 ) >= 6
dds <- dds[idx,]
nrow(dds) #15068 genes in this analysis--would we expect more?

write.csv (as.data.frame(counts(dds)), file="./06_differential_peak_analysis/2021.04.23_gene_level_countsV2.csv")

# Differential expression analysis
dds <- DESeq(dds) #creates the analysis

# Plot counts - check counts of individual genes
counts_RNA <- (counts(dds, normalize=T))
res_ab = results(dds, contrast=c("status","MUT","WT"))
res_ab
summary(res_ab)

write.csv(as.data.frame(res_ab), file="/Users/fitzsimmonscm/Desktop/m6AIP_262counts/peak_counts_final_list/gene_level_counts/2021.03.09_gene-level-analysis_MutvsWT_deseq2-results.csv")

res_stats <- as.data.frame(res_ab)

# Plotting Things in volcano plot to see what our peaks look like. 
testvolcano <-  res_stats %>%
  mutate("logDESeqP" = -1*log10(padj)) %>%
  mutate(
    "direction" = case_when(
      log2FoldChange <= -1.5 & padj < 0.05 ~ "downreg",
      log2FoldChange >= 1.5 & padj < 0.05 ~ "upreg",
      TRUE ~ "no_change")) %>%
  dplyr::mutate(direction = factor(direction, levels = c("upreg", "downreg", "no_change")))

m6A_scatter <- ggplot(testvolcano, aes(x = log2FoldChange, y = logDESeqP, color=direction)) +
  geom_point() +
  theme_bw() +
  geom_vline(xintercept = -1.5, linetype="dashed") + geom_vline(xintercept = 1.5, linetype = "dashed") +
  geom_hline(yintercept = 1.30, linetype = "dashed") +
  scale_color_manual(values=c("red", "blue", "#999999")) +
  labs(x = "Log2FoldChange IP (Mut/WT)", y = "-Log10(P)", title = "Gene-level input Analysis of Mut vs WT") +
  xlim(-4,4)

m6A_scatter




