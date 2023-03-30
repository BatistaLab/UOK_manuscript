##some of the resouces that I have used
##https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/05_Annotation_and_Visualisation.nb.html#overview
##http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

#####################################################################################################
#libraries required for this session
library(DESeq2)
library(biomaRt)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")
library("gplots")
#####################################################################################################

#DESeq2 analysis from HTseq tablesZ
setwd("/Users/batistapj/Desktop/r2_yes.counts")

sampleFiles <- list.files(pattern="*.count")       #load files, make sure no other files include the word count
sampleFiles                                        # this checks the order of the files
status <- factor(c(rep("Min",3), rep("Mrpf",3), rep("Win",3), rep("Wrpf",3)))
sampleTable <- data.frame (sampleName = sampleFiles, fileName = sampleFiles, status=status)
directory <- "/Users/batistapj/Desktop/r2_yes.counts"
head(sampleTable)
des <- formula(~status)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design= des)

#filter the dataset (this is more stringent than removing all 0)
#filters out genes that don't have at least 5 normalized reads in 3 samples
#removed the filter using the normalized reads - values on the table were odd, almost every row had the same value

nrow(ddsHTSeq)
keep <- rowSums(counts(ddsHTSeq)) >= 1
dds <- ddsHTSeq[keep,]
#counts_RNA <- (counts(dds, normalize=FALSE))
#write.csv (as.data.frame(counts_RNA), file="idx3.csv")
nrow(dds)
head(dds)

#####################################################################################################
# Diferential expression analysis
dds <- DESeq(dds) #creates the analysis

#####################################################################################################
#####################################################################################################

# Extracting transformed values - for heatmaps and scatter plots

rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds)
rlogMat <- assay (rld)
vstMat <- assay (vsd)
#write.csv (as.data.frame(rlogMat), file="ribo2.csv")

counts_RNA <- (counts(dds, normalize=T))
write.csv (as.data.frame(counts_RNA), file="ribo2_cntRNA_.csv")

# Sample comparisons
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
distsRL <- dist(t(assay(rld)))
mat <- as.matrix (distsRL) 
rownames(mat) <- colnames (mat) <- with (colData (dds), paste (status, sep = " : "))
hc <- hclust (distsRL)

heatmap.2(mat, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col = rev (hmcol))

plotPCA (rld,intgroup=c("status"))

plotCounts(dds, gene="ENSG00000210112.1", intgroup="status")



##########################################################################################################################################################################################################
##########################################################################################################################################################################################################

setwd("/Users/batistapj/Desktop/r2_yes.counts")
counts <- read.table("countstable_R.txt", header=T, row.names=1)

head(counts)
dim(counts)

assay <- factor(c(rep("Input",3), rep("RPF",3), rep("Input",3), rep("RPF",3)))
condition <- factor(c(rep("WT",6), rep("Mut",6)))
coldata <- data.frame(row.names=colnames(counts), condition, assay)
coldata

dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design= ~ assay + condition + assay:condition)
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

# Diferential expression analysis

dds <- DESeq(dds, test="LRT", reduced= ~ assay + condition) #creates the analysis
res <- results(dds)
summary(res)

resOrderedlfc <- res[order (res $ log2FoldChange),]
head(resOrderedlfc)
# Exploring and exporting results
plotMA(res, ylim=c(-10,10))
write.csv(as.data.frame(res), file="multi_rib.csv")

##########################################################################################################################################################################################################
##########################################################################################################################################################################################################

# Sample comparisons
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
distsRL <- dist(t(assay(rld)))
mat <- as.matrix (distsRL) 
rownames(mat) <- colnames (mat) <- with (colData (dds), paste (status, sep = " : "))
hc <- hclust (distsRL)

heatmap.2(mat, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col = rev (hmcol))

plotPCA (rld,intgroup=c("status"))

plotCounts(dds, gene="ENSG00000163209.15", intgroup="assay")


##########################################################################################################################################################################################################
##########################################################################################################################################################################################################
