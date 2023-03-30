# Analysis of differential gene expression at the peak IP level with gene-level correction
# Christina Fitzsimmons
# Last updated 2021-06-10

#### 1. Read in the data ####
# Set working directory 
setwd ("/Users/fitzsimmonscm/Documents/Projects_and_Data/Batista_Lab/FH_SDHB_Project/Rprojects/UOK_manuscript/m6A_IP/")
library(DESeq2)
library(ggpubr)
library(MASS)
library(plotly)
library(tidyverse)
library(viridis)
library(ggExtra)

#importing the data for conversions
name_conversion <- read.csv(file = "/Users/fitzsimmonscm/Documents/Projects_and_Data/Batista_Lab/FH_SDHB_Project/Rprojects/UOK_manuscript/longest_isoform/gProfiler_hsapiens_2021_03_10_conversion.csv", header = TRUE, sep = ",")%>%
  dplyr::select(ENSG_ID, ENST_ID, gene_name)

# table to convert peak IDs and ENSG_IDs
peak_conversion_table <- read.delim("./04_bedtools_merge_mutant-wildtype/2021.04/2021.04.05_mutant-wildtype_unique-shared_peaks.txt", header = TRUE, sep = "\t")%>%
  dplyr::mutate(unique_ID = case_when(
    cell1 == 'WT' & gtf_score == 1 ~"Wildtype",
    cell1 == 'M' & gtf_score ==1 ~ "Mutant",
    gtf_score == 2 ~ "Both")) %>%
  dplyr::select(chr, start, end, ENST1, strand1, attribute, unique_ID) %>%
  dplyr::rename(peak_ID = attribute, strand = strand1, ENSGID = ENST1) %>%
  separate(col = ENSGID, into = c("ENSG_ID", NA), sep = "[.]") # Need to include brackets around period or it means 'any character'

# Importing the gene level expression analysis
# This is mapped against the gencodeV35 GTF and is in ENSG namespace
gene_level <- read.delim(file = "./06_differential_peak_analysis/2021.03/2021.03.09_gene-level-analysis_MutvsWT_deseq2-results.csv", header = TRUE, sep = ",") %>%
  dplyr::rename(ENSGID = X, gene_l2fc = log2FoldChange, gene_padj = padj) %>%
  separate(col = ENSGID, into = c("ENSG_ID", NA), sep = "[.]") 

# importing the data for the different analyses--DESEQ2
multifactor2 <- read.delim(file = "./06_differential_peak_analysis/2021.04.05/2021.04.06_multifactor_2count-filter_counts_above-0_Mut-vs-WT_deseq2_results.csv", header = TRUE, sep = ",") %>%
  dplyr::rename(peak_ID = X, peak_l2fc = log2FoldChange, DEseq_padj = padj)

#Importing the QNB data
# The genes supplied to QNB were "pre-filtered" with DESeq2 before we ran the analysis
QNB_filter <- read.csv(file ="./06_differential_peak_analysis/2021.04.05/2021.04.06_QNB-results_multifactor-prefiltered_2counts-peaklist.csv", header = TRUE) %>%
  dplyr::rename(peak_ID = X, QNB_padj = padj) %>%
  dplyr::select(peak_ID, QNB_padj)



#### 2. Combine data tables ####
# Goal = Making a data table to compare the gene level and peak level data

# merging peak conversion and name conversion
master_conversion <- inner_join(name_conversion, peak_conversion_table, by = "ENSG_ID") %>%
  dplyr::select(ENSG_ID, gene_name, peak_ID, chr, start, end, strand, unique_ID )


# Annotating the gene_level table (ENSG names) with the conversion table
gene_level_conversion_table <- merge(master_conversion, gene_level, by="ENSG_ID") %>%
  dplyr::select(ENSG_ID, gene_name, peak_ID, chr, start, end, strand, unique_ID, gene_l2fc, gene_padj )

# Joining the gene_level table with our peak_l2fc table
counts2_table <- dplyr::inner_join(gene_level_conversion_table, multifactor2, by ="peak_ID") %>%
  dplyr::select(ENSG_ID, gene_name, peak_ID, chr, start, end, strand, unique_ID, gene_l2fc, gene_padj, peak_l2fc, DEseq_padj) %>%
  dplyr::mutate (diff_l2fc = peak_l2fc - gene_l2fc)


## Adding the QNB data to the data table
GTF_count_table <- inner_join(counts2_table, QNB_filter)


#### Adding the log10 (padj) column and annotating the column based on l2fc ####
# Preparing the data to be plotted
master_volcano <- GTF_count_table %>%
  arrange() %>%
  mutate("logDEseqP" = -1*log10(DEseq_padj)) %>%
  mutate("logQNBP" = -1*log10(QNB_padj)) %>%
  mutate("padj_sig" = case_when(
    (DEseq_padj > 0.05 & QNB_padj > 0.05) ~"not_sig",
    (DEseq_padj < 0.05 & QNB_padj > 0.05) ~ "DESeq",
    (DEseq_padj > 0.05 & QNB_padj < 0.05) ~ "QNB",
    (DEseq_padj < 0.05 & QNB_padj < 0.05) ~"Both", 
    TRUE ~ "other")) %>%
  mutate(direction = case_when(
      diff_l2fc <= -1.5 & DEseq_padj < 0.05 ~ "diff_l2fc < -1.5",
      diff_l2fc >= 1.5 & DEseq_padj < 0.05 ~ "diff_l2fc > 1.5",
      TRUE ~ "no_change")) %>%
  dplyr::mutate(direction = factor(direction, levels = c("no_change","diff_l2fc < -1.5", "diff_l2fc > 1.5"))) 


#distinct number of genes in this dataframe
volcano_stats <-master_volcano %>%
  dplyr::filter(direction == 'diff_l2fc > 1.5')

#### 3. Generating plots #### 
# Making the volcano plot
plot2 <- ggplot(master_volcano, aes(x = diff_l2fc, y = logDEseqP, color = direction, label = gene_name)) +
  geom_vline(xintercept = -1.5, linetype="dashed", color = "red") + geom_vline(xintercept = 1.5, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype="dashed", color = "black") +
  geom_hline(yintercept = 1.30, linetype = "dashed", color = "red") +
  geom_point() +
  theme_bw() +
  xlim(-5,5) +
  theme(legend.position = "none") +
  scale_color_manual(values=c("black", "black", "black")) +
  labs(x = "Log2FoldChange IP (Mut/WT)", y = "-Log10(DEseq2 padj)", 
       title = "Differential m6A-expression in UOK262 cells") 
plot2

plot3 <- ggExtra::ggMarginal(plot2, type = c("histogram"), margins = 'x', xparams = list( bins = 50))
plot3

# plotly interactive version of the above volcano
test1 <- ggplotly(plot2)
test1


# Volcano plot colored by density
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

master_volcano_omit <- na.omit(master_volcano)

master_volcano_omit$density <- get_density(master_volcano_omit$diff_l2fc, master_volcano_omit$logDEseqP, n = 100)

plot_den <- ggplot(master_volcano_omit, aes(x = diff_l2fc, y = logDEseqP, color = density)) +
  geom_point() +
  theme_bw() +
  xlim(-5,5) +
  geom_vline(xintercept = -1.5, linetype="dashed") + geom_vline(xintercept = 1.5, linetype = "dashed") +
  geom_hline(yintercept = 1.30, linetype = "dashed") +
  scale_color_viridis()+
  labs(x = "Log2FoldChange IP (Mut/WT)", y = "-Log10(DEseq2 padj)", 
       title = "Differential m6A-expression in UOK262 cells") 
plot_den


#### 4. How to the QNB and DEseq values compare to each other?####
scatter_padj <- ggplot(master_volcano, aes(x = DEseq_padj, y = QNB_padj, color = padj_sig)) +
  geom_point() +
  theme_bw()+
  labs(title = "Scatter plot of QNB and DEseq2 padj values")+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")
scatter_padj


padj_table <- master_volcano %>%
  dplyr::filter(padj_sig != 'not_sig') %>%
  na.omit

quantile(padj_table$QNB_padj)
quantile(padj_table$DEseq_padj)

scatter_padj <- ggplot(padj_table, aes(x = DEseq_padj, y = QNB_padj, color = padj_sig)) +
  geom_point() +
  theme_bw()+
  labs(title = "Scatter plot of QNB and DEseq2 padj values", subtitle = "Peaks with padj > 0.05 in both removed")+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0.05) + geom_hline(yintercept = 0.05)
scatter_padj


##### 5. Subset data for HOMER motif analysis####

HOMER_all_peaks <- master_volcano %>% dplyr::select(peak_ID, chr, start, end, strand)
#write_delim(HOMER_all_peaks, "./09_homer_analysis/2021.04.08_ALL-peaks.txt", delim = "\t", col_names = FALSE)

WT_unique <- master_volcano %>% dplyr::filter(unique_ID == "Wildtype") %>% dplyr::select(peak_ID, chr, start, end, strand)
#write_delim(WT_unique, "./09_homer_analysis/2021.04.08_WT-unique-peaks.txt", delim = "\t", col_names = FALSE)

Mutant_unique <- master_volcano %>% dplyr::filter(unique_ID == "Mutant") %>% dplyr::select(peak_ID, chr, start, end, strand)
#write_delim(Mutant_unique, "./09_homer_analysis/2021.04.08_Mutant-unique-peaks.txt", delim = "\t", col_names = FALSE)

Both_celltype <- master_volcano %>% dplyr::filter(unique_ID == "Both") %>% dplyr::select(peak_ID, chr, start, end, strand)
#write_delim(Both_celltype, "./09_homer_analysis/2021.04.08_both-celltype-peaks.txt", delim = "\t", col_names = FALSE)

upregulated_peaks <- master_volcano %>% dplyr::filter (direction == "diff_l2fc > 1.5") %>% dplyr::select(ENSG_ID, chr, start, end, strand, diff_l2fc)
#write_delim(upregulated_peaks, "./10_GO.Analysis/2021.04.08_upregulated-peaks.txt", delim = "\t", col_names = TRUE)

downregulated_peaks <- master_volcano %>% dplyr::filter (direction == "diff_l2fc < -1.5") %>% dplyr::select(ENSG_ID, peak_ID, chr, start, end, strand, diff_l2fc)
#write_delim(downregulated_peaks, "./10_GO.Analysis/2021.04.08_downregulated-peaks.txt", delim = "\t", col_names = TRUE)


upregulated_mutant_other <- master_volcano_fomit %>% 
  dplyr::filter (direction == "diff_l2fc > 1.5") %>%
  dplyr::filter(unique_ID == 'Mutant') %>%
  dplyr::filter(peak_location == 'other_location') %>%
  dplyr::select(ENSG_ID, chr, start, end, strand, diff_l2fc)
#write_delim(upregulated_mutant_other, "./10_GO.Analysis/2021.05.04_upregulated-peaks_mutantpeaks_location-not-first100.txt", delim = "\t", col_names = TRUE)

upregulated_mutant_first100 <- master_volcano_fomit %>% 
  dplyr::filter (direction == "diff_l2fc > 1.5") %>%
  dplyr::filter(unique_ID == 'Mutant') %>%
  dplyr::filter(peak_location == 'first_100') %>%
  dplyr::select(ENSG_ID, chr, start, end, strand, diff_l2fc)
#write_delim(upregulated_mutant_first100, "./10_GO.Analysis/2021.05.04_upregulated-peaks_mutantpeaks_location-first100.txt", delim = "\t", col_names = TRUE)


upregulated_both_other <- master_volcano_fomit %>% 
  dplyr::filter (direction == "diff_l2fc > 1.5") %>%
  dplyr::filter(unique_ID == 'Both') %>%
  dplyr::filter(peak_location == 'other_location') %>%
  dplyr::select(ENSG_ID, chr, start, end, strand, diff_l2fc)
#write_delim(upregulated_both_other, "./10_GO.Analysis/2021.05.04_upregulated-peaks_bothcells_location-not-first100.txt", delim = "\t", col_names = TRUE)

upregulated_both_first100 <- master_volcano_fomit %>% 
  dplyr::filter (direction == "diff_l2fc > 1.5") %>%
  dplyr::filter(unique_ID == 'Both') %>%
  dplyr::filter(peak_location == 'first_100') %>%
  dplyr::select(ENSG_ID, chr, start, end, strand, diff_l2fc)
#write_delim(upregulated_both_first100, "./10_GO.Analysis/2021.05.04_upregulated-peaks_bothcells_location-first100.txt", delim = "\t", col_names = TRUE)


downreg_WT_other <- master_volcano_fomit %>% 
  dplyr::filter (direction == "diff_l2fc < -1.5") %>%
  dplyr::filter(unique_ID == 'Wildtype') %>%
  dplyr::filter(peak_location == 'other_location') %>%
  dplyr::select(ENSG_ID, chr, start, end, strand, diff_l2fc)
write_delim(upregulated_both_other, "./10_GO.Analysis/2021.05.04_downreg-peaks_WT_location-not-first100.txt", delim = "\t", col_names = TRUE)

downreg_WT_first100 <- master_volcano_fomit %>% 
  dplyr::filter (direction == "diff_l2fc < -1.5") %>%
  dplyr::filter(unique_ID == 'Wildtype') %>%
  dplyr::filter(peak_location == 'first_100') %>%
  dplyr::select(ENSG_ID, chr, start, end, strand, diff_l2fc)
write_delim(upregulated_both_first100, "./10_GO.Analysis/2021.05.09_dwonreg-peaks_WT_location-first100.txt", delim = "\t", col_names = TRUE)





#### 6. Bringing in the Jaffrey halflife data ####
Jaffrey_data <- read.delim(file = "/Users/fitzsimmonscm/Documents/Projects_and_Data/Batista_Lab/FH_SDHB_Project/Rprojects/UOK_manuscript/SLAMseq/collapsed_reads/2021.04.08_Jaffrey_6hr_12hr_halflife.txt",
                           header=TRUE, sep = "\t") %>%
  dplyr::select(gene_name, variant, half6, half12) %>%
  dplyr::rename(ENSG_ID = gene_name) %>%
  pivot_wider(names_from = variant, values_from = c(half6, half12)) %>%
  dplyr::mutate(half6_ratio = half6_Mutant/half6_Wildtype, half12_ratio = half12_Mutant/half12_Wildtype) 


NLS_data <- read.delim(file = "/Users/fitzsimmonscm/Documents/Projects_and_Data/Batista_Lab/FH_SDHB_Project/Rprojects/UOK_manuscript/SLAMseq/2021.04.20_simple-NLS-halflife.txt",
                       header = TRUE,
                       sep = "\t") %>%
  dplyr::rename(ENSG_ID = gene_name, NLS_Mhalf = Mutant, NLS_WThalf = Wildtype)


master_volcano2 <- left_join(master_volcano, Jaffrey_data)
master_volcano3 <- left_join(master_volcano2, NLS_data, by=)
#write_delim(master_volcano2, file = "./06_differential_peak_analysis/2021.04.08_TABLE-S3_preliminary_m6A_master_table.txt", delim = "\t", col_names = TRUE)


#### 7. Annotating the peaks with HOMER motif data ####

motif_annotation <- read.delim(file = "./09_homer_analysis/2021.04/len5_rna_plus-bkdg/all-peaks/peak_annotation_motif1/2021.04.09_all-peaks_RRACH_annotation_locations.txt", 
                               header = TRUE, sep = "\t") %>%
  dplyr::rename(peak_ID = contains(match = "annotatePeaks.pl"), annotate_info = contains(match = "X1"), start = Start, end = End, chr = Chr, strand = Strand) %>%
  dplyr::select(chr, start, end, peak_ID, strand, annotate_info) 


##### Analysis of possible m6Am peaks--annotated the first 100 nt ####
first100_fomit <- read.delim(file = "./08_first100nt_intersection/2021.04.12_gencodeV35_first100nt_intersection_fomit_V2.txt", 
                             header = FALSE, sep = "\t") %>%
  dplyr::filter(V13 > 1) %>%
  dplyr::select(V1, V2, V3, V4, V6, V13) %>%
  dplyr::rename (chr = V1, start = V2, end = V3, peak_ID = V4, strand = V6, overlap = V13)

master_volcano_fomit_semi <- semi_join(master_volcano3, first100_fomit, by = c("peak_ID")) %>% dplyr::mutate(peak_location = 'first_100')
master_volcano_fomit_anti <- anti_join(master_volcano3, first100_fomit, by = c("peak_ID")) %>% dplyr::mutate(peak_location = 'other_location')
master_volcano_fomit <- bind_rows(master_volcano_fomit_semi, master_volcano_fomit_anti) %>%
  dplyr::mutate(peak_location = factor(peak_location, levels = c("first_100", "other_location")))


# Ploting the data from the f25 bedtools intersection
plot3 <- ggplot() +
  geom_point(data = master_volcano_fomit_anti, aes(x = diff_l2fc, y = logDEseqP, color = peak_location)) +
  geom_point(data = master_volcano_fomit_semi, aes(x = diff_l2fc, y = logDEseqP, color = peak_location)) +
  theme_bw() +
  xlim(-5,5) +
  geom_vline(xintercept = -1.5, linetype="dashed") + geom_vline(xintercept = 1.5, linetype = "dashed") +
  geom_hline(yintercept = 1.30, linetype = "dashed") +
  scale_color_manual(values=c("purple", "grey")) +
  labs(x = "Log2FoldChange IP (Mut/WT)", y = "-Log10(DEseq2_padj)", 
       title = "Differential m6A-expression in UOK262 cells", 
       subtitle = "Possible m6Am peaks: f = omit") 
plot3


#### 8. Comparing our list of first 100 nt peaks to PCIF1 papers from Jaffrey and Shi labs ####

jaffrey_PCIF1 <- readxl::read_xlsx(path = "./08_first100nt_intersection/Jaffrey_m6Am_peaklist.xlsx", sheet = 1, col_names = TRUE) %>%
  dplyr::filter(confident_gene_assignment != ".") %>%
  dplyr::rename(gene_name = confident_gene_assignment)%>%
  dplyr::distinct()

Shi_PCIF1 <- readxl::read_xlsx(path = "./08_first100nt_intersection/Shi_m6Am_peaklist.xlsx", sheet = 1, col_names = TRUE, skip = 1) %>%
  dplyr::rename(ENSG_ID = 'Ensembl Gene ID', gene_name = Name) %>%
  dplyr::select(ENSG_ID, gene_name) %>%
  dplyr::distinct()

# what common targets are shared by the Jaffrey and Shi lists?
PCIF1_shared <- dplyr::inner_join(jaffrey_PCIF1, Shi_PCIF1, by = 'gene_name') %>% dplyr::mutate(PCIF_target = 'Both') # 657 shared m6Am genes
PCIF1_jaffrey_unique <- dplyr::anti_join(jaffrey_PCIF1, Shi_PCIF1, by = 'gene_name') %>% dplyr::mutate(PCIF_target = 'Jaffrey') #1261
PCIF1_Shi_unique <- dplyr::anti_join(Shi_PCIF1, jaffrey_PCIF1, by = 'gene_name') %>% dplyr::mutate(PCIF_target = 'Shi') #886 


PCIF_masterlist <- bind_rows(PCIF1_shared, PCIF1_jaffrey_unique, PCIF1_Shi_unique)  %>%
  dplyr::select(gene_name, PCIF_target) #2804 genes on the combined list

batista_m6Am_targets <- master_volcano_fomit %>%
  dplyr::filter(peak_location == 'first_100') %>%
  dplyr::inner_join(PCIF_masterlist, by = 'gene_name')








