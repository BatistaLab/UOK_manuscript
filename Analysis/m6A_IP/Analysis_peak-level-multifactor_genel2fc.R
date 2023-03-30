# Analysis of differential gene expression at the peak IP level with gene-level correction
# Christina Fitzsimmons
# Last updated 2021-06-10

#### Read in the data ####
# Set working directory 
setwd ("/Users/fitzsimmonscm/Documents/Projects_and_Data/Batista_Lab/FH_SDHB_Project/Rprojects/UOK_manuscript/m6A_IP/")
library(DESeq2)
library(ggpubr)
library(MASS)
library(plotly)
library(tidyverse)
library(viridis)

#importing the data for conversions
#importing the data for conversions
name_conversion <- read.csv(file = "/Users/fitzsimmonscm/Documents/Projects_and_Data/Batista_Lab/FH_SDHB_Project/Rprojects/UOK_manuscript/longest_isoform/gProfiler_hsapiens_2021_03_10_conversion.csv", header = TRUE, sep = ",")%>%
  dplyr::select(ENSG_ID, ENST_ID, gene_name)

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
QNB_filter <- read.csv(file ="./06_differential_peak_analysis/2021.04.05/2021.04.06_QNB-results_multifactor-prefiltered_2counts-peaklist.csv", header = TRUE) %>%
  dplyr::rename(peak_ID = X, QNB_padj = padj) %>%
  dplyr::select(peak_ID, QNB_padj)



#### Combine data tables ####
# Making a data table to compare the gene level and peak level data

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


## Adding the QNB data
GTF_count_table <- inner_join(counts2_table, QNB_filter)


### Adding the log10 (padj) column and annotating the column based on l2fc
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


# Making the volcano plot
plot2 <- ggplot(master_volcano, aes(x = diff_l2fc, y = logDEseqP, color = direction, label = gene_name)) +
  geom_point() +
  theme_bw() +
  xlim(-5,5) +
  geom_vline(xintercept = -1.5, linetype="dashed") + geom_vline(xintercept = 1.5, linetype = "dashed") +
  geom_hline(yintercept = 1.30, linetype = "dashed") +
  scale_color_manual(values=c("black", "black", "black")) +
  labs(x = "Log2FoldChange IP (Mut/WT)", y = "-Log10(DEseq2 padj)", 
       title = "Differential m6A-expression in UOK262 cells") 
plot2


# plotly version of above volcano
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


#### How to the QNB and DEseq values compare to each other?####
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


##### Subset data for HOMER ####

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




#### Bringing in the Jaffrey halflife data ####
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


halflife_plots <- master_volcano3 %>%
  na.omit() %>%
  group_by(diff_l2fc, half12_ratio) %>%
  dplyr::arrange(desc(diff_l2fc, half12_ratio))

half_volcano <- halflife_plots %>%
  dplyr::filter (direction == "diff_l2fc > 1.5") %>%
  dplyr::filter (half12_Mutant <= 48 & half12_Wildtype <= 48)


plot4 <- ggplot(data = half_volcano, aes(x = half12_Wildtype, y = half12_Mutant, label = gene_name)) +
  geom_point() +
  theme_bw() +
  labs(title = 'Scatterplot of halflife', subtitle = 'Peaks with padj < 0.05 and diff_l2fc > 1.5') +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed')
plot4

test4 <- ggplotly(plot4)
test4


#### Annotating the peaks with HOMER motif data ####

motif_annotation <- read.delim(file = "./09_homer_analysis/2021.04/len5_rna_plus-bkdg/all-peaks/peak_annotation_motif1/2021.04.09_all-peaks_RRACH_annotation_locations.txt", 
                               header = TRUE, sep = "\t") %>%
  dplyr::rename(peak_ID = contains(match = "annotatePeaks.pl"), annotate_info = contains(match = "X1"), start = Start, end = End, chr = Chr, strand = Strand) %>%
  dplyr::select(chr, start, end, peak_ID, strand, annotate_info) 


##### Analysis of possible m6Am peaks
first100_fomit <- read.delim(file = "./08_first100nt_intersection/2021.04.12_gencodeV35_first100nt_intersection_fomit_V2.txt", 
                             header = FALSE, sep = "\t") %>%
  dplyr::filter(V13 > 1) %>%
  dplyr::select(V1, V2, V3, V4, V6, V13) %>%
  dplyr::rename (chr = V1, start = V2, end = V3, peak_ID = V4, strand = V6, overlap = V13)

master_volcano_fomit_semi <- semi_join(master_volcano3, first100_fomit, by = c("peak_ID")) %>% dplyr::mutate(peak_location = 'first_100')
master_volcano_fomit_anti <- anti_join(master_volcano3, first100_fomit, by = c("peak_ID")) %>% dplyr::mutate(peak_location = 'other_location')
master_volcano_fomit <- bind_rows(master_volcano_fomit_semi, master_volcano_fomit_anti) %>%
  dplyr::mutate(peak_location = factor(peak_location, levels = c("first_100", "other_location")))



# f25 plot
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


#### Comparing our list of first 100 nt peaks to PCIF1 papers from Jaffrey and Shi labs ####

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


plot5 <- ggplot(batista_m6Am_targets, aes(x = PCIF_target, y = half12_ratio)) +
  geom_violin(aes(color = PCIF_target)) +
  theme_bw() 
plot5

plot6 <- ggplot() +
  geom_point(data = batista_m6Am_targets, aes(x = diff_l2fc, y = logDEseqP, color = PCIF_target)) +
  theme_bw() +
  xlim(-5,5) +
  geom_vline(xintercept = -1.5, linetype="dashed") + geom_vline(xintercept = 1.5, linetype = "dashed") +
  geom_hline(yintercept = 1.30, linetype = "dashed") +
  labs(x = "Log2FoldChange IP (Mut/WT)", y = "-Log10(DEseq2_padj)", 
       title = "Differential m6A-expression in UOK262 cells", 
       subtitle = "Possible m6Am peaks colored by publication") 
plot6






#### Importing EMT gene signatures ####
EMT_sig <- readxl::read_xlsx(path = "./11_EMT_intersection/EMT_gene_signature_dataset2.xlsx", col_names = TRUE) %>%
  dplyr::rename(gene_name = Gene)

# how many EMT genes are on our sig upregulated list?
EMT_m6A_up <- master_volcano_fomit %>%
  dplyr::filter(direction == 'diff_l2fc > 1.5') %>%
  inner_join(EMT_sig, by = c('gene_name'))


#### RRACH motif annotation ####
RRACH_annotation <- read.delim(file = './09_homer_analysis/2021.04/len5_rna_plus-bkdg/all-peaks/peak_annotation_motif1/2021.04.09_all-peaks_RRACH_annotation_locations.txt',
                               header = TRUE, 
                               sep = '\t') %>%
  dplyr::rename(peak_ID = 1, RRACH_motif = 22) %>%
  dplyr::select(peak_ID, RRACH_motif)


master_volcano_export <- left_join(master_volcano_fomit, RRACH_annotation, by="peak_ID")
#write_delim(master_volcano_export, file = "./2021.04.21_m6A_master_peaklist.txt", col_names = TRUE, delim = "\t")


#### Comparing half-life of peaks between mutant and WT ####
half_analysis <- master_volcano_export  %>%
  na.omit() %>%
  dplyr::filter(peak_location == 'other_location')

scatter12 <- ggplot(half_analysis, aes(x = half12_Wildtype, y = NLS_WThalf)) +
  geom_point()+
  theme_bw() +
  ylim(0, 50) +
  xlim(0, 50) +
  geom_abline(slope = 1, intercept = 0) +
  labs(title = 'Scatterplot of Half-life calculations', x = 'Jaffrey 12 hr half-life', y = 'NLS Wildtype half-life')
scatter12


# peaks that are unique to mutant only
stat1 <- master_volcano_fomit %>%
  dplyr::filter(unique_ID == 'Mutant') %>%
  dplyr::filter(peak_location == 'first_100') %>%
  dplyr::filter(direction == 'diff_l2fc > 1.5')


stat2 <- master_volcano_fomit %>%
  dplyr::filter(unique_ID == 'Mutant') %>%
  dplyr::filter(peak_location == 'first_100') %>%
  dplyr::filter(direction == 'diff_l2fc < -1.5')


stat3 <- master_volcano_fomit %>%
  dplyr::filter(unique_ID == 'Mutant') %>%
  dplyr::filter(peak_location == 'other_location') %>%
  dplyr::filter(direction == 'diff_l2fc > 1.5')


stat4 <- master_volcano_fomit %>%
  dplyr::filter(unique_ID == 'Mutant') %>%
  dplyr::filter(peak_location == 'other_location') %>%
  dplyr::filter(direction == 'diff_l2fc < -1.5')


# Peaks that are found in Both celltypes
stat5 <- master_volcano_fomit %>%
  dplyr::filter(unique_ID == 'Both') %>%
  dplyr::filter(peak_location == 'first_100') %>%
  dplyr::filter(direction == 'diff_l2fc > 1.5')


stat6 <- master_volcano_fomit %>%
  dplyr::filter(unique_ID == 'Both') %>%
  dplyr::filter(peak_location == 'first_100') %>%
  dplyr::filter(direction == 'diff_l2fc < -1.5')


stat7 <- master_volcano_fomit %>%
  dplyr::filter(unique_ID == 'Both') %>%
  dplyr::filter(peak_location == 'other_location') %>%
  dplyr::filter(direction == 'diff_l2fc > 1.5')


stat8 <- master_volcano_fomit %>%
  dplyr::filter(unique_ID == 'Both') %>%
  dplyr::filter(peak_location == 'other_location') %>%
  dplyr::filter(direction == 'diff_l2fc < -1.5')


# Peaks that are found in Wildtype only
stat9 <- master_volcano_fomit %>%
  dplyr::filter(unique_ID == 'Wildtype') %>%
  dplyr::filter(peak_location == 'first_100') %>%
  dplyr::filter(direction == 'diff_l2fc > 1.5')


stat10 <- master_volcano_fomit %>%
  dplyr::filter(unique_ID == 'Wildtype') %>%
  dplyr::filter(peak_location == 'first_100') %>%
  dplyr::filter(direction == 'diff_l2fc < -1.5')


stat11 <- master_volcano_fomit %>%
  dplyr::filter(unique_ID == 'Wildtype') %>%
  dplyr::filter(peak_location == 'other_location') %>%
  dplyr::filter(direction == 'diff_l2fc > 1.5')


stat12 <- master_volcano_fomit %>%
  dplyr::filter(unique_ID == 'Wildtype') %>%
  dplyr::filter(peak_location == 'other_location') %>%
  dplyr::filter(direction == 'diff_l2fc < -1.5')


#### Filtering joins to look for unique genes ####

join1.2 <- anti_join(stat1, stat2, by = "ENSG_ID")
join1.3 <- anti_join(stat1, stat3, by = "ENSG_ID")
join1.4 <- anti_join(stat1, stat4, by = "ENSG_ID")
join1.5 <- anti_join(stat1, stat5, by = "ENSG_ID") %>% dplyr::mutate(join = 'join1.5')
join1.6 <- anti_join(stat1, stat6, by = "ENSG_ID")


join5.1 <- anti_join(stat5, stat1, by = "ENSG_ID") %>% dplyr::mutate(join = 'join5.1')
join1.10 <- anti_join(stat1, stat10, by = "ENSG_ID") %>% dplyr::mutate(join = 'join1.10')
join10.1 <- anti_join(stat10, stat1, by = "ENSG_ID") %>% dplyr::mutate(join = 'join10.1')
join3.7 <- anti_join(stat3, stat7, by = "ENSG_ID") %>% dplyr::mutate (join = 'join3.7')
join7.3 <- anti_join(stat7, stat3, by = "ENSG_ID") %>% dplyr::mutate(join = 'join7.3')
join3.12 <- anti_join(stat3, stat12, by = "ENSG_ID") %>% dplyr::mutate(join = 'join3.12')
join12.3 <- anti_join(stat12, stat3, by = "ENSG_ID") %>% dplyr::mutate (join = 'join12.3')


df_antt_half <- bind_rows(join3.7, join3.12, join7.3, join12.3) %>%
  na.omit()




# boxplot of counts
half_plot_joined <- ggplot(df_antt_half, aes(x = join, y = half6_ratio, color = join)) +
  geom_boxplot()+
  theme_bw()
  labs( title= 'Boxplot of half-life ratio (6 hr)', x = 'unique genes', y = 'half-life ratio (Mut/WT)')
half_plot_joined

interact_box <- ggplotly(half_plot_joined)
interact_box


(data = half_analysis_significant, aes(x = peakcount_percentile, y = half12_ratio, color = peakcount_percentile, label = 'gene_name')) +
  geom_boxplot() +
  theme_bw() +
  geom_hline(yintercept = 1, linetype = 'dashed')+
  labs(title = 'Boxplot of half-life Ratio', x = 'Number of m6A peaks in CDS/3UTR', y = 'Ratio Mut/WT halflife')

half_plot6



Jaffrey_plot <- ggplot() +
  geom_boxplot(data = master_volcano_export, aes(x = half6_Mutant)) +
  geom_boxplot(data = master_volcano_export, aes(x = half6_Wildtype)) +
  geom_point()+
  theme_bw() 

Jaffrey_plot


# Histogram of counts
halfplot_1 <- ggplot(half_analysis, aes(x = peakcount_exon3UTR)) +
  geom_histogram() +
  theme_bw() +
  labs(title = 'Histogram showing number of peaks per gene', subtitle = 'Peak location: not located in first 100 nt')
halfplot_1

quantile(half_analysis$peakcount_exon3UTR)

# ECDF of counts
half_plot2 <- ggplot(data = half_analysis, aes(x = half12_ratio, color = peakcount_percentile)) +
  stat_ecdf() +
  theme_bw() +
  labs(title = 'ECDF plot of 12-hr half-life ratio', subtitle = 'Colored by peak percentile', y = 'cumulative fraction', x = ' Ratio Mut half / WT half')
half_plot2

Jaffrey_long <- read.delim(file = "/Users/fitzsimmonscm/Documents/Projects_and_Data/Batista_Lab/FH_SDHB_Project/Rprojects/UOK_manuscript/SLAMseq/collapsed_reads/2021.04.08_Jaffrey_6hr_12hr_halflife.txt",
                           header=TRUE, sep = "\t") %>%
  dplyr::select(gene_name, variant, half6, half12) %>%
  dplyr::rename(ENSG_ID = gene_name)

Jaffrey_long_m6AIntersect <- left_join(Jaffrey_long, master_volcano_fomit) %>%
  na.omit() %>%
  dplyr::filter(peak_location == 'other_location') %>%
  add_count(gene_name, name = 'peakcount_exon3UTR') %>%
  dplyr::mutate('peakcount_percentile' =
                  case_when(
                    peakcount_exon3UTR == 2 ~ 'peakcount=1',
                    peakcount_exon3UTR == 4 ~ 'peakcount=2',
                    peakcount_exon3UTR == 6 ~ 'peakcount=3',
                    peakcount_exon3UTR >= 8 ~ 'peakcount >= 4',
                    TRUE ~ 'other')) %>%
  dplyr::mutate(peakcount_percentile = factor(peakcount_percentile, levels = c("peakcount=1","peakcount=2", "peakcount=3", "peakcount >= 4")))
  
  
# ECDF of counts
half_plot3 <- ggplot(data = Jaffrey_long_m6AIntersect, aes(x = half12, color = peakcount_percentile)) +
  stat_ecdf() +
  theme_bw() +
  xlim(0,24) +
  facet_wrap(~variant, ncol = 4) +
  labs(title = 'ECDF plot of 12-hr half-life', x = 'half-life', y = 'cumulative fraction', subtitle = 'Facet by Celltype')
half_plot3

# boxplot of counts
half_plot4 <- ggplot(data = half_analysis, aes(x = peakcount_percentile, y = half12_ratio, color = peakcount_percentile)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = 'Boxplot of half-life Ratio', x = 'Number of m6A peaks in CDS/3UTR', y = 'Ratio Mut/WT halflife')
half_plot4


#### Filtering the half-life for diff l2fc ####

Jaffrey_long_significant <- Jaffrey_long_m6AIntersect %>%
  dplyr::filter(direction == 'diff_l2fc > 1.5')

half_analysis_significant <- half_analysis %>%
  dplyr::filter(direction == 'diff_l2fc > 1.5')

ks.test(half_analysis_significant$half12_Mutant, half_analysis_significant$half12_Wildtype, alternative = c('two.sided'))



half_plot5 <- ggplot(data = Jaffrey_long_significant, aes(x = half12, color = variant)) +
  stat_ecdf() +
  theme_bw() +
  xlim(0,24) +
  labs(title = 'ECDF plot of 12-hr half-life', x = 'half-life', y = 'cumulative fraction', subtitle = 'Filtered for padj and l2fc')
half_plot5

# boxplot of counts
half_plot6 <- ggplot(data = half_analysis_significant, aes(x = peakcount_percentile, y = half12_ratio, color = peakcount_percentile, label = 'gene_name')) +
  geom_boxplot() +
  theme_bw() +
  geom_hline(yintercept = 1, linetype = 'dashed')+
  labs(title = 'Boxplot of half-life Ratio', x = 'Number of m6A peaks in CDS/3UTR', y = 'Ratio Mut/WT halflife')
  
half_plot6


#### Generating a list for GO-analysis ####
GO_termlist <- half_analysis_significant %>% 
  dplyr::filter(half12_ratio < 1)

#write_delim(GO_termlist, file = "2021.04.18_mastertable_halflife_GOterms.txt", delim = "\t", col_names = TRUE)

GO_table <- read.delim(file = "/Users/fitzsimmonscm/Desktop/pantherChart.txt", sep = "\t", header = FALSE) %>%
  dplyr::rename(Category = V2, number_genes = V3, gene_percent_gene_total = V4, percent_gene_process_total = V5) %>%
  arrange(desc(gene_percent_gene_total))

GO_plot <- ggplot(GO_table, aes(x = reorder(Category, gene_percent_gene_total), y = gene_percent_gene_total)) +
  geom_col() +
  theme_bw() +
  coord_flip() +
  labs (title = 'GO-terms m6A peaks', subtitle = 'padj < 0.05, diff_l2fc >= 1.5, half_ratio < 1', y = 'Gene percent of total gene count', x = 'GO-term category')

GO_plot



#### Both cell types + RAS genes ####
Ras_GO <- read.delim(file = "./10_GO.Analysis/upregulated/Ras_signaling_pathway_celltype-both_p-0.05.csv", sep = ",", header = TRUE) %>%
  dplyr::rename(ENSG_ID = converted_alias)
Ras_filter_semi <- semi_join(master_volcano_fomit, Ras_GO) %>% dplyr::mutate(GO_term = 'Reg. of Ras')
Ras_filter_anti <- anti_join(master_volcano_fomit, Ras_GO) %>% dplyr::mutate(GO_term = 'other')

ACTIN_GO <- read.delim(file="./10_GO.Analysis/upregulated/gProfiler_hsapiens_5-5-2021_mutant_not100_upreg_actinGO.csv", sep = ',', header = TRUE) %>%
  dplyr::rename(ENSG_ID = converted_alias)

actin_filter_semi <- semi_join(master_volcano_fomit, ACTIN_GO ) %>% dplyr::mutate(GO_term = 'Cytoskeleton')
actin_filter_anti <- anti_join(master_volcano_fomit, ACTIN_GO) %>% dplyr::mutate(GO_term = 'other')



Ras_plot <- ggplot() +
  geom_point(data = actin_filter_anti, aes(x = diff_l2fc, y = logDEseqP, color = GO_term)) +
  geom_point(data = actin_filter_semi, aes(x = diff_l2fc, y = logDEseqP, color = GO_term)) +
  theme_bw() +
  xlim(-5,5) +
  geom_vline(xintercept = -1.5, linetype="dashed") + geom_vline(xintercept = 1.5, linetype = "dashed") +
  geom_hline(yintercept = 1.30, linetype = "dashed") +
  scale_color_manual(values=c("red", "grey")) +
  labs(x = "Log2FoldChange IP (Mut/WT)", y = "-Log10(DEseq2_padj)", 
       title = "Differential m6A-expression in UOK262 cells", 
       subtitle = "GO: Actin Cytoskeleton Organization") 
Ras_plot


#### ALL upregulated GO-term annotation ####
ALL_upregulated <- readxl::read_xlsx(path = "./10_GO.Analysis/upregulated/METASCAPE_ALL_UPREGULATED/2021.05.14_GO.BiologicalProcesses_Upregulated/metascape_result.xlsx", 
                                     col_names = TRUE) %>%
  dplyr::rename(ENSG_ID = MyList)

GO_celljxn <- ALL_upregulated %>% 
  filter_at(8, all_vars(. == '1.0'))


actin_filter_semi <- semi_join(master_volcano_fomit, GO_celljxn ) %>% dplyr::mutate(GO_term = 'cell_junction') %>%
  dplyr::filter(direction == 'diff_l2fc > 1.5')

## Count the number of unique reads on the master qPCR validation list
qPCR_list <- readxl::read_xlsx(path = "/Users/fitzsimmonscm/Desktop/2021.05.28_qPCR_master.xlsx", 
                               col_names = TRUE)

qPCR_list_count <- qPCR_list %>%
  dplyr::add_count(gene_name, name = 'list_count')





#### Mutant upregulated GO-terms ####
Mutup_GO <- readxl::read_excel(path ="./10_GO.Analysis/upregulated/mutantunique_upregulated_GOterm_gene_list.xlsx", 
                               col_names = TRUE, col_types = NULL) %>%
  dplyr::rename(ENSG_ID = MyList)


GO0030036 <- Mutup_GO %>% 
  filter_at(16, all_vars(. == '1.0'))

actin_filter_semi <- semi_join(master_volcano_fomit, GO0030036 ) %>% dplyr::mutate(GO_term = 'Actin Cytoskeleton') %>%
  dplyr::filter(direction == 'diff_l2fc > 1.5')


actin_filter_anti <- anti_join(master_volcano_fomit, GO0030036) %>% dplyr::mutate(GO_term = 'Other')


Actin_plot <- ggplot() +
  geom_point(data = actin_filter_anti, aes(x = diff_l2fc, y = logDEseqP, color = GO_term)) +
  geom_point(data = actin_filter_semi, aes(x = diff_l2fc, y = logDEseqP, color = GO_term)) +
  theme_bw() +
  xlim(-5,5) +
  geom_vline(xintercept = -1.5, linetype="dashed") + geom_vline(xintercept = 1.5, linetype = "dashed") +
  geom_hline(yintercept = 1.30, linetype = "dashed") +
  scale_color_manual(values=c("red", "grey")) +
  labs(x = "Log2FoldChange IP (Mut/WT)", y = "-Log10(DEseq2_padj)", 
       title = "Differential m6A-expression in UOK262 cells", 
       subtitle = "GO: Actin Cytoskeleton Organization") 
Actin_plot

# Boxplot to compare direction of change
p6 <- ggplot(actin_filter_semi, aes(x = join, y = aom_log2)) +
  geom_violin(aes(fill=join)) +
  theme_bw() +
  labs(x = "Unique peaks and direction of change", 
       y = "Log2 Half-life Ratio (Mut/WT)", 
       title = "Half-Life of m6A Transcripts not in first 100 nt") +
  ggpubr::rotate_x_text(angle = 45) +
  stat_compare_means(ref.group = "non-target", method = "t.test", aes(label=..p.adj..))
p6

actin_filter_vertical <- aomisc_vertical %>%
  dplyr::semi_join(nontarget5)


#### mRNA processing ####

GO0006397 <- Mutup_GO %>% 
  filter_at(17, all_vars(. == '1.0'))

RNAprocess_filter_semi <- semi_join(master_volcano_fomit, GO0006397  ) %>% dplyr::mutate(GO_term = 'mRNA Processing') %>%
  dplyr::filter(direction == 'diff_l2fc > 1.5')


RNAprocess_filter_anti <- anti_join(master_volcano_fomit, GO0006397 ) %>% dplyr::mutate(GO_term = 'Other')


RNAprocess_plot <- ggplot() +
  geom_point(data = RNAprocess_filter_anti, aes(x = diff_l2fc, y = logDEseqP, color = GO_term)) +
  geom_point(data = RNAprocess_filter_semi, aes(x = diff_l2fc, y = logDEseqP, color = GO_term)) +
  theme_bw() +
  xlim(-5,5) +
  geom_vline(xintercept = -1.5, linetype="dashed") + geom_vline(xintercept = 1.5, linetype = "dashed") +
  geom_hline(yintercept = 1.30, linetype = "dashed") +
  scale_color_manual(values=c("blue", "grey")) +
  labs(x = "Log2FoldChange IP (Mut/WT)", y = "-Log10(DEseq2_padj)", 
       title = "Differential m6A-expression in UOK262 cells", 
       subtitle = "GO: mRNA Processing") 
RNAprocess_plot













#### response to growth factor ####
WTdown_GO <- readxl::read_excel(path ="./10_GO.Analysis/WT_downreg_not100/metascape_result.xlsx", 
                               col_names = TRUE, col_types = NULL) %>%
  dplyr::rename(ENSG_ID = MyList)


GO0070848 <- WTdown_GO %>% 
  filter_at(7, all_vars(. == '1.0'))

growthfactor_filter_semi <- semi_join(master_volcano_fomit, GO0070848) %>% dplyr::mutate(GO_term = 'Response to Growth Factor')
growthfactor_filter_anti <- anti_join(master_volcano_fomit, GO0070848) %>% dplyr::mutate(GO_term = 'Other')


growthfactor_plot <- ggplot() +
  geom_point(data = growthfactor_filter_anti, aes(x = diff_l2fc, y = logDEseqP, color = GO_term)) +
  geom_point(data = growthfactor_filter_semi, aes(x = diff_l2fc, y = logDEseqP, color = GO_term)) +
  theme_bw() +
  xlim(-5,5) +
  geom_vline(xintercept = -1.5, linetype="dashed") + geom_vline(xintercept = 1.5, linetype = "dashed") +
  geom_hline(yintercept = 1.30, linetype = "dashed") +
  scale_color_manual(values=c("grey", "green")) +
  labs(x = "Log2FoldChange IP (Mut/WT)", y = "-Log10(DEseq2_padj)", 
       title = "Differential m6A-expression in UOK262 cells", 
       subtitle = "GO: Response to Growth Factor") 
growthfactor_plot 



WTdown_GO <- readxl::read_excel(path ="./10_GO.Analysis/WT_downreg_not100/metascape_result.xlsx", 
                                col_names = TRUE, col_types = NULL) %>%
  dplyr::rename(ENSG_ID = MyList)

#### blood vessel growth ####
GO0001568 <- WTdown_GO %>% 
  filter_at(8, all_vars(. == '1.0'))

blood_filter_semi <- semi_join(master_volcano_fomit, GO0001568) %>% dplyr::mutate(GO_term = 'Blood Vessel Development')
blood_filter_anti <- anti_join(master_volcano_fomit, GO0001568) %>% dplyr::mutate(GO_term = 'Other')


growthfactor_plot <- ggplot() +
  geom_point(data = blood_filter_anti, aes(x = diff_l2fc, y = logDEseqP, color = GO_term)) +
  geom_point(data = blood_filter_semi, aes(x = diff_l2fc, y = logDEseqP, color = GO_term)) +
  theme_bw() +
  xlim(-5,5) +
  geom_vline(xintercept = -1.5, linetype="dashed") + geom_vline(xintercept = 1.5, linetype = "dashed") +
  geom_hline(yintercept = 1.30, linetype = "dashed") +
  scale_color_manual(values=c("purple", "grey")) +
  labs(x = "Log2FoldChange IP (Mut/WT)", y = "-Log10(DEseq2_padj)", 
       title = "Differential m6A-expression in UOK262 cells", 
       subtitle = "GO: Blood Vessel Development") 
growthfactor_plot 

