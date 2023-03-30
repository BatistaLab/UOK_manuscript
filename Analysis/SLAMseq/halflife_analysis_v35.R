# Analysis of RNA half-life using first-order decay model
# Christina Fitzsimmons
# Last updated 2022-05-05

# Note: halflife analysis will crash if run in .Rmd format. Instead, use .R scripts from the command line. 

#### Read in the data ####
# Set working directory 

setwd ("/Users/fitzsimmonscm/Documents/Projects_and_Data/Batista_Lab/FH_SDHB_Project/Rprojects/UOK_manuscript/SLAMseq/")
base_input <- "./collapsed_reads/raw/"
base_output <- "./collapsed_reads/analysis/2022.03/"

# Libraries
library (tidyverse)
library (aomisc)
library (broom)
library (drc)
library (ggpubr)
library (matrixTests)
library (matrixStats)
library (minpack.lm)
library (plotly)
library (purrr)
library (tibble)


##### DEFINING THE IMPORT FUNCTION ####
read_NGS_output_file <- function(filepath) {
  data <- read_tsv(filepath, col_types = list(
    gene_name = col_character(),
    length = col_integer(),
    readsCPM = col_double(), 
    conversionRate = col_double(),    
    Tcontent = col_double(),
    coverageOnTs = col_double(),
    conversionsOnTs = col_double(),
    readCount = col_double(),
    tcReadCount = col_double(),
    multimapCount = col_double()
  ))
  return(data)
}

theme_set(theme_bw())
##### Import Mutant collapsed files ####
MA_omit <- read_NGS_output_file(filepath = "./collapsed_reads/raw/Mut_omitA_S1_L008_R1_001.fastq_slamdunk_mapped_filtered_tcount_collapsed.csv") %>%
  dplyr::mutate ('variant' = 'Mutant', replicate = 'A', timepoint = 24 , short='MAomit') %>%
  dplyr::select(gene_name, conversionRate, conversionsOnTs, readCount, variant, replicate, timepoint, short) 

MB_omit <- read_NGS_output_file(filepath = "./collapsed_reads/raw/Mut_omitB_S2_L008_R1_001.fastq_slamdunk_mapped_filtered_tcount_collapsed.csv") %>%
  dplyr::mutate ('variant' = 'Mutant', replicate = 'B', timepoint = 24 , short='MBomit') %>%
  dplyr::select(gene_name, conversionRate, conversionsOnTs, readCount, variant, replicate, timepoint, short) 

MC_omit <- read_NGS_output_file(filepath = "./collapsed_reads/raw/Mut_omitC_S3_L008_R1_001.fastq_slamdunk_mapped_filtered_tcount_collapsed.csv") %>%
  dplyr::mutate ('variant' = 'Mutant', replicate = 'C', timepoint = 24 , short='MComit') %>%
  dplyr::select(gene_name, conversionRate, conversionsOnTs, readCount, variant, replicate, timepoint, short)  


MA_0hr <- read_NGS_output_file(filepath = "./collapsed_reads/raw/Mut_0hrA_S4_L008_R1_001.fastq_slamdunk_mapped_filtered_tcount_collapsed.csv") %>% 
  dplyr::mutate ('variant' = 'Mutant', replicate = 'A', timepoint = 0, short='MA0') %>%
  dplyr::select(gene_name, conversionRate, conversionsOnTs, readCount, variant, replicate, timepoint, short)

MB_0hr <- read_NGS_output_file(filepath = "./collapsed_reads/raw/Mut_0hrB_S5_L008_R1_001.fastq_slamdunk_mapped_filtered_tcount_collapsed.csv") %>%
  dplyr::mutate ('variant' = 'Mutant', replicate = 'B', timepoint = 0, short='MB0') %>%
  dplyr::select(gene_name, conversionRate, conversionsOnTs, readCount, variant, replicate, timepoint, short)

MC_0hr <- read_NGS_output_file(filepath = "./collapsed_reads/raw/Mut_0hrC_S6_L008_R1_001.fastq_slamdunk_mapped_filtered_tcount_collapsed.csv") %>%
  dplyr::mutate ('variant' = 'Mutant', replicate = 'C', timepoint = 0, short = 'MC0') %>%
  dplyr::select(gene_name, conversionRate, conversionsOnTs, readCount, variant, replicate, timepoint, short)


MA_6hr <- read_NGS_output_file(filepath = "./collapsed_reads/raw/Mut_6hrA_S7_L008_R1_001.fastq_slamdunk_mapped_filtered_tcount_collapsed.csv") %>%
  dplyr::mutate ('variant' = 'Mutant', replicate = 'A', timepoint = 6, short='MA6') %>%
  dplyr::select(gene_name, conversionRate, conversionsOnTs, readCount, variant, replicate, timepoint, short)

MB_6hr <- read_NGS_output_file(filepath = "./collapsed_reads/raw/Mut_6hrB_S8_L008_R1_001.fastq_slamdunk_mapped_filtered_tcount_collapsed.csv") %>%
  dplyr::mutate ('variant' = 'Mutant', replicate = 'B', timepoint = 6, short = 'MB6') %>%
  dplyr::select(gene_name, conversionRate, conversionsOnTs, readCount, variant, replicate, timepoint, short)

MC_6hr <- read_NGS_output_file(filepath = "./collapsed_reads/raw/Mut_6hrC_S9_L008_R1_001.fastq_slamdunk_mapped_filtered_tcount_collapsed.csv") %>%
  dplyr::mutate ('variant' = 'Mutant', replicate = 'C', timepoint = 6, short = 'MC6') %>%
  dplyr::select(gene_name, conversionRate, conversionsOnTs, readCount, variant, replicate, timepoint, short)


MA_12hr <- read_NGS_output_file(filepath = "./collapsed_reads/raw/Mut_12hrA_S10_L008_R1_001.fastq_slamdunk_mapped_filtered_tcount_collapsed.csv") %>%
  dplyr::mutate ('variant' = 'Mutant', replicate = 'A', timepoint = 12, short ='MA12') %>%
  dplyr::select(gene_name, conversionRate, conversionsOnTs, readCount, variant, replicate, timepoint, short)

MB_12hr <- read_NGS_output_file(filepath = "./collapsed_reads/raw/Mut_12hrB_S11_L008_R1_001.fastq_slamdunk_mapped_filtered_tcount_collapsed.csv") %>%
  dplyr::mutate ('variant' = 'Mutant', replicate = 'B', timepoint = 12, short = 'MB12') %>%
  dplyr::select(gene_name, conversionRate, conversionsOnTs, readCount, variant, replicate, timepoint, short)

MC_12hr <- read_NGS_output_file(filepath = "./collapsed_reads/raw/Mut_12hrC_S12_L008_R1_001.fastq_slamdunk_mapped_filtered_tcount_collapsed.csv") %>%
  dplyr::mutate ('variant' = 'Mutant', replicate = 'C', timepoint = 12, short = 'MC12') %>%
  dplyr::select(gene_name, conversionRate, conversionsOnTs, readCount, variant, replicate, timepoint, short)



##### Import Wildtype Collapsed Files ####
WA_omit <- read_NGS_output_file(filepath = "./collapsed_reads/raw/WT_omitA_S13_L008_R1_001.fastq_slamdunk_mapped_filtered_tcount_collapsed.csv") %>%
  dplyr::mutate ('variant' = 'Wildtype', replicate = 'A', timepoint = 24 , short='WAomit') %>%
  dplyr::select(gene_name, conversionRate, conversionsOnTs, readCount, variant, replicate, timepoint, short) 

WB_omit <- read_NGS_output_file(filepath = "./collapsed_reads/raw/WT_omitB_S14_L008_R1_001.fastq_slamdunk_mapped_filtered_tcount_collapsed.csv") %>%
  dplyr::mutate ('variant' = 'Wildtype', replicate = 'B', timepoint = 24 , short='WBomit') %>%
  dplyr::select(gene_name, conversionRate, conversionsOnTs, readCount, variant, replicate, timepoint, short) 

WC_omit <- read_NGS_output_file(filepath = "./collapsed_reads/raw/WT_omitC_S15_L008_R1_001.fastq_slamdunk_mapped_filtered_tcount_collapsed.csv") %>%
  dplyr::mutate ('variant' = 'Wildtype', replicate = 'C', timepoint = 24 , short='WComit') %>%
  dplyr::select(gene_name, conversionRate, conversionsOnTs, readCount, variant, replicate, timepoint, short) 



WA_0hr <- read_NGS_output_file(filepath = "./collapsed_reads/raw/WT_0hrA_S16_L008_R1_001.fastq_slamdunk_mapped_filtered_tcount_collapsed.csv") %>%
  dplyr::mutate ('variant' = 'Wildtype', replicate = 'A', timepoint = 0, short ='WA0') %>%
  dplyr::select(gene_name, conversionRate, conversionsOnTs, readCount, variant, replicate, timepoint, short)
  
WB_0hr <- read_NGS_output_file(filepath = "./collapsed_reads/raw/WT_0hrB_S17_L008_R1_001.fastq_slamdunk_mapped_filtered_tcount_collapsed.csv") %>%
  dplyr::mutate ('variant' = 'Wildtype', replicate = 'B', timepoint = 0, short = 'WB0') %>%
  dplyr::select(gene_name, conversionRate, conversionsOnTs, readCount, variant, replicate, timepoint, short)

WC_0hr <- read_NGS_output_file(filepath = "./collapsed_reads/raw/WT_0hrC_S18_L008_R1_001.fastq_slamdunk_mapped_filtered_tcount_collapsed.csv") %>%
  dplyr::mutate ('variant' = 'Wildtype', replicate = 'C', timepoint = 0, short = 'WC0') %>%
  dplyr::select(gene_name, conversionRate, conversionsOnTs, readCount, variant, replicate, timepoint, short)


WA_6hr <- read_NGS_output_file(filepath = "./collapsed_reads/raw/WT_6hrA_S19_L008_R1_001.fastq_slamdunk_mapped_filtered_tcount_collapsed.csv") %>%
  dplyr::mutate ('variant' = 'Wildtype', replicate = 'A', timepoint = 6, short = 'WA6') %>%
  dplyr::select(gene_name, conversionRate, conversionsOnTs, readCount, variant, replicate, timepoint, short)

WB_6hr <- read_NGS_output_file(filepath = "./collapsed_reads/raw/WT_6hrB_S20_L008_R1_001.fastq_slamdunk_mapped_filtered_tcount_collapsed.csv") %>%
  dplyr::mutate ('variant' = 'Wildtype', replicate = 'B', timepoint = 6, short = 'WB6') %>%
  dplyr::select(gene_name, conversionRate, conversionsOnTs, readCount, variant, replicate, timepoint, short)

WC_6hr <- read_NGS_output_file(filepath = "./collapsed_reads/raw/WT_6hrC_S21_L008_R1_001.fastq_slamdunk_mapped_filtered_tcount_collapsed.csv") %>%
  dplyr::mutate ('variant' = 'Wildtype', replicate = 'C', timepoint = 6, short = 'WC6') %>%
  dplyr::select(gene_name, conversionRate, conversionsOnTs, readCount, variant, replicate, timepoint, short)


WA_12hr <- read_NGS_output_file(filepath = "./collapsed_reads/raw/WT_12hrA_S22_L008_R1_001.fastq_slamdunk_mapped_filtered_tcount_collapsed.csv") %>%
  dplyr::mutate ('variant' = 'Wildtype', replicate = 'A', timepoint = 12, short = 'WA12') %>%
  dplyr::select(gene_name, conversionRate, conversionsOnTs, readCount, variant, replicate, timepoint, short)

WB_12hr <- read_NGS_output_file(filepath = "./collapsed_reads/raw/WT_12hrB_S23_L008_R1_001.fastq_slamdunk_mapped_filtered_tcount_collapsed.csv") %>%
  dplyr::mutate ('variant' = 'Wildtype', replicate = 'B', timepoint = 12, short = 'WB12') %>%
  dplyr::select(gene_name, conversionRate, conversionsOnTs, readCount, variant, replicate, timepoint, short)

WC_12hr <- read_NGS_output_file(filepath = "./collapsed_reads/raw/WT_12hrC_S24_L008_R1_001.fastq_slamdunk_mapped_filtered_tcount_collapsed.csv") %>%
  dplyr::mutate ('variant' = 'Wildtype', replicate = 'C', timepoint = 12, short = 'WC12') %>%
  dplyr::select(gene_name, conversionRate, conversionsOnTs, readCount, variant, replicate, timepoint, short)

pre_filtered <- bind_rows(MA_omit, MA_0hr, MA_6hr, MA_12hr, 
                          MB_omit, MB_0hr, MB_6hr, MB_12hr, 
                          MC_omit, MC_0hr, MC_6hr, MC_12hr, 
                          WA_omit, WA_0hr, WA_6hr, WA_12hr, 
                          WB_omit, WB_0hr, WB_6hr, WB_12hr, 
                          WC_omit, WC_0hr, WC_6hr, WC_12hr)

#write out pre-filtered data set
write_delim(pre_filtered, file = "./collapsed_reads/analysis/2022.03/2022.03.23_gencodeV35_prefiltered_data.txt", delim = "\t", col_names = TRUE)

##### Preliminary filtering for genes that are seen in all timepoints ####

#How many genes do we start with?
unique_genes <- pre_filtered %>%
  group_by(gene_name) %>%
  dplyr::distinct(gene_name)

fulldata <- pre_filtered %>%
  dplyr::add_count(gene_name, name = 'replicate_count') %>%
  dplyr::filter (replicate_count == 24) %>%
  ungroup()

#### Background subtraction of omit timepoints ####
bkgd_corrected <- fulldata %>%
  dplyr::select(gene_name, variant, timepoint, replicate, conversionRate) %>%
  pivot_wider(names_from = c(timepoint),
               values_from = conversionRate) %>%
  dplyr::mutate(bkgd_0 = `0`-`24`, bkgd_6 = `6`-`24`, bkgd_12 = `12`-`24`) %>% # generate background-subtracted columns
  dplyr::mutate(bkgd_0 = if_else(bkgd_0 < 0, 0, bkgd_0),
                bkgd_6 = if_else(bkgd_6 < 0, 0, bkgd_6),
                bkgd_12 = if_else(bkgd_12 < 0, 0, bkgd_12)) %>% # get rid of negative counts. Can only have positive counts
  dplyr::filter(bkgd_0 >= bkgd_6 & bkgd_6 >= bkgd_12) %>% # filter so that conversion values are always decreasing
  dplyr::mutate(row_sum = bkgd_0 + bkgd_6 + bkgd_12) %>% 
  dplyr::filter(row_sum != 0 ) %>% # remove rows with zeros across all time points. 
  dplyr::mutate(norm_0 = bkgd_0/bkgd_0, norm_6 = bkgd_6/bkgd_0, norm_12 = bkgd_12/bkgd_0) %>%
  dplyr::select(gene_name, variant, replicate, norm_0, norm_6, norm_12) %>% # select only columns of interest
  dplyr::rename(`0` = norm_0, `6` = norm_6, `12` = norm_12) %>% #convert 'bkgd' back to numbers only for conversion to long format. 
  pivot_longer(cols = c(`0`, `6`, `12`), names_to = 'timepoint', 
               values_to = "conversionRate", names_transform = list(timepoint = as.double))

unique_genes <- bkgd_corrected %>%
  group_by(gene_name) %>%
  dplyr::distinct(gene_name)

# write out background corrected data set
write_delim(bkgd_corrected, file = "./collapsed_reads/analysis/2022.03/2022.05.02_gencodeV35_background_subtracted_data.txt", delim = "\t", col_names = TRUE)

#### dataframe to test map function ####
test_replicate <- c("A", "A", "A", "B", "B", "B", "C", "C", "C")
test_timepoint <- c(0, 6, 12, 0, 6, 12, 0, 6, 12)
test_conversion <- c(29, 24, 13, 62, 34, 27, 89, 15, 5)

test_df <- data.frame(replicate, timepoint, conversion)

test_fit <- nlsLM(conversion ~ Plat*(exp(-k * timepoint)), 
                  data = test_df,
                  start = list(Plat = , k = .693),
                  upper = c(Inf, Inf), lower = c(0, 0),
                  control = nls.lm.control(maxiter = 1000))

plot1 <- ggplot(test_df, aes(x = timepoint, y = conversion)) +
  geom_point()
plot1


#### 1. Fitting a 1-phase decay function with NLS and some constraints ####

# fit 1 = orginal slam-seq paper: https://github.com/breichholf/smRNAseq/blob/master/analysis/final.R
# fit 2 = self-starting fxn recommended here: https://www.statforbiology.com/2020/stat_nls_usefulfunctions/#exponential-function 
poss_nlsLM <- possibly(.f = nlsLM, otherwise = NULL) # wrapper for NLS function. If NLS raises error, it will print NULL and continue
poss_NLS_SS <- possibly(.f = nls, otherwise = NULL)


fitted <- bkgd_corrected %>%
  nest(data = c(-gene_name, -variant)) %>%
  dplyr::mutate(
    fit1 = (map(data, ~poss_nlsLM(conversionRate ~ Plat * (exp(-k * timepoint)),
                                data = .,
                                start = list(Plat = 0.05, k = 0.5),
                                upper = c(Inf, Inf), lower = c(0, 0),
                                control = nls.lm.control(maxiter = 1000))))) %>%
  dplyr::mutate(
    fit2 =map(data, ~poss_NLS_SS(conversionRate~NLS.expoDecay(timepoint, a, k), data = .)))


#### Removing NULL values and selecting which function fit better ####
fit.info <-fitted %>% 
  dplyr::mutate(fit_select = case_when(
    is.null(fit1) & is.null(fit2) ~ "none", 
    is.null(fit1) ~ "NLS_SS",
    is.null(fit2) ~ "NLS_LM", 
    TRUE ~ "unknown"))
  dplyr::filter(fit_select != "none")

#### unpacking fit 1 to plot some basic summary statistics ####
unpack1 <- fitted %>%
  dplyr::mutate(tidied = map(fit1, tidy), 
                augmented = map(fit1, augment), 
                glanced = map(fit1, glance)) %>%
  unnest(data) %>%
  unnest(tidied) %>% 
  unnest(augmented, names_repair = "unique") %>% 
  unnest(glanced)



simple_unpack1 <- unpack1 %>%
  dplyr::select(gene_name, replicate, variant, term, estimate, std.error, statistic, p.value, deviance) %>%
  dplyr::filter(term =='k') %>%
  dplyr::distinct() %>%
  mutate(half1 = 0.693/(estimate)) %>%
  dplyr::filter (half1 > 0.035) %>% # removes rows with mostly 0/1 very low counts but that were not removed above
  group_by(gene_name, variant) %>%
  dplyr::add_count(gene_name, name="replicate_counts") %>%
  dplyr::filter(replicate_counts > 1) # makes sure we have at least 2 replicates for each timepoint that we are considering. 


write_csv(simple_unpack1, file = "./2022.05.20_SLAMseq_halflife_model_fit1.csv")

#### Unpacking fit2 and writing out ####

unpack2 <- fitted %>%
  dplyr::mutate(tidied = map(fit2, tidy), 
                augmented = map(fit2, augment), 
                glanced = map(fit2, glance)) %>%
  unnest(tidied) %>% 
  unnest(augmented) %>% 
  unnest(glanced)


simple_unpack2 <- unpack2 %>%
  dplyr::select(gene_name, variant, term, estimate, std.error, statistic, p.value, .fitted, .resid, deviance) %>%
  dplyr::filter(term =='k') %>%
  dplyr::distinct() %>%
  mutate(half2 = 0.693/(estimate)) %>%
  group_by(gene_name, variant)


write_csv(simple_unpack2, file = "./2022.05.05_SLAMseq_halflife_model_fit2.csv")