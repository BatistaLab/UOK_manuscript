library(tidyverse)
library(ggpubr)
library(reshape2)
library(forcats)

utr_rates <- read.delim(file = "/Users/fitzsimmonscm/Documents/Projects_and_Data/Batista_Lab/FH_SDHB_Project/Rprojects/UOK_manuscript/SLAMseq/collapsed_reads/slamdunk_utrratesplot.tsv", header = TRUE, sep = "\t")

# Full plot of all timepoints
utr_melt <- melt(data=utr_rates, id=c("Category")) %>%
  dplyr::rename(replicate = Category, conversion = variable, frequency = value) %>%
  dplyr::mutate(replicate = factor(replicate, levels = c(
    "WT_12hrC", "WT_12hrB","WT_12hrA",
    "WT_06hrC", "WT_06hrB","WT_06hrA",
    "WT_0hrC", "WT_0hrB", "WT_0hrA",
    "WT_omitC", "WT_omitB", "WT_omitA",
    "Mut_12hrC", "Mut_12hrB", "Mut_12hrA",
    "Mut_06hrC", "Mut_06hrB", "Mut_06hrA", 
    "Mut_0hrC", "Mut_0hrB", "Mut_0hrA", 
    "Mut_omitC", "Mut_omitB", "Mut_omitA")))

utr_rate_plot <- ggplot(utr_melt, aes(x=replicate, y=frequency, fill = conversion)) +
  geom_bar(position ="stack", stat = "identity") +
  theme_bw() +
  coord_flip() +
  labs(title = "Mutation rate per 3'UTR base", x = "Replicate", y = "Conversion Rate")

utr_rate_plot


# Simplified plot of just T to C conversions
utr_simplified <- utr_rates %>%
  dplyr::select(1:2) %>%
  separate(Category, into = c("celltype", "replicates"), sep = "_") %>%
  dplyr::mutate(timepoint = c("0","0","0","12","12","12","6","6","6","omit","omit","omit",
                              "0","0","0","12","12","12","6","6","6","omit","omit","omit"), 
                replicate = as.numeric(c("1","2","3","1","2","3",
                                         "1","2","3","1","2","3",
                                         "1","2","3","1","2","3",
                                         "1","2","3","1","2","3"))) %>%
  dplyr::select(celltype, timepoint, replicate, T.C) %>%
  dplyr::rename(conversion_rate = T.C) %>%
  dplyr::mutate(timepoint = factor(timepoint, levels=c("omit", "0", "6", "12"))) %>%
  dplyr::mutate(celltype = factor(celltype, levels = c("Mut", "WT")))

simple_utr_rate_plot <- ggplot() +
  geom_bar(data = utr_simplified, aes(x=timepoint, y=conversion_rate, fill = celltype),
           show.legend = FALSE, stat = "summary", fun = "mean", position = "dodge" ) +
  geom_jitter(data = utr_simplified, aes(x=timepoint, y=conversion_rate), 
              show.legend = FALSE, stat="identity", position ="dodge") +
  theme_bw() +
  facet_grid(~celltype) +
  labs(title = "Mutation rate per 3'UTR base", y= "T.C Conversion Rate", x= "Timepoint (hr)")

simple_utr_rate_plot


