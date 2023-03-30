
#Set the directory for where your results (summary.txt) are from ImageJ
dir <- "/Users/lungerjc/Documents/Phenotype Assays/Invasion Assays/Transwell Invasion/2_23_2021_invasion_10k_48hr/invaded_cells/Cropped/Results/"

#Set the sequence of sample information that you are going to use to label each well
seq <- c("Cell", "Serum", "Replicate")

#Read summary.txt and add the sample information in the order by which imageJ read your samples (look at summary.txt if you are unsure of the order)
#Name your samples using information in the order of "seq", separating each piece of information with an underscore 
#For example, my seq is "Cell", "Serum", "Replicate", so my first well is A1, which is "WT_+_1"
summary <- read_tsv(file = paste(dir, "summary.txt", sep = "")) %>%
  mutate(sample = c("WT_+_1", "WT_+_2", "WT_+_3", "WT_+_4","WT_-_1","MUT_+_1", "MUT_+_2", "MUT_+_3", "MUT_+_4","MUT_-_1","PAR_+_1", "PAR_+_2", "PAR_+_3", "PAR_+_4","PAR_-_1","HK2_+_1", "HK2_+_2", "HK2_+_3", "HK2_+_4","HK2_-_1" )) %>%
  separate(`sample`, into = seq, sep = "_",remove = FALSE) #This step separates sample information into individual columns

#Order your cell types in the order you want them graphed
summary$Cell <- factor(summary$Cell, levels = c("HK2", "PAR", "WT", "MUT"))

#Melt the dataframe for ggplot
summary2 <- summary %>%
  melt(id = c("Cell", "Serum", "Replicate", "Count") )

#Plot each point individually
ggplot(summary, aes(x = Cell, y = Count,  color = Cell, group = Serum)) +
  geom_point(position=position_dodge(width=0.5))+
  theme_minimal()+
  scale_color_manual(values = c("lightgrey", "darkorange","gray", "olivedrab3")) +
  labs(title = "Transwell Invasion", y = "Count")

#Plot boxplots
ggplot(summary, aes(x= Serum, y = Count,fill=Cell)) +
  geom_boxplot()+
  theme_minimal() +
  scale_fill_manual(values = c("lightgrey", "darkorange","gray", "olivedrab3")) +
  labs(title = "Transwell Invasion", y = "Count")