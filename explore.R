library(tidyverse)
library(factoextra)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/CysticFibrosisGroup/ExoCF/CFRD_EV_biomarker/"
setwd(base_dir)


#read transcriptomics data (miRNA)
umi_counts <- read.csv("data/formatted/umi_counts.csv", row.names = 1)
colnames(umi_counts) <- gsub("^X", "", colnames(umi_counts))

umi_counts <- data.frame(t(umi_counts))

#read metadata
meta_data <- read.csv("data/formatted/meta_data.csv")



#create pca plots

#create metadata subset and then merge with miRNA data
meta_data_subset <- meta_data %>%
  select(sample_long_name, quant_batch) %>%
  rename(c("condition" = "quant_batch")) %>%
  arrange(condition)
plot_title <- "Variation with quantification batch"
legend_title <- "quantification batch"
plot_file_name <- "plots/pca_quant_batch.png"


pca_df <- umi_counts[meta_data_subset$sample_long_name,]
pca_matrix <- prcomp(pca_df)

fviz_pca_ind(pca_matrix, axes = c(1,2),
             habillage = meta_data_subset$condition,
             geom = c("point", "text"),
             pointsize = 4,
             addEllipses = T, ellipse.level = 0.8,
             alpha.ind = 1,
             label="none",
             legend.title = legend_title, 
             title = plot_title)
# 
# fviz_pca_var(pca_matrix) 
ggsave(plot_file_name)
