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
plot_pca_for_condition <- function(umi_counts, meta_data_subset,
                                   plot_title, legend_title, plot_file_name){
  pca_df <- umi_counts[meta_data_subset$sample_long_name,]
  pca_matrix <- prcomp(pca_df)
  
  fviz_pca_ind(pca_matrix, axes = c(1,2),
               habillage = meta_data_subset$condition,
               geom = c("point", "text"),
               pointsize = 2,
               addEllipses = T, ellipse.level = 0.8,
               alpha.ind = 1,
               label="none",
               legend.title = legend_title, 
               title = plot_title,
               ggtheme = theme_classic())
  # 
  # fviz_pca_var(pca_matrix) 
  ggsave(plot_file_name)  
}

options(scipen = 999) #removes scientific notation
options(scipen = 0) #resets scientific notation

#pca plot quant batch
meta_data_subset <- meta_data %>%
  select(sample_long_name, quant_batch) %>%
  rename(c("condition" = "quant_batch")) %>%
  arrange(condition)
plot_title <- "Variation with quantification batch"
legend_title <- "quantification batch"
plot_file_name <- "plots/pca_quant_batch.png"
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title, legend_title, plot_file_name)


#pca plot quality
meta_data_subset <- meta_data %>%
  select(sample_long_name, seq_miR_library_quality) %>%
  rename(c("condition" = "seq_miR_library_quality")) %>%
  arrange(condition)
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "Variation with quality", 
                       legend_title = "library quality", 
                       plot_file_name = "plots/pca_lib_quality.png")


#pca plot condition
meta_data_subset <- meta_data %>%
  filter(condition != "CF_pre_post_modulator") %>%
  select(sample_long_name, condition) %>%
  arrange(condition)
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "Variation with condition", 
                       legend_title = "condition", 
                       plot_file_name = "plots/pca_condition.png")


#pca plot CFVsHC
meta_data_subset <- meta_data %>%
  filter(condition != "CF_pre_post_modulator") %>%
  select(sample_long_name, condition) %>%
  mutate(condition = case_when(condition == "HC" ~ condition,
                               TRUE ~ "CF")) %>%
  arrange(condition)
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "CF Vs HC", 
                       legend_title = "condition", 
                       plot_file_name = "plots/pca_CFVsHC.png")


#pca plot CFVsHC good quality
meta_data_subset <- meta_data %>%
  filter(condition != "CF_pre_post_modulator", seq_miR_library_quality == "Good") %>%
  select(sample_long_name, condition) %>%
  mutate(condition = case_when(condition == "HC" ~ condition,
                               TRUE ~ "CF")) %>%
  arrange(condition)
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "CF Vs HC good quality samples", 
                       legend_title = "condition", 
                       plot_file_name = "plots/pca_CFVsHC_good.png")

#pca plot CFVsHC good quality AU
meta_data_subset <- meta_data %>%
  filter(condition != "CF_pre_post_modulator", 
         seq_miR_library_quality == "Good",
         country == "AU") %>%
  select(sample_long_name, condition) %>%
  mutate(condition = case_when(condition == "HC" ~ condition,
                               TRUE ~ "CF")) %>%
  arrange(condition)
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "CF Vs HC good quality samples from AU", 
                       legend_title = "condition", 
                       plot_file_name = "plots/pca_CFVsHC_good_AU.png")

#pca plot CFVsHC good quality DK
meta_data_subset <- meta_data %>%
  filter(condition != "CF_pre_post_modulator", 
         seq_miR_library_quality == "Good",
         country == "DK") %>%
  select(sample_long_name, condition) %>%
  mutate(condition = case_when(condition == "HC" ~ condition,
                               TRUE ~ "CF")) %>%
  arrange(condition)
#only CF samples, so cant find pca


#pca plot CFVsHC good quality AU adult
meta_data_subset <- meta_data %>%
  filter(condition != "CF_pre_post_modulator", 
         seq_miR_library_quality == "Good",
         country == "AU",
         age_group == "adult") %>%
  select(sample_long_name, condition) %>%
  mutate(condition = case_when(condition == "HC" ~ condition,
                               TRUE ~ "CF")) %>%
  arrange(condition)
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "CF Vs HC good quality adult samples from AU", 
                       legend_title = "condition", 
                       plot_file_name = "plots/pca_CFVsHC_good_AU_adult.png")

#pca plot CFVsHC good quality AU child
meta_data_subset <- meta_data %>%
  filter(condition != "CF_pre_post_modulator", 
         seq_miR_library_quality == "Good",
         country == "AU",
         age_group == "child") %>%
  select(sample_long_name, condition) %>%
  mutate(condition = case_when(condition == "HC" ~ condition,
                               TRUE ~ "CF")) %>%
  arrange(condition)
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "CF Vs HC good quality child samples from AU", 
                       legend_title = "condition", 
                       plot_file_name = "plots/pca_CFVsHC_good_AU_child.png")


#pca plot CFRDVsIGTVsNGT
meta_data_subset <- meta_data %>%
  filter(condition %in% c("CFRD", "IGT", "NGT")) %>%
  select(sample_long_name, condition) %>%
  arrange(condition)
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "CF Vs IGT Vs NGT", 
                       legend_title = "condition", 
                       plot_file_name = "plots/pca_CFRDVsIGTVsNGT.png")


#pca plot CFRDVsIGTVsNGT good quality samples
meta_data_subset <- meta_data %>%
  filter(condition %in% c("CFRD", "IGT", "NGT"), seq_miR_library_quality == "Good") %>%
  select(sample_long_name, condition) %>%
  arrange(condition)
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "CF Vs IGT Vs NGT good quality samples", 
                       legend_title = "condition", 
                       plot_file_name = "plots/pca_CFRDVsIGTVsNGT_good.png")


#pca plot CFRDVsIGTVsNGT good quality AU samples
meta_data_subset <- meta_data %>%
  filter(condition %in% c("CFRD", "IGT", "NGT"), 
         seq_miR_library_quality == "Good",
         country == "AU") %>%
  select(sample_long_name, condition) %>%
  arrange(condition)
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "CF Vs IGT Vs NGT good quality samples from AU", 
                       legend_title = "condition", 
                       plot_file_name = "plots/pca_CFRDVsIGTVsNGT_good_AU.png")

#pca plot CFRDVsIGTVsNGT good quality DK samples
meta_data_subset <- meta_data %>%
  filter(condition %in% c("CFRD", "IGT", "NGT"), 
         seq_miR_library_quality == "Good",
         country == "DK") %>%
  select(sample_long_name, condition) %>%
  arrange(condition)
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "CF Vs IGT Vs NGT good quality samples from DK", 
                       legend_title = "condition", 
                       plot_file_name = "plots/pca_CFRDVsIGTVsNGT_good_DK.png")


#pca plot  good quality AU samples
meta_data_subset <- meta_data %>%
  filter(condition %in% c("CFRD", "IGT", "NGT"), 
         seq_miR_library_quality == "Good",
         country == "AU") %>%
  select(sample_long_name, condition) %>%
  arrange(condition)
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "CF Vs IGT Vs NGT good quality samples from AU", 
                       legend_title = "condition", 
                       plot_file_name = "plots/pca_CFRDVsIGTVsNGT_good_AU.png")



#pca plot effect of modulator good quality samples
meta_data_subset <- meta_data %>%
  filter(!is.na(pre_post_modulator), seq_miR_library_quality == "Good") %>%
  select(sample_long_name, pre_post_modulator) %>%
  rename(c("condition" = "pre_post_modulator")) %>%
  arrange(condition)
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "pre post modulator good quality samples", 
                       legend_title = "pre post modulator", 
                       plot_file_name = "plots/pca_prepostmodulator_good.png")


#pca plot effect of modulator just for samples meant for pre_post_modulator comparison good quality samples
meta_data_subset <- meta_data %>%
  filter(!is.na(pre_post_modulator), 
         seq_miR_library_quality == "Good",
         condition == "CF_pre_post_modulator") %>%
  select(sample_long_name, pre_post_modulator) %>%
  rename(c("condition" = "pre_post_modulator")) %>%
  arrange(condition)
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "pre post modulator good quality samples specific", 
                       legend_title = "pre post modulator", 
                       plot_file_name = "plots/pca_prepostmodulator_good_specific.png")


#pca plot effect of modulator on good quality adult samples
meta_data_subset <- meta_data %>%
  filter(!is.na(pre_post_modulator), 
         seq_miR_library_quality == "Good",
         age_group == "adult") %>%
  select(sample_long_name, pre_post_modulator) %>%
  rename(c("condition" = "pre_post_modulator")) %>%
  arrange(condition)
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "pre post modulator good quality adult samples", 
                       legend_title = "pre post modulator", 
                       plot_file_name = "plots/pca_prepostmodulator_good_adult.png")


#pca plot effect of modulator on good quality adult AU samples
meta_data_subset <- meta_data %>%
  filter(!is.na(pre_post_modulator), 
         seq_miR_library_quality == "Good",
         age_group == "adult",
         country == "AU") %>%
  select(sample_long_name, pre_post_modulator) %>%
  rename(c("condition" = "pre_post_modulator")) %>%
  arrange(condition)
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "pre post modulator good quality adult AU samples", 
                       legend_title = "pre post modulator", 
                       plot_file_name = "plots/pca_prepostmodulator_good_adult_AU.png")


#pca plot effect of modulator on adult AU samples
meta_data_subset <- meta_data %>%
  filter(!is.na(pre_post_modulator), 
         age_group == "adult",
         country == "AU") %>%
  select(sample_long_name, pre_post_modulator) %>%
  rename(c("condition" = "pre_post_modulator")) %>%
  arrange(condition)
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "pre post modulator adult AU samples", 
                       legend_title = "pre post modulator", 
                       plot_file_name = "plots/pca_prepostmodulator_adult_AU.png")

#so selecting good quality samples only doesn't seem to make significant difference
