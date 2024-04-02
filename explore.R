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

sum(is.na(meta_data$seq_plate))

#create pca plots

#create metadata subset and then merge with miRNA data
plot_pca_for_condition <- function(umi_counts, meta_data_subset,
                                   plot_title, legend_title, plot_file_name){
  pca_df <- umi_counts[meta_data_subset$sample_long_name,]
  pca_matrix <- prcomp(pca_df)
  
  count_df <- data.frame(count = summary(factor(meta_data_subset$condition))) %>%
    rownames_to_column("condition")

  meta_data_subset <- meta_data_subset %>%
    mutate(condition = factor(condition)) %>%
    inner_join(count_df) %>%
    mutate(condition = paste0(condition, " (", count, ")"))
  
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
  filter(is.na(pre_post_modulator) | pre_post_modulator != 1) %>%
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
  filter(is.na(pre_post_modulator) | pre_post_modulator != 1) %>%
  filter(seq_miR_library_quality == "Good") %>%
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
  filter(is.na(pre_post_modulator) | pre_post_modulator != 1) %>%
  filter(seq_miR_library_quality == "Good", country == "AU") %>%
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
  filter(is.na(pre_post_modulator) | pre_post_modulator != 1) %>%
  filter(seq_miR_library_quality == "Good", country == "DK") %>%
  select(sample_long_name, condition) %>%
  mutate(condition = case_when(condition == "HC" ~ condition,
                               TRUE ~ "CF")) %>%
  arrange(condition)
#only CF samples, so cant find pca


#pca plot CFVsHC good quality AU adult
meta_data_subset <- meta_data %>%
  filter(is.na(pre_post_modulator) | pre_post_modulator != 1) %>%
  filter(seq_miR_library_quality == "Good", country == "AU", age_group == "adult") %>%
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
  filter(is.na(pre_post_modulator) | pre_post_modulator != 1) %>%
  filter(seq_miR_library_quality == "Good", country == "AU", age_group == "child") %>%
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
  filter(is.na(pre_post_modulator) | pre_post_modulator != 1) %>%
  filter(condition %in% c("CFRD", "IGT", "NGT")) %>%
  select(sample_long_name, condition) %>%
  arrange(condition)
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "CFRD Vs IGT Vs NGT", 
                       legend_title = "condition", 
                       plot_file_name = "plots/pca_CFRDVsIGTVsNGT.png")


#pca plot CFRDVsIGTVsNGT good quality samples
meta_data_subset <- meta_data %>%
  filter(is.na(pre_post_modulator) | pre_post_modulator != 1) %>%
  filter(condition %in% c("CFRD", "IGT", "NGT"), seq_miR_library_quality == "Good") %>%
  select(sample_long_name, condition) %>%
  arrange(condition)
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "CFRD Vs IGT Vs NGT good quality samples", 
                       legend_title = "condition", 
                       plot_file_name = "plots/pca_CFRDVsIGTVsNGT_good.png")


#pca plot CFRDVsIGTVsNGT good quality AU samples
meta_data_subset <- meta_data %>%
  filter(is.na(pre_post_modulator) | pre_post_modulator != 1) %>%
  filter(condition %in% c("CFRD", "IGT", "NGT"), 
         seq_miR_library_quality == "Good",
         country == "AU") %>%
  select(sample_long_name, condition) %>%
  arrange(condition)
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "CFRD Vs IGT Vs NGT good quality samples from AU", 
                       legend_title = "condition", 
                       plot_file_name = "plots/pca_CFRDVsIGTVsNGT_good_AU.png")

#pca plot CFRDVsIGTVsNGT good quality DK samples
meta_data_subset <- meta_data %>%
  filter(is.na(pre_post_modulator) | pre_post_modulator != 1) %>%
  filter(condition %in% c("CFRD", "IGT", "NGT"), 
         seq_miR_library_quality == "Good",
         country == "DK") %>%
  select(sample_long_name, condition) %>%
  arrange(condition)
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "CFRD Vs IGT Vs NGT good quality samples from DK", 
                       legend_title = "condition", 
                       plot_file_name = "plots/pca_CFRDVsIGTVsNGT_good_DK.png")


#####effect of modulator

#pca plot effect of modulator
meta_data_subset <- meta_data %>%
  filter(!is.na(pre_post_modulator)) %>%
  select(sample_long_name, pre_post_modulator) %>%
  rename(c("condition" = "pre_post_modulator")) %>%
  mutate(condition = case_when(condition == 0 ~ "pre-modulator",
                                 condition == 1 ~ "post-modulator")) %>%
  mutate(condition = factor(condition, levels = c("pre-modulator", "post-modulator"))) %>%
  arrange(condition)
#note : specifying factor levels did not change the ordering of legend

plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "Pre Vs Post modulator samples", 
                       legend_title = "pre post modulator", 
                       plot_file_name = "plots/pca_prepostmodulator.png")


#pca plot effect of modulator good quality samples
meta_data_subset <- meta_data %>%
  filter(!is.na(pre_post_modulator), seq_miR_library_quality == "Good") %>%
  select(sample_long_name, pre_post_modulator) %>%
  rename(c("condition" = "pre_post_modulator")) %>%
  mutate("condition" = case_when(condition == 0 ~ "pre-modulator",
                                 condition == 1 ~ "post-modulator")) %>%
  arrange(condition)
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "Pre Vs Post modulator good quality samples", 
                       legend_title = "pre post modulator", 
                       plot_file_name = "plots/pca_prepostmodulator_good.png")


#pca plot effect of modulator just for samples meant for pre_post_modulator comparison good quality samples
meta_data_subset <- meta_data %>%
  filter(!is.na(pre_post_modulator), 
         seq_miR_library_quality == "Good",
         condition == "CF_pre_post_modulator") %>%
  select(sample_long_name, pre_post_modulator) %>%
  rename(c("condition" = "pre_post_modulator")) %>%
  mutate("condition" = case_when(condition == 0 ~ "pre-modulator",
                                 condition == 1 ~ "post-modulator")) %>%
  arrange(condition)
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "Pre Vs Post modulator good quality samples specific", 
                       legend_title = "pre post modulator", 
                       plot_file_name = "plots/pca_prepostmodulator_good_specific.png")


#pca plot effect of modulator on good quality adult samples
meta_data_subset <- meta_data %>%
  filter(!is.na(pre_post_modulator), 
         seq_miR_library_quality == "Good",
         age_group == "adult") %>%
  select(sample_long_name, pre_post_modulator) %>%
  rename(c("condition" = "pre_post_modulator")) %>%
  mutate("condition" = case_when(condition == 0 ~ "pre-modulator",
                                 condition == 1 ~ "post-modulator")) %>%
  arrange(condition)
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "Pre Vs Post modulator good quality adult samples", 
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
  mutate("condition" = case_when(condition == 0 ~ "pre-modulator",
                                 condition == 1 ~ "post-modulator")) %>%
  arrange(condition)
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "Pre Vs Post modulator good quality adult AU samples", 
                       legend_title = "pre post modulator", 
                       plot_file_name = "plots/pca_prepostmodulator_good_adult_AU.png")


#pca plot effect of modulator on adult AU samples
meta_data_subset <- meta_data %>%
  filter(!is.na(pre_post_modulator), 
         age_group == "adult",
         country == "AU") %>%
  select(sample_long_name, pre_post_modulator) %>%
  rename(c("condition" = "pre_post_modulator")) %>%
  mutate("condition" = case_when(condition == 0 ~ "pre-modulator",
                                 condition == 1 ~ "post-modulator")) %>%
  arrange(condition)
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "Pre Vs Post modulator adult AU samples", 
                       legend_title = "pre post modulator", 
                       plot_file_name = "plots/pca_prepostmodulator_adult_AU.png")

#so selecting good quality samples only doesn't seem to make significant difference


meta_data_subset <- meta_data %>%
  filter(!is.na(pre_post_modulator), 
         age_group == "child",
         country == "AU") %>%
  select(sample_long_name, pre_post_modulator) %>%
  rename(c("condition" = "pre_post_modulator")) %>%
  mutate("condition" = case_when(condition == 0 ~ "pre-modulator",
                                 condition == 1 ~ "post-modulator")) %>%
  arrange(condition)
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "Pre Vs Post modulator child AU samples", 
                       legend_title = "pre post modulator", 
                       plot_file_name = "plots/pca_prepostmodulator_child_AU.png")



#pre-post modulator for CFRD
meta_data_subset <- meta_data %>%
  filter(!is.na(pre_post_modulator)) %>%
  filter(condition == "CFRD") %>%
  select(sample_long_name, pre_post_modulator) %>%
  rename(c("condition" = "pre_post_modulator")) %>%
  mutate(condition = case_when(condition == 0 ~ "pre-modulator",
                               condition == 1 ~ "post-modulator")) %>%
  arrange(condition)

plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "Pre Vs Post modulator CFRD samples", 
                       legend_title = "pre post modulator", 
                       plot_file_name = "plots/pca_prepostmodulator_CFRD.png")



#pre-post modulator for IGT
meta_data_subset <- meta_data %>%
  filter(!is.na(pre_post_modulator)) %>%
  filter(condition == "IGT") %>%
  select(sample_long_name, pre_post_modulator) %>%
  rename(c("condition" = "pre_post_modulator")) %>%
  mutate(condition = case_when(condition == 0 ~ "pre-modulator",
                               condition == 1 ~ "post-modulator")) %>%
  arrange(condition)

plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "Pre Vs Post modulator IGT samples", 
                       legend_title = "pre post modulator", 
                       plot_file_name = "plots/pca_prepostmodulator_IGT.png")



#pre-post modulator for NGT
meta_data_subset <- meta_data %>%
  filter(!is.na(pre_post_modulator)) %>%
  filter(condition == "NGT") %>%
  select(sample_long_name, pre_post_modulator) %>%
  rename(c("condition" = "pre_post_modulator")) %>%
  mutate(condition = case_when(condition == 0 ~ "pre-modulator",
                               condition == 1 ~ "post-modulator")) %>%
  arrange(condition)

plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title = "Pre Vs Post modulator NGT samples", 
                       legend_title = "pre post modulator", 
                       plot_file_name = "plots/pca_prepostmodulator_NGT.png")




############# multiple factors in pca


plot_title <- "CFRD Vs IGT Vs NGT"
legend_title <- "condition"
plot_file_name <- "plots/pca_after_norm/pca_CFRDVsIGTVsNGT.png"
plot_pca_for_condition(umi_counts, meta_data_subset,
                       plot_title, legend_title, plot_file_name)




plot_pca_2factors <- function(umi_counts, meta_data_subset,
                              plot_title, legend_title, plot_file_name, norm){
  
  count_df <- data.frame(count = summary(factor(meta_data_subset$condition))) %>%
    rownames_to_column("condition")
  meta_data_subset <- meta_data_subset %>%
    mutate(condition = factor(condition)) %>%
    inner_join(count_df) %>%
    mutate(condition = paste0(condition, " (", count, ")")) %>%
    select(-c(count))
  
  count_df <- data.frame(count = summary(factor(meta_data_subset$age_group))) %>%
    rownames_to_column("age_group")
  meta_data_subset <- meta_data_subset %>%
    mutate(age_group = factor(age_group)) %>%
    inner_join(count_df) %>%
    mutate(age_group = paste0(age_group, " (", count, ")")) %>%
    select(-c(count))
  
  
  pca_df <- umi_counts[meta_data_subset$sample_long_name,]
  pca_df <- as.data.frame(t(pca_df))

  
  if(norm == "znorm_logCPM"){

    pca_df <- edgeR::cpm(pca_df, log=TRUE)
    pca_df <- as.data.frame(t(pca_df))
    
    pca_df <- scale(pca_df)
    # apply(pca_df, FUN = mean, MARGIN = 2)
    # apply(pca_df, FUN = sd, MARGIN = 2)
  } else if (norm == "znorm_logTMM"){
    
    group <- sapply(strsplit(meta_data_subset$condition, split = " "), 
                    FUN = function(x){
                      return (x[[1]])
                    })
    
    dge <- edgeR::DGEList(counts = pca_df, group = group)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    tmm <- edgeR::cpm(dge, log = TRUE)
    
    pca_df <- as.data.frame(t(tmm))
    pca_df <- scale(pca_df)

  } else{
    pca_df <- as.data.frame(t(pca_df))
  }
  
  
  pca_matrix <- prcomp(pca_df)
  
  # propor <- data.frame(summary(pca_matrix)$importance)
  # rowSums(propor)
  
  
  dim_red_df <- data.frame(pca_matrix$x[, c(1, 2)]) %>%
    rownames_to_column("sample_long_name") %>%
    inner_join(meta_data_subset)

  ggplot2::ggplot(dim_red_df,
                  ggplot2::aes(x = PC1, y = PC2, colour = condition, shape = age_group)) +
    ggplot2::geom_point() +
    ggplot2::labs(title = plot_title)

  ggsave(plot_file_name)  
}



#all

meta_data_subset <- meta_data %>%
  filter(is.na(pre_post_modulator) | pre_post_modulator != 1) %>%
  filter(condition %in% c("CFRD", "IGT", "NGT")) %>%
  select(sample_long_name, condition, age_group, sex, country) %>%
  arrange(condition)



plot_pca_2factors(umi_counts = umi_counts,
                  meta_data_subset = meta_data_subset,
                  plot_title = "CFRD Vs IGT Vs NGT znorm_logCPM",
                  legend_title = "condition",
                  norm = "znorm_logCPM",
                  plot_file_name = "plots/pca_after_norm/pca_CFRDVsIGTVsNGT_znorm_logCPM.png")
plot_pca_2factors(umi_counts = umi_counts,
                  meta_data_subset = meta_data_subset,
                  plot_title = "CFRD Vs IGT Vs NGT znorm_logTMM",
                  legend_title = "condition",
                  norm = "znorm_logTMM",
                  plot_file_name = "plots/pca_after_norm/pca_CFRDVsIGTVsNGT_znorm_logTMM.png")
plot_pca_2factors(umi_counts = umi_counts,
                  meta_data_subset = meta_data_subset,
                  plot_title = "CFRD Vs IGT Vs NGT no norm",
                  legend_title = "condition",
                  norm = "",
                  plot_file_name = "plots/pca_after_norm/pca_CFRDVsIGTVsNGT_nonorm.png")




#AU

meta_data_subset <- meta_data %>%
  filter(is.na(pre_post_modulator) | pre_post_modulator != 1) %>%
  filter(condition %in% c("CFRD", "IGT", "NGT")) %>%
  filter(country == "AU") %>%
  select(sample_long_name, condition, age_group, sex, country) %>%
  arrange(condition)


plot_pca_2factors(umi_counts = umi_counts,
                  meta_data_subset = meta_data_subset,
                  plot_title = "CFRD Vs IGT Vs NGT AU znorm_logCPM",
                  legend_title = "condition",
                  norm = "znorm_logCPM",
                  plot_file_name = "plots/pca_after_norm/pca_CFRDVsIGTVsNGT_AU_znorm_logCPM.png")
plot_pca_2factors(umi_counts = umi_counts,
                  meta_data_subset = meta_data_subset,
                  plot_title = "CFRD Vs IGT Vs NGT AU znorm_logTMM",
                  legend_title = "condition",
                  norm = "znorm_logTMM",
                  plot_file_name = "plots/pca_after_norm/pca_CFRDVsIGTVsNGT_AU_znorm_logTMM.png")
plot_pca_2factors(umi_counts = umi_counts,
                  meta_data_subset = meta_data_subset,
                  plot_title = "CFRD Vs IGT Vs NGT AU no norm",
                  legend_title = "condition",
                  norm = "",
                  plot_file_name = "plots/pca_after_norm/pca_CFRDVsIGTVsNGT_AU_nonorm.png")


#DK

meta_data_subset <- meta_data %>%
  filter(is.na(pre_post_modulator) | pre_post_modulator != 1) %>%
  filter(condition %in% c("CFRD", "IGT", "NGT")) %>%
  filter(country == "DK") %>%
  select(sample_long_name, condition, age_group, sex, country) %>%
  arrange(condition)


plot_pca_2factors(umi_counts = umi_counts,
                  meta_data_subset = meta_data_subset,
                  plot_title = "CFRD Vs IGT Vs NGT DK znorm_logCPM",
                  legend_title = "condition",
                  norm = "znorm_logCPM",
                  plot_file_name = "plots/pca_after_norm/pca_CFRDVsIGTVsNGT_DK_znorm_logCPM.png")
plot_pca_2factors(umi_counts = umi_counts,
                  meta_data_subset = meta_data_subset,
                  plot_title = "CFRD Vs IGT Vs NGT DK znorm_logTMM",
                  legend_title = "condition",
                  norm = "znorm_logTMM",
                  plot_file_name = "plots/pca_after_norm/pca_CFRDVsIGTVsNGT_DK_znorm_logTMM.png")
plot_pca_2factors(umi_counts = umi_counts,
                  meta_data_subset = meta_data_subset,
                  plot_title = "CFRD Vs IGT Vs NGT DK no norm",
                  legend_title = "condition",
                  norm = "",
                  plot_file_name = "plots/pca_after_norm/pca_CFRDVsIGTVsNGT_DK_nonorm.png")




####################################

# 2024 Feb 16
# check if all post modulator samples have pre modulator samples

phenotype.prot <- read.table("data/formatted/prot_phenotype_333_2024Jan.txt", header=TRUE, sep="\t")
phenotype.tra <- read.table("data/formatted/tra_phenotype_2024Jan.txt", header=TRUE, sep="\t")

postmod.prot <- phenotype.prot %>%
  dplyr::filter(PreModulatorVsPostModulator == "PostModulator") %>%
  dplyr::select(c(individual_id, sample_name, condition, batch_name))
premod.prot <- phenotype.prot %>%
  dplyr::filter(PreModulatorVsPostModulator == "PreModulator") %>%
  dplyr::select(c(individual_id, sample_name, condition, batch_name))

summary(factor(premod.prot$condition))
# CFRD                 HC                IGT                NGT UNKNOWN TO PREDICT 
# 38                 13                 38                 56                 22 

summary(factor(postmod.prot$condition))
# CFRD                IGT                NGT UNKNOWN TO PREDICT 
# 46                 33                 75                 12 

#how many post mod do not have pre mod
prot_matching_data <- postmod.prot %>%
  left_join(premod.prot, by = "individual_id", suffix = c("_postmod", "_premod"),
            relationship = "many-to-many")
sum(is.na(prot_matching_data$sample_name_premod))
#73
length(unique(prot_matching_data$sample_name_postmod))
#162
length(prot_matching_data$sample_name_postmod)
#188

#there are some replicate samples - so total number of samples should be obtained as
length(postmod.prot$sample_name)
#166

#73 out of 166 postmod do not have corresponding premod i.e. for same individual

prot_matching_data.nopremod <- prot_matching_data %>%
  dplyr::filter(is.na(sample_name_premod)) %>%
  dplyr::select(c(1:4))
write.csv(prot_matching_data.nopremod, "data/formatted/post_mod_with_no_premod.prot.csv", row.names = FALSE)



postmod.tra <- phenotype.tra %>%
  dplyr::filter(PreModulatorVsPostModulator == "PostModulator") %>%
  dplyr::select(c(individual_id, sample_name, condition, batch_name))
premod.tra <- phenotype.tra %>%
  dplyr::filter(PreModulatorVsPostModulator == "PreModulator") %>%
  dplyr::select(c(individual_id, sample_name, condition, batch_name))

summary(factor(premod.tra$condition))
# CFRD                 HC                IGT                NGT UNKNOWN TO PREDICT 
# 34                  6                 44                 62                 15 

summary(factor(postmod.tra$condition))
# CFRD                IGT                NGT UNKNOWN TO PREDICT 
# 48                 36                 76                 13 


#how many post mod do not have pre mod
tra_matching_data <- postmod.tra %>%
  left_join(premod.tra, by = "individual_id", suffix = c("_postmod", "_premod"),
            relationship = "many-to-many")
sum(is.na(tra_matching_data$sample_name_premod))
#66
length(unique(tra_matching_data$sample_name_postmod))
#173
length(tra_matching_data$sample_name_postmod)
#197

#66 out of 173 postmod do not have corresponding premod i.e. for same individual

tra_matching_data.nopremod <- tra_matching_data %>%
  dplyr::filter(is.na(sample_name_premod)) %>%
  dplyr::select(c(1:4))
write.csv(tra_matching_data.nopremod, "data/formatted/post_mod_with_no_premod.tra.csv", row.names = FALSE)


####################################

#check how many samples - pre and post have FEV values

phenotype.prot <- read.table("data/formatted/prot_phenotype_333_2024Jan.txt", header=TRUE, sep="\t") %>%
  mutate(FEV1 = gsub(" ", "", FEV1)) %>%
  mutate(FEV1 = ifelse(FEV1 == "NA", NA, FEV1)) %>%
  mutate(FEV1 = as.numeric(FEV1))
phenotype.tra <- read.table("data/formatted/tra_phenotype_2024Jan.txt", header=TRUE, sep="\t") %>%
  mutate(FEV1 = gsub(" ", "", FEV1)) %>%
  mutate(FEV1 = ifelse(FEV1 == "NA", NA, FEV1)) %>%
  mutate(FEV1 = as.numeric(FEV1))

sum(!is.na(phenotype.prot$FEV1))
#184
length(phenotype.prot$FEV1)
#333

sum(!is.na(phenotype.tra$FEV1))
#186
length(phenotype.tra$FEV1)
#334

phenotype.prot <- phenotype.prot %>%
  filter(condition %in% c("CFRD", "IGT", "NGT"))
phenotype.tra <- phenotype.tra %>%
  filter(condition %in% c("CFRD", "IGT", "NGT"))

sum(!is.na(phenotype.prot$FEV1))
#184
length(phenotype.prot$FEV1)
#286

sum(!is.na(phenotype.tra$FEV1))
#186
length(phenotype.tra$FEV1)
#300


phenotype.prot.pre <- phenotype.prot %>%
  dplyr::filter(PreModulatorVsPostModulator == "PreModulator")
phenotype.prot.post <- phenotype.prot %>%
  dplyr::filter(PreModulatorVsPostModulator == "PostModulator")

phenotype.tra.pre <- phenotype.tra %>%
  dplyr::filter(PreModulatorVsPostModulator == "PreModulator")
phenotype.tra.post <- phenotype.tra %>%
  dplyr::filter(PreModulatorVsPostModulator == "PostModulator")


sum(!is.na(phenotype.prot.pre$FEV1))
#109
length(phenotype.prot.pre$FEV1)
#132

sum(!is.na(phenotype.prot.post$FEV1))
#75
length(phenotype.prot.post$FEV1)
#154



sum(!is.na(phenotype.tra.pre$FEV1))
#108
length(phenotype.tra.pre$FEV1)
#140

sum(!is.na(phenotype.tra.post$FEV1))
#78
length(phenotype.tra.post$FEV1)
#160
