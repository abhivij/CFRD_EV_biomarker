library(tidyverse)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/CysticFibrosisGroup/ExoCF/CFRD_EV_biomarker/"
setwd(base_dir)

meta_data <- read.csv("data/formatted/meta_data.csv")

phenotype <- meta_data %>%
  rename("Sample" = "sample_long_name") %>%
  mutate(Sample = paste0("X", Sample))

phenotype <- phenotype %>%
  mutate("CFRDVsIGT" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ NA_character_,
                                 condition == "CFRD" ~ "CFRD",
                                 condition == "IGT" ~ "IGT",
                                 TRUE ~ NA_character_))
summary(factor(phenotype$CFRDVsIGT))
phenotype <- phenotype %>%
  mutate("CFRDVsNGT" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ NA_character_,
                                 condition == "CFRD" ~ "CFRD",
                                 condition == "NGT" ~ "NGT",
                                 TRUE ~ NA_character_))
summary(factor(phenotype$CFRDVsNGT))

phenotype <- phenotype %>%
  mutate("IGTVsNGT" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ NA_character_,
                                condition == "IGT" ~ "IGT",
                                condition == "NGT" ~ "NGT",
                                TRUE ~ NA_character_))
summary(factor(phenotype$IGTVsNGT))

phenotype <- phenotype %>%
  mutate("CFVsHC" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ NA_character_,
                                condition == "HC" ~ "HC",
                                TRUE ~ "CF"))
summary(factor(phenotype$CFVsHC))

phenotype <- phenotype %>%
  mutate("PreModulatorVsPostModulator" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ "PostModulator",
                              TRUE ~ "PreModulator"))

#two samples are not prepended with X in the data file. Replicating that here
phenotype <- phenotype %>%
  mutate(Sample = sub("XCPH10_S271", "CPH10_S271", Sample)) %>%
  mutate(Sample = sub("XCPH22_S272", "CPH22_S272", Sample))


write.table(phenotype, 
            file = "data/formatted/phenotype.txt", sep="\t", row.names=FALSE)

p2 <- read.table("data/formatted/phenotype.txt", header = TRUE)

all.equal(phenotype, p2)



#test pipeline read

data <- read.table("data/formatted/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                   nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")

all.equal(umi_counts, data)

sum(is.na(data))
#data format : (transcripts x samples)
data[is.na(data)] <- 0
phenotype <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t")

filter = expression(country == 'AU')
classification_criteria <- "CFRDVsIGT"

extracted_samples <- phenotype %>% subset(eval(filter))
extracted_samples <- extracted_samples[!is.na(extracted_samples[classification_criteria]), ]
extracted_samples$Sample <- factor(extracted_samples$Sample)

filtered_samples_read_count <- data %>% dplyr::select(extracted_samples$Sample)

#from the extracted_samples, select the 'Sample' column and classification_criteria column
filtered_samples_output_labels <- extracted_samples[, c('Sample', classification_criteria)]
colnames(filtered_samples_output_labels) <- c("Sample", "Label")




#below code for creating new phenotype outdated after adding this information in meta data
# ###################### creating new phenotype file with updating CFRD/IGT/NGT status for some child samples
# 
# results <- read_excel("prediction_pipeline/sch_pred_with_clinical_results.xlsx")
# results <- results[-c(1), c(2, 6, 8, 9, 10, 11)]
# colnames(results) <- c("Sample", "pre_post_mod", "prediction", "matching_dates", "OGTT", "actual")
# new_label_info <- results %>%
#   filter(matching_dates == "P") %>%
#   filter(!is.na(actual)) %>%
#   dplyr::select(c(Sample, actual, pre_post_mod)) %>%
#   mutate(Sample = paste0("X", Sample))
# 
# phenotype <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t") %>%
#   left_join(new_label_info)
# phenotype_new <- phenotype %>%
#   mutate(condition = case_when(!is.na(actual) ~ actual,
#                                TRUE ~ condition))
# 
# sum(phenotype_new$condition != phenotype$condition)
# sum(new_label_info$Sample %in% phenotype$Sample)
# new_label_info$Sample[!(new_label_info$Sample %in% phenotype$Sample)]
# 
# 
# phenotype_new <- phenotype_new %>%
#   dplyr::select(-c(actual, pre_post_mod))
# 
# phenotype_new <- phenotype_new %>%
#   mutate("CFRDVsIGT" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ NA_character_,
#                                  condition == "CFRD" ~ "CFRD",
#                                  condition == "IGT" ~ "IGT",
#                                  TRUE ~ NA_character_))
# summary(factor(phenotype_new$CFRDVsIGT))
# phenotype_new <- phenotype_new %>%
#   mutate("CFRDVsNGT" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ NA_character_,
#                                  condition == "CFRD" ~ "CFRD",
#                                  condition == "NGT" ~ "NGT",
#                                  TRUE ~ NA_character_))
# summary(factor(phenotype_new$CFRDVsNGT))
# 
# phenotype_new <- phenotype_new %>%
#   mutate("IGTVsNGT" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ NA_character_,
#                                 condition == "IGT" ~ "IGT",
#                                 condition == "NGT" ~ "NGT",
#                                 TRUE ~ NA_character_))
# summary(factor(phenotype_new$IGTVsNGT))
# 
# 
# 
# summary(factor(phenotype$CFRDVsIGT))
# summary(factor(phenotype$CFRDVsNGT))
# summary(factor(phenotype$IGTVsNGT))
# 
# 
# new_label_info_sub <- new_label_info %>%
#   filter(pre_post_mod == 0)
# 
# 
# write.table(phenotype_new, 
#             file = "data/formatted/phenotype_new.txt", sep="\t", row.names=FALSE)




#check counts of different groups
meta_data <- read.csv("data/formatted/meta_data.csv") %>%
  mutate(country_condition_ag = paste(country, age_group, condition))

non_modulator <- meta_data %>% filter(is.na(pre_post_modulator) | pre_post_modulator == 0)
modulator <- meta_data %>% filter(!is.na(pre_post_modulator) & pre_post_modulator == 1)

summary(factor(non_modulator$country))
summary(factor(non_modulator$country_condition_ag))

summary(factor(modulator$country))
summary(factor(modulator$country_condition_ag))


phenotype <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t") %>%
  mutate(country_condition_ag = paste(country, age_group, condition))

non_modulator <- phenotype %>% filter(is.na(pre_post_modulator) | pre_post_modulator == 0)
modulator <- phenotype %>% filter(!is.na(pre_post_modulator) & pre_post_modulator == 1)

summary(factor(non_modulator$country))
summary(factor(non_modulator$country_condition_ag))

summary(factor(modulator$country))
summary(factor(modulator$country_condition_ag))

#######################################################
#######################################################
#######################################################
#updated phenotype
library(tidyverse)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/CysticFibrosisGroup/ExoCF/CFRD_EV_biomarker/"
setwd(base_dir)

meta_data <- read.csv("data/formatted/meta_data_updated.csv") %>%
  rename(c("condition" = "condition_updated"))

summary(factor(meta_data$condition))

phenotype <- meta_data %>%
  rename("Sample" = "sample_long_name") %>%
  mutate(Sample = paste0("X", Sample))

phenotype <- phenotype %>%
  mutate("CFRDVsIGT" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ NA_character_,
                                 condition == "CFRD" ~ "CFRD",
                                 condition == "IGT" ~ "IGT",
                                 TRUE ~ NA_character_))
summary(factor(phenotype$CFRDVsIGT))
phenotype <- phenotype %>%
  mutate("CFRDVsNGT" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ NA_character_,
                                 condition == "CFRD" ~ "CFRD",
                                 condition == "NGT" ~ "NGT",
                                 TRUE ~ NA_character_))
summary(factor(phenotype$CFRDVsNGT))

phenotype <- phenotype %>%
  mutate("IGTVsNGT" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ NA_character_,
                                condition == "IGT" ~ "IGT",
                                condition == "NGT" ~ "NGT",
                                TRUE ~ NA_character_))
summary(factor(phenotype$IGTVsNGT))

phenotype <- phenotype %>%
  mutate("CFVsHC" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ NA_character_,
                              condition == "HC" ~ "HC",
                              TRUE ~ "CF"))
summary(factor(phenotype$CFVsHC))

phenotype <- phenotype %>%
  mutate("PreModulatorVsPostModulator" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ "PostModulator",
                                                   TRUE ~ "PreModulator"))

#two samples are not prepended with X in the data file. Replicating that here
phenotype <- phenotype %>%
  mutate(Sample = sub("XCPH10_S271", "CPH10_S271", Sample)) %>%
  mutate(Sample = sub("XCPH22_S272", "CPH22_S272", Sample))


write.table(phenotype, 
            file = "data/formatted/phenotype.txt", sep="\t", row.names=FALSE)

p2 <- read.table("data/formatted/phenotype.txt", header = TRUE)

all.equal(phenotype, p2)


#######################################################
#phenotype for proteomics

meta_data <- read.csv("data/proteomics/prot_metadata_updated.csv")
phenotype <- meta_data %>%
  dplyr::rename(c("Sample" = "label")) %>%
  mutate(Sample = paste0("lfq_", Sample))

phenotype <- phenotype %>%
  mutate("CFRDVsIGT" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ NA_character_,
                                 condition == "CFRD" ~ "CFRD",
                                 condition == "IGT" ~ "IGT",
                                 TRUE ~ NA_character_))
summary(factor(phenotype$CFRDVsIGT))
phenotype <- phenotype %>%
  mutate("CFRDVsNGT" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ NA_character_,
                                 condition == "CFRD" ~ "CFRD",
                                 condition == "NGT" ~ "NGT",
                                 TRUE ~ NA_character_))
summary(factor(phenotype$CFRDVsNGT))

phenotype <- phenotype %>%
  mutate("IGTVsNGT" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ NA_character_,
                                condition == "IGT" ~ "IGT",
                                condition == "NGT" ~ "NGT",
                                TRUE ~ NA_character_))
summary(factor(phenotype$IGTVsNGT))

phenotype <- phenotype %>%
  mutate("CFVsHC" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ NA_character_,
                              condition == "HC" ~ "HC",
                              TRUE ~ "CF"))
summary(factor(phenotype$CFVsHC))

phenotype <- phenotype %>%
  mutate("PreModulatorVsPostModulator" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ "PostModulator",
                                                   TRUE ~ "PreModulator"))

phenotype_only_main <- phenotype %>%
  filter(mq_batch == "main")


#check data read
data.imputed.combined <- read.csv("data/proteomics/imputed_combined.csv", row.names = 1) %>%
  dplyr::select(phenotype$Sample)

data.imputed.main <- read.csv("data/proteomics/imputed_main.csv", row.names = 1) %>%
  dplyr::select(phenotype_only_main$Sample)


write.table(phenotype, 
            file = "data/formatted/prot_phenotype.txt", sep="\t", row.names=FALSE)
write.table(phenotype_only_main, 
            file = "data/formatted/prot_phenotype_only_main.txt", sep="\t", row.names=FALSE)


#######################################################
#phenotype for proteomics with 333 samples - includes the 2023 DK samples

meta_data <- read.csv("data/proteomics/prot_metadata_all_2023Aug.csv")
phenotype <- meta_data %>%
  mutate(label = gsub(pattern = "-", replacement = ".", label, fixed = TRUE)) %>%
  dplyr::rename(c("Sample" = "label")) %>%
  mutate(Sample = paste0("lfq_", Sample))

phenotype <- phenotype %>%
  mutate("CFRDVsIGT" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ NA_character_,
                                 condition == "CFRD" ~ "CFRD",
                                 condition == "IGT" ~ "IGT",
                                 TRUE ~ NA_character_))
summary(factor(phenotype$CFRDVsIGT))
phenotype <- phenotype %>%
  mutate("CFRDVsNGT" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ NA_character_,
                                 condition == "CFRD" ~ "CFRD",
                                 condition == "NGT" ~ "NGT",
                                 TRUE ~ NA_character_))
summary(factor(phenotype$CFRDVsNGT))

phenotype <- phenotype %>%
  mutate("IGTVsNGT" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ NA_character_,
                                condition == "IGT" ~ "IGT",
                                condition == "NGT" ~ "NGT",
                                TRUE ~ NA_character_))
summary(factor(phenotype$IGTVsNGT))

phenotype <- phenotype %>%
  mutate("CFVsHC" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ NA_character_,
                              condition == "HC" ~ "HC",
                              TRUE ~ "CF"))
summary(factor(phenotype$CFVsHC))

phenotype <- phenotype %>%
  mutate("PreModulatorVsPostModulator" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ "PostModulator",
                                                   TRUE ~ "PreModulator"))
summary(factor(phenotype$PreModulatorVsPostModulator))

phenotype_no_other <- phenotype %>%
  filter(batch_name != "other")


#check data read
data.imputed.mf <- read.csv("data/proteomics/data_333samples_imputed_mf.csv", row.names = 1) %>%
  dplyr::select(phenotype$Sample)
data.imputed.zero <- read.csv("data/proteomics/data_333samples_imputed_zero.csv", row.names = 1) %>%
  dplyr::select(phenotype$Sample)

data_no_other.imputed.mf <- read.csv("data/proteomics/data_no_other_315samples_imputed_mf.csv", 
                                     row.names = 1) %>%
  dplyr::select(phenotype_no_other$Sample)
data_no_other.imputed.zero <- read.csv("data/proteomics/data_no_other_315samples_imputed_zero.csv", 
                                       row.names = 1) %>%
  dplyr::select(phenotype_no_other$Sample)


write.table(phenotype, 
            file = "data/formatted/prot_phenotype_333.txt", sep="\t", row.names=FALSE)
write.table(phenotype_no_other, 
            file = "data/formatted/prot_phenotype_no_other_315.txt", sep="\t", row.names=FALSE)




#########################################
#phenotype for 334 samples RNA samples

meta_data <- read.csv("data/formatted/tra_metadata_all_2023Oct.csv")
phenotype <- meta_data %>%
  rename("Sample" = "sample_long_name") %>%
  mutate(Sample = paste0("X", Sample))

phenotype <- phenotype %>%
  mutate("CFRDVsIGT" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ NA_character_,
                                 condition == "CFRD" ~ "CFRD",
                                 condition == "IGT" ~ "IGT",
                                 TRUE ~ NA_character_))
summary(factor(phenotype$CFRDVsIGT))
phenotype <- phenotype %>%
  mutate("CFRDVsNGT" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ NA_character_,
                                 condition == "CFRD" ~ "CFRD",
                                 condition == "NGT" ~ "NGT",
                                 TRUE ~ NA_character_))
summary(factor(phenotype$CFRDVsNGT))

phenotype <- phenotype %>%
  mutate("IGTVsNGT" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ NA_character_,
                                condition == "IGT" ~ "IGT",
                                condition == "NGT" ~ "NGT",
                                TRUE ~ NA_character_))
summary(factor(phenotype$IGTVsNGT))

phenotype <- phenotype %>%
  mutate("CFVsHC" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ NA_character_,
                              condition == "HC" ~ "HC",
                              TRUE ~ "CF"))
summary(factor(phenotype$CFVsHC))

phenotype <- phenotype %>%
  mutate("PreModulatorVsPostModulator" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ "PostModulator",
                                                   TRUE ~ "PreModulator"))
summary(factor(phenotype$PreModulatorVsPostModulator))
#two samples are not prepended with X in the data file. Replicating that here
phenotype <- phenotype %>%
  mutate(Sample = sub("XCPH10_S271", "CPH10_S271", Sample)) %>%
  mutate(Sample = sub("XCPH22_S272", "CPH22_S272", Sample))

phenotype <- phenotype %>%
  arrange(condition, batch_name, Sample)

write.table(phenotype, 
            file = "data/formatted/tra_phenotype_2024Jan.txt", sep="\t", row.names=FALSE)

#test reading data
data <- read.csv("data/formatted/rna_all/umi_counts_filter90.csv", row.names = 1)

data <- data %>%
  dplyr::select(c(phenotype$Sample))
#test reading data end


#######################################################
#phenotype for proteomics with 333 samples - includes the 2023 DK samples with tra sample names

meta_data <- read.csv("data/proteomics/prot_metadata_all_2023Oct.csv")
phenotype <- meta_data %>%
  mutate(label = gsub(pattern = "-", replacement = ".", label, fixed = TRUE)) %>%
  dplyr::rename(c("Sample" = "label")) %>%
  mutate(Sample = paste0("lfq_", Sample))

phenotype <- phenotype %>%
  mutate("CFRDVsIGT" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ NA_character_,
                                 condition == "CFRD" ~ "CFRD",
                                 condition == "IGT" ~ "IGT",
                                 TRUE ~ NA_character_))
summary(factor(phenotype$CFRDVsIGT))
phenotype <- phenotype %>%
  mutate("CFRDVsNGT" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ NA_character_,
                                 condition == "CFRD" ~ "CFRD",
                                 condition == "NGT" ~ "NGT",
                                 TRUE ~ NA_character_))
summary(factor(phenotype$CFRDVsNGT))

phenotype <- phenotype %>%
  mutate("IGTVsNGT" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ NA_character_,
                                condition == "IGT" ~ "IGT",
                                condition == "NGT" ~ "NGT",
                                TRUE ~ NA_character_))
summary(factor(phenotype$IGTVsNGT))

phenotype <- phenotype %>%
  mutate("CFVsHC" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ NA_character_,
                              condition == "HC" ~ "HC",
                              TRUE ~ "CF"))
summary(factor(phenotype$CFVsHC))

phenotype <- phenotype %>%
  mutate("PreModulatorVsPostModulator" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ "PostModulator",
                                                   TRUE ~ "PreModulator"))
summary(factor(phenotype$PreModulatorVsPostModulator))

phenotype <- phenotype %>%
  mutate(sample_long_name_t = ifelse(is.na(sample_long_name_t) | grepl("^X.*", sample_long_name_t, fixed = FALSE),
                                     sample_long_name_t,
                                     paste0("X", sample_long_name_t)))
phenotype <- phenotype %>%
  mutate(sample_long_name_t = sub("XCPH10_S271", "CPH10_S271", sample_long_name_t)) %>%
  mutate(sample_long_name_t = sub("XCPH22_S272", "CPH22_S272", sample_long_name_t))

phenotype_no_other <- phenotype %>%
  filter(batch_name != "other")


#check data read
data.imputed.mf <- read.csv("data/proteomics/data_333samples_imputed_mf.csv", row.names = 1) %>%
  dplyr::select(phenotype$Sample)
data.imputed.zero <- read.csv("data/proteomics/data_333samples_imputed_zero.csv", row.names = 1) %>%
  dplyr::select(phenotype$Sample)

data_no_other.imputed.mf <- read.csv("data/proteomics/data_no_other_315samples_imputed_mf.csv", 
                                     row.names = 1) %>%
  dplyr::select(phenotype_no_other$Sample)
data_no_other.imputed.zero <- read.csv("data/proteomics/data_no_other_315samples_imputed_zero.csv", 
                                       row.names = 1) %>%
  dplyr::select(phenotype_no_other$Sample)


write.table(phenotype, 
            file = "data/formatted/prot_phenotype_333_2024Jan.txt", sep="\t", row.names=FALSE)
write.table(phenotype_no_other, 
            file = "data/formatted/prot_phenotype_no_other_315_2024Jan.txt", sep="\t", row.names=FALSE)

phenotype <- read.table("data/formatted/prot_phenotype_333_2024Jan.txt", header = TRUE, sep = "\t")
phenotype.old <- read.table("data/formatted/prot_phenotype_333.txt", header = TRUE, sep = "\t")

all.equal(phenotype.old, phenotype)
#not equal
all.equal(phenotype.old %>% dplyr::select(-c(sample_long_name_t)), 
          phenotype %>% dplyr::select(-c(sample_long_name_t)))
#TRUE


phenotype.prot <- read.table("data/formatted/prot_phenotype_333_2024Jan.txt", header = TRUE, sep = "\t")

phenotype.tra <- read.table("data/formatted/tra_phenotype_2024Jan.txt", header = TRUE, sep = "\t")

tra_sample_names <- unique(phenotype.prot$sample_long_name_t[!is.na(phenotype.prot$sample_long_name_t)])
tra_sample_names_all <- unique(phenotype.tra$Sample)

missing_in_prot <- setdiff(tra_sample_names_all, tra_sample_names)
#33

missing_in_tra <- setdiff(tra_sample_names, tra_sample_names_all) 
#0

prot_sample_names_all <- unique(phenotype.prot$rawfile)
#334 - but includes replicates



#check batches in cfrd, igt
phenotype <- read.table("data/formatted/prot_phenotype_333_2024Jan.txt", header = TRUE, sep = "\t")
summary(factor(phenotype$batch_name))
# main   new other 
# 252    63    18   

phenotype_sub <- phenotype %>%
  filter(!is.na(CFRDVsIGT))
summary(factor(phenotype_sub$batch_name))
# main other 
# 67     9 

phenotype_sub_2 <- phenotype %>%
  filter(batch_name == "new")
summary(factor(phenotype_sub_2$condition))
# CFRD  IGT  NGT 
# 16   16   31 

summary(factor(phenotype_sub_2$PreModulatorVsPostModulator))
# PostModulator 
# 63 

#verified that all samples from new batch are post modulator