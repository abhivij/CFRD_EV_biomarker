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

