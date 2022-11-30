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

write.table(phenotype, 
            file = "data/formatted/phenotype.txt", quote=FALSE, sep="\t", row.names=FALSE)

p2 <- read.table("data/formatted/phenotype.txt", header = TRUE)

all.equal(phenotype, p2)
#types of age, FEV1 differ - that's okay

all.equal(p2, p2_new %>% select(-c(CFVsHC, PreModulatorVsPostModulator)))




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
