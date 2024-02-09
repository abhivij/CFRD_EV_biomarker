library(tidyverse)
library(readxl)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/CysticFibrosisGroup/ExoCF/CFRD_EV_biomarker/"
setwd(base_dir)

phenotype_tra <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t") %>%
  dplyr::select(c(Sample)) %>%
  dplyr::rename(c("tra_sample_name" = "Sample"))
phenotype_prot <- read.table("data/formatted/prot_phenotype_333.txt", header=TRUE, sep="\t") %>%
  dplyr::select(c(Sample, sample_long_name_t)) %>%
  dplyr::rename(c("prot_sample_name" = "Sample",
                  "tra_sample_name" = "sample_long_name_t"))
sum(!is.na(phenotype_prot$prot_sample_name))  
sum(!is.na(phenotype_prot$tra_sample_name))  

sum(!is.na(unique(phenotype_prot$prot_sample_name)))  
sum(!is.na(unique(phenotype_prot$tra_sample_name))) 

phenotype_combined <- phenotype_prot %>%
  full_join(phenotype_tra) %>%
  arrange(tra_sample_name)

sum(!is.na(phenotype_combined$prot_sample_name))  
#333
sum(!is.na(phenotype_combined$tra_sample_name)) 
#284

length(unique(phenotype_combined$tra_sample_name))
#this includes NA
#275

#so 274

length(unique(phenotype_combined$prot_sample_name))
#this includes NA
#334

#so 333

sum(is.na(phenotype_combined$prot_sample_name))
#35

sum(is.na(phenotype_combined$tra_sample_name))
#84

sum(is.na(phenotype_combined$prot_sample_name) & is.na(phenotype_combined$tra_sample_name))
#0

#add a new column for sample name
#can later change this sample name as per display requirement

#currently the sample name is transcriptomic sample name if sample is present otherwise proteomics
# (since - proteomics sample name is longer)
#but if proteomic replicates is present, then proteomics sample name is used 

replicate_prot_samples <- phenotype_combined %>%
  inner_join(phenotype_combined %>%
               filter(!is.na(tra_sample_name)) %>%
               group_by(tra_sample_name) %>%
               summarize(n = n()) %>%
               filter(n > 1))

#also later recreate this file after including new set of transcriptomic samples
phenotype_combined <- phenotype_combined %>%
  mutate(sample_name = case_when(prot_sample_name %in% replicate_prot_samples$prot_sample_name ~ prot_sample_name,
                                 is.na(tra_sample_name) ~ prot_sample_name,
                                 TRUE ~ tra_sample_name)) %>%
  arrange(tra_sample_name) %>%
  mutate(sample_name = sub("lfq_haronW.", "", sample_name, fixed = TRUE)) %>%
  mutate(sample_name = sub("lfq_W.", "", sample_name, fixed = TRUE)) %>%
  mutate(sample_name = sub("^X", "", sample_name))

write.table(replicate_prot_samples, 
            file = "data/formatted/replicate_prot_samples.txt", sep="\t", row.names=FALSE)

write.table(phenotype_combined, 
            file = "data/formatted/pt_mapping_and_common_sample_name.txt", sep="\t", row.names=FALSE)

phenotype_combined2 <- read.table("data/formatted/pt_mapping_and_common_sample_name.txt", header = TRUE)

all.equal(phenotype_combined, phenotype_combined2)




###########################
#reecreating mapping file using 2024 Jan phenotype files



phenotype_tra <- read.table("data/formatted/tra_phenotype_2024Jan.txt", header=TRUE, sep="\t") %>%
  dplyr::select(c(Sample)) %>%
  dplyr::rename(c("tra_sample_name" = "Sample"))
phenotype_prot <- read.table("data/formatted/prot_phenotype_333_2024Jan.txt", header=TRUE, sep="\t") %>%
  dplyr::select(c(Sample, sample_long_name_t)) %>%
  dplyr::rename(c("prot_sample_name" = "Sample",
                  "tra_sample_name" = "sample_long_name_t"))
sum(!is.na(phenotype_prot$prot_sample_name)) 
# 333
sum(!is.na(phenotype_prot$tra_sample_name))  
# 311

sum(!is.na(unique(phenotype_prot$prot_sample_name)))
# 333

sum(!is.na(unique(phenotype_prot$tra_sample_name))) 
# 301

phenotype_combined <- phenotype_prot %>%
  full_join(phenotype_tra) %>%
  arrange(tra_sample_name)

sum(!is.na(phenotype_combined$prot_sample_name))  
#333
sum(!is.na(phenotype_combined$tra_sample_name)) 
#344

length(unique(phenotype_combined$tra_sample_name))
#this includes NA
#335

#so 334

length(unique(phenotype_combined$prot_sample_name))
#this includes NA
#334

#so 333

sum(is.na(phenotype_combined$prot_sample_name))
#33

sum(is.na(phenotype_combined$tra_sample_name))
#22

sum(is.na(phenotype_combined$prot_sample_name) & is.na(phenotype_combined$tra_sample_name))
#0

#add a new column for sample name
#can later change this sample name as per display requirement

#currently the sample name is transcriptomic sample name if sample is present otherwise proteomics
# (since - proteomics sample name is longer)
#but if proteomic replicates is present, then proteomics sample name is used 

replicate_prot_samples <- phenotype_combined %>%
  inner_join(phenotype_combined %>%
               filter(!is.na(tra_sample_name)) %>%
               group_by(tra_sample_name) %>%
               summarize(n = n()) %>%
               filter(n > 1))

phenotype_combined <- phenotype_combined %>%
  mutate(sample_name = case_when(prot_sample_name %in% replicate_prot_samples$prot_sample_name ~ prot_sample_name,
                                 is.na(tra_sample_name) ~ prot_sample_name,
                                 TRUE ~ tra_sample_name)) %>%
  arrange(tra_sample_name) %>%
  mutate(sample_name = sub("lfq_haronW.", "", sample_name, fixed = TRUE)) %>%
  mutate(sample_name = sub("lfq_W.", "", sample_name, fixed = TRUE)) %>%
  mutate(sample_name = sub("^X", "", sample_name))

write.table(replicate_prot_samples, 
            file = "data/formatted/replicate_prot_samples_Jan2024.txt", sep="\t", row.names=FALSE)

write.table(phenotype_combined, 
            file = "data/formatted/pt_mapping_and_common_sample_name_Jan2024.txt", sep="\t", row.names=FALSE)

phenotype_combined2 <- read.table("data/formatted/pt_mapping_and_common_sample_name_Jan2024.txt", header = TRUE)

all.equal(phenotype_combined, phenotype_combined2)
