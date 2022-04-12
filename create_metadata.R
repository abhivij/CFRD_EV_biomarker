library(tidyverse)
library(readxl)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/CysticFibrosisGroup/ExoCF/CFRD_EV_biomarker/"
setwd(base_dir)

sample_info <- read.table("data/WAT10166-WAT10180-WAT9739_md5sum.txt", sep = " ") %>%
  select(V3) %>%
  separate(V3, into = c(NA, NA, NA, "file_name"), sep = "/") %>%
  separate("file_name", into = c("sample_long_name", NA), sep = "_L") %>%
  unique() %>%
  separate("sample_long_name", into = c("sample_name", "illumina_sample_number"), sep = "_S", remove = FALSE) %>%
  mutate(illumina_sample_number = strtoi(illumina_sample_number)) %>%
  arrange(sample_name)

quantification_batch <- read.csv("data/formatted/quantification_batch.csv")
sample_info <- sample_info %>%
  inner_join(quantification_batch)

sample_info <- sample_info %>%
  mutate(individual_id = case_when(
    grepl("2020_t0_CPH", sample_name, fixed = TRUE) ~ gsub(pattern = "2020_t0_", 
                                                           replacement = "",
                                                           x = sample_name, 
                                                           fixed = TRUE),
    grepl("2017_CPH", sample_name, fixed = TRUE) ~ gsub(pattern = "2017_", 
                                                           replacement = "",
                                                           x = sample_name, 
                                                           fixed = TRUE),
    grepl("CPH", sample_name, fixed = TRUE) ~ sample_name,
    grepl("17_", sample_name, fixed = TRUE) ~ substr(sample_name, 7, 11),  # 17_20_#####
                                                                             ##### is the indiv_id 
    grepl("11_", sample_name, fixed = TRUE) ~ substr(sample_name, 7, 11),
    grepl("14E_", sample_name, fixed = TRUE) ~ substr(sample_name, 8, 12),
    TRUE ~ NA_character_))
sample_info <- sample_info %>%
  mutate(year = case_when(
    grepl("2020_t0_CPH", sample_name, fixed = TRUE) ~ "2020",
    grepl("2017_CPH", sample_name, fixed = TRUE) ~ "2017",
    grepl("CPH", sample_name, fixed = TRUE) ~ "unknown",
    grepl("17_", sample_name, fixed = TRUE) ~ paste0("20", substr(sample_name, 4, 5)),  # 17_20_#####
    grepl("11_", sample_name, fixed = TRUE) ~ paste0("20", substr(sample_name, 4, 5)),
    grepl("14E_", sample_name, fixed = TRUE) ~ paste0("20", substr(sample_name, 5, 6)),
    TRUE ~ NA_character_))

#from CopenhagenClinicalData.xlsx and from proteomics metadata,
#   CPH10, CPH22 are of year 2020
sample_info <- sample_info %>%
  mutate(year = ifelse((sample_name %in% c("CPH10", "CPH22")), "2020", year))

sample_info <- sample_info %>%
  select(sample_long_name, illumina_sample_number, quant_batch, 
         sample_name, individual_id, year) %>%
  mutate(biotype = "serum", .after = quant_batch)

write.csv(sample_info, "data/formatted/sample_info_temp.txt", row.names = FALSE)

proteomic_metadata <- read.table("data/proteomics_sinfo.txt", header = TRUE, sep = "\t") %>%
  mutate(individualid = gsub("CPH0", "CPH", individualid, fixed = TRUE)) %>%
  filter(!individualid %in% c("glufib", "QC")) %>%
  select(c(rawfile, label, individualid, condition, cohort, 
           agegroup, modulator, biotype, year, technicalreplicate, FEV1)) 

proteomic_metadata <- proteomic_metadata %>%
  filter(technicalreplicate == "t1")

# proteomic_metadata_multiple_entries <- proteomic_metadata %>%
#   group_by(individualid, year) %>%
#   summarise(n = n()) %>%
#   arrange(desc(n))  
# 
# write.csv(proteomic_metadata, "data/proteomic_processing.txt", row.names = FALSE)
# 
# write.csv(proteomic_metadata_multiple_entries, "data/proteomic_multiple_entries.txt", row.names = FALSE)

# proteomic_metadata_single_entry <- proteomic_metadata_multiple_entries %>%
#   filter(n == 1)
# #222
# 
# proteomic_copenhagen <- proteomic_metadata %>%
#   filter(cohort == "DK")
#ignore proteomic metadata - use master sheet



# length(unique(sample_info$sample_name))
# #272
# length(unique(sample_info$individual_id))
# #180
# 
# length(unique(proteomic_metadata$individualid))
# #174


# sample_info_multiple_entries <- sample_info %>%
#   group_by(individual_id, year) %>%
#   summarise(n = n()) %>%
#   arrange(desc(n))  
# sample_info_single_entry <- sample_info_multiple_entries %>%
#   filter(n == 1)
#194



####obtain other data from master sheet

#copenhagen
sample_info.cph <- sample_info %>%
  filter(grepl("CPH", sample_name, fixed = TRUE))
length(unique(sample_info.cph$individual_id))

mastersheet_info.cph <- read_xlsx("data/ExoCF mastersheet_RPA Copenhagen SCH.xlsx",
                              sheet = 3)
mastersheet_info.cph <- mastersheet_info.cph %>%
  select(-c(8:18))
mastersheet_info.cph <- mastersheet_info.cph[2:373,]
mastersheet_info.cph <- mastersheet_info.cph %>%
  filter(`EXOSOME ISOLATION` != "Not isolated")

colnames(mastersheet_info.cph)[1] <- "diabetes_status"
which(mastersheet_info.cph$diabetes_status == "4=CFRD")
which(mastersheet_info.cph$diabetes_status == "3=IGT")

mastersheet_info.cph[c(which(mastersheet_info.cph$diabetes_status == "4=CFRD"):
                         (which(mastersheet_info.cph$diabetes_status == "3=IGT") - 1)), 
                     "diabetes_status"] <- "CFRD"
mastersheet_info.cph[c(which(mastersheet_info.cph$diabetes_status == "3=IGT"):
                         (which(mastersheet_info.cph$diabetes_status == "2=Indeterminate") - 1)), 
                     "diabetes_status"] <- "IGT"
mastersheet_info.cph[c(which(mastersheet_info.cph$diabetes_status == "2=Indeterminate"):
                         (which(mastersheet_info.cph$diabetes_status == "1=NGT") - 1)), 
                     "diabetes_status"] <- "IND"
mastersheet_info.cph[c(which(mastersheet_info.cph$diabetes_status == "1=NGT"):
                         dim(mastersheet_info.cph)[1]), 
                     "diabetes_status"] <- "NGT"

mastersheet_info.cph[c(2:3), 2] <- 9
mastersheet_info.cph[c(135:137), 3] <- 59

length(unique(mastersheet_info.cph$`Study ID (2020)`))
#93
length(unique(mastersheet_info.cph$`CF number (2017)`))
#93

colnames(mastersheet_info.cph)[4] <- "intake"
mastersheet_info.cph <- mastersheet_info.cph %>%
  filter(!intake %in% c("2020 Serum t=120", "2017 Plasma"))

colnames(mastersheet_info.cph)[6] <- "is_sample_available"
mastersheet_info.cph <- mastersheet_info.cph %>%
  filter(is_sample_available == 1)

colnames(mastersheet_info.cph)[c(2,3)] <- c("study_id", "CF_number")
mastersheet_info.cph <- mastersheet_info.cph %>%
  mutate(individual_id = paste0("CPH", study_id), .after = study_id) %>%
  separate(intake, into = c("year", NA, NA), sep = " ", remove = FALSE)
mastersheet_info.cph <- mastersheet_info.cph %>%
  select(c(1:17), -c(is_sample_available, `EXOSOME ISOLATION`)) 

mastersheet_info.cph <- mastersheet_info.cph %>%
  mutate(sample_name = case_when(individual_id %in% c("CPH10", "CPH22") ~ individual_id,
                                 year == "2017" ~ paste("2017", individual_id, sep = "_"),
                                 year == "2020" ~ paste("2020_t0", individual_id, sep = "_"),
                                 TRUE ~ NA_character_), 
         .after = year) %>%
  select(-c(individual_id, year))

# sample_info_sub <- sample_info.cph %>%
#   select(sample_name, quant_batch)
# 
# cph_info <- mastersheet_info.cph %>%
#   left_join(sample_info_sub)
#2020_t0_CPH24 entry not present in 272 samples


cph_info <- sample_info.cph %>%
  inner_join(mastersheet_info.cph)

summary(factor(cph_info$`Sex 1=male, 2=female`))

cph_info <- cph_info[,-c(9:11, 16:19)]
colnames(cph_info)[9:13] <- c("pre_post_modulator", "cftr_modulator_2020", "sex", "age", "FEV1")

cph_info <- cph_info %>%
  mutate(sex = case_when(sex == 1 ~ "M",
                         sex == 2 ~ "F",
                         TRUE ~ NA_character_))
  
summary(factor(cph_info$sex))

write.csv(cph_info, "data/formatted/sample_info_cph.csv", row.names = FALSE)
