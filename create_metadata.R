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

write.csv(sample_info, "data/formatted/sample_info_temp.txt", row.names = FALSE)

proteomic_metadata <- read.table("data/proteomics_sinfo.txt", header = TRUE, sep = "\t") %>%
  mutate(individualid = gsub("CPH0", "CPH", individualid, fixed = TRUE)) %>%
  filter(!individualid %in% c("glufib", "QC")) %>%
  select(c(individualid, condition, cohort, agegroup, modulator, biotype, year, FEV1)) %>%
  unique() 

write.csv(proteomic_metadata, "data/proteomic_processing.txt", row.names = FALSE)

proteomic_metadata <- proteomic_metadata %>%
  group_by(individualid, year) %>%
  summarise(n = n()) %>%
  arrange(desc(n))
write.csv(proteomic_metadata, "data/proteomic_multiple_entries.txt", row.names = FALSE)

length(unique(proteomic_metadata$individualid))
