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

#special cases - samples named differently in master sheet
sample_info <- sample_info %>%
  mutate(sample_name = sub("034MJC", "034MC", sample_name, fixed = TRUE)) %>%
  mutate(sample_name = sub("050JMC", "050JM", sample_name, fixed = TRUE))

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



# write.csv(sample_info, "data/formatted/sample_info_temp.txt", row.names = FALSE)
# 
# proteomic_metadata <- read.table("data/proteomics_sinfo.txt", header = TRUE, sep = "\t") %>%
#   mutate(individualid = gsub("CPH0", "CPH", individualid, fixed = TRUE)) %>%
#   filter(!individualid %in% c("glufib", "QC")) %>%
#   select(c(rawfile, label, individualid, condition, cohort, 
#            agegroup, modulator, biotype, year, technicalreplicate, FEV1)) 
# 
# proteomic_metadata <- proteomic_metadata %>%
#   filter(technicalreplicate == "t1")

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
  select(c(1:17, 19, 20), -c(is_sample_available, `EXOSOME ISOLATION`)) 
colnames(mastersheet_info.cph)[16:17] <- c("OGTT_2020", "OGTT_2017")

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

cph_info <- cph_info[,-c(9:11, 18:19)]
colnames(cph_info)[9:15] <- c("pre_post_modulator", "cftr_modulator_2020", "sex", "age", 
                              "mutation1", "mutation2",
                              "FEV1")

cph_info <- cph_info %>%
  mutate(sex = case_when(sex == 1 ~ "M",
                         sex == 2 ~ "F",
                         TRUE ~ NA_character_))

cph_info <- cph_info %>%
  relocate(FEV1, .before = mutation1)
  
summary(factor(cph_info$sex))

length(unique(cph_info$individual_id))
#93

cph_info <- cph_info %>%
  mutate(sample_intake = year, .after = year)

cph_info <- cph_info %>%
  mutate(OGTT = case_when(year == "2017" ~ OGTT_2017,
                          year == "2020" ~ OGTT_2020,
                          TRUE ~ NA_real_)) %>%
  dplyr::select(-c(OGTT_2017, OGTT_2020))

write.csv(cph_info, "data/formatted/sample_info_cph.csv", row.names = FALSE)



#RPA/NSW
sample_info.rpa <- sample_info %>%
  filter(grepl("^17_", sample_name))
length(unique(sample_info.rpa$individual_id))

mastersheet_info.rpa_sch_14e <- read_xlsx("data/ExoCF mastersheet_RPA Copenhagen SCH.xlsx",
                                          sheet = 2)
mastersheet_info.rpa_sch_14e <- mastersheet_info.rpa_sch_14e %>%
  select(-c(9:23, 25:32, 34, 36, 42:44))
mastersheet_info.rpa_sch_14e <- mastersheet_info.rpa_sch_14e[3:83,]
colnames(mastersheet_info.rpa_sch_14e)[1:2] <- c("diabetes_status", "sample_name")
colnames(mastersheet_info.rpa_sch_14e)[c(7, 9, 10, 11)] <- c("age", "rna_extraction", "sex", "FEV1")
colnames(mastersheet_info.rpa_sch_14e)[c(12, 13)] <- c("mutation1", "mutation2")

for(i in c(2:dim(mastersheet_info.rpa_sch_14e)[1])){
  print(i)
  if(!is.na(mastersheet_info.rpa_sch_14e[i, "sex"]) & 
     mastersheet_info.rpa_sch_14e[i, "sex"] == "See above"){
    mastersheet_info.rpa_sch_14e[i, "sex"] <- mastersheet_info.rpa_sch_14e[(i-1), "sex"]
  }
  if(!is.na(mastersheet_info.rpa_sch_14e[i, "FEV1"]) &
     mastersheet_info.rpa_sch_14e[i, "FEV1"] == "See above"){
    mastersheet_info.rpa_sch_14e[i, "FEV1"] <- mastersheet_info.rpa_sch_14e[(i-1), "FEV1"]
  }
  if(!is.na(mastersheet_info.rpa_sch_14e[i, "mutation1"]) &
     mastersheet_info.rpa_sch_14e[i, "mutation1"] == "See above"){
    mastersheet_info.rpa_sch_14e[i, "mutation1"] <- mastersheet_info.rpa_sch_14e[(i-1), "mutation1"]
  }
  if(!is.na(mastersheet_info.rpa_sch_14e[i, "mutation2"]) &
     mastersheet_info.rpa_sch_14e[i, "mutation2"] == "See above"){
    mastersheet_info.rpa_sch_14e[i, "mutation2"] <- mastersheet_info.rpa_sch_14e[(i-1), "mutation2"]
  }
}
colnames(mastersheet_info.rpa_sch_14e)[6] <- "sample_intake"

mastersheet_info.rpa <- mastersheet_info.rpa_sch_14e %>%
  filter(grepl("^17-", sample_name)) %>%
  select(-c("rna_extraction"))
mastersheet_info.rpa <- mastersheet_info.rpa %>%
  mutate(sample_name = sub("On insulin from SCH 26/4/17", "", sample_name, fixed = TRUE)) %>%
  mutate(sample_name = gsub("-", "_", sample_name, fixed = TRUE)) %>%
  mutate(sample_name = gsub("[[:space:]]", "", sample_name))

mastersheet_info.rpa[c(which(mastersheet_info.rpa$diabetes_status == "CFRD"):
                         (which(mastersheet_info.rpa$diabetes_status == "IGT") - 1)), 
                     "diabetes_status"] <- "CFRD"
mastersheet_info.rpa[c(which(mastersheet_info.rpa$diabetes_status == "IGT"):
                         (which(mastersheet_info.rpa$diabetes_status == "NGT") - 1)), 
                     "diabetes_status"] <- "IGT"
mastersheet_info.rpa[c(which(mastersheet_info.rpa$diabetes_status == "NGT"):
                         dim(mastersheet_info.rpa)[1]), 
                     "diabetes_status"] <- "NGT"

options(digits = 3)
mastersheet_info.rpa <- mastersheet_info.rpa %>%
  mutate(age = as.double(age))

mastersheet_info.rpa <- mastersheet_info.rpa %>%
  separate(FEV1, into = c("FEV1", NA), sep = " ")

# rpa_info <- sample_info.rpa %>%
#   left_join(mastersheet_info.rpa)
# sum(is.na(rpa_info$diabetes_status))
#0

missing_samples_rpa <- mastersheet_info.rpa %>%
  anti_join(sample_info.rpa) %>%
  mutate(sample_name = gsub("_", "-", sample_name, fixed = TRUE))
write.csv(missing_samples_rpa, "data/missing_samples_rpa.csv", row.names = FALSE)

rpa_info <- sample_info.rpa %>%
  inner_join(mastersheet_info.rpa)
sum(is.na(rpa_info$diabetes_status))

rpa_info <- rpa_info %>%
  select(-c(10, 11, 14, 19:21))
colnames(rpa_info)[9] <- "pre_post_modulator"

rpa_info <- rpa_info %>%
  separate(pre_post_modulator, into = c("pre_post_modulator", "modulator"), sep = "\\[") %>%
  mutate(modulator = gsub("]", "", modulator, fixed = TRUE))

#modifying modulators as per mail thread with Sheila and also for sample 17_18_020DH
rpa_info <- rpa_info %>%
  mutate(pre_post_modulator = case_when(sample_name %in% c("17_18_011AS", "17_18_014DM", 
                                                           "17_18_034MC", "17_18_056SB", 
                                                           "17_19_064WC") ~ "0",
                                        TRUE ~ pre_post_modulator)) %>%
  mutate(modulator = case_when(sample_name %in% c("17_18_011AS", "17_18_014DM", 
                                                  "17_18_034MC", "17_18_056SB", 
                                                  "17_19_064WC", "17_18_020DH") ~ NA_character_,
                               TRUE ~ modulator))



#SCH/NSW
sample_info.sch <- sample_info %>%
  filter(grepl("^11_", sample_name))
length(unique(sample_info.sch$individual_id))

mastersheet_info.sch <- mastersheet_info.rpa_sch_14e %>%
  mutate(sample_name = gsub("(SCH)", "", sample_name, fixed = TRUE)) %>%
  mutate(sample_name = gsub("[[:space:]]", "", sample_name)) 

mastersheet_info.sch[c(which(mastersheet_info.sch$diabetes_status == "CFRD"):
                         (which(mastersheet_info.sch$diabetes_status == "IGT") - 1)), 
                     "diabetes_status"] <- "CFRD"
mastersheet_info.sch[c(which(mastersheet_info.sch$diabetes_status == "IGT"):
                         (which(mastersheet_info.sch$diabetes_status == "NGT") - 1)), 
                     "diabetes_status"] <- "IGT"
mastersheet_info.sch[c(which(mastersheet_info.sch$diabetes_status == "NGT"):
                         (which(mastersheet_info.sch$diabetes_status == "HC") - 1)), 
                     "diabetes_status"] <- "NGT"
mastersheet_info.sch[c(which(mastersheet_info.sch$diabetes_status == "HC"):
                         dim(mastersheet_info.sch)[1]), 
                     "diabetes_status"] <- "HC"

mastersheet_info.sch <- mastersheet_info.sch %>%
  filter(grepl("^11-", sample_name)) 

mastersheet_info.sch <- mastersheet_info.sch %>%
  mutate(sample_name = sub("OninsulinfromSCH26/4/17", "", sample_name, fixed = TRUE))

mastersheet_info.sch <- mastersheet_info.sch %>%
  select(c(1:3, 6, 7, 9, 10, 11))
colnames(mastersheet_info.sch)[3] <- "pre_post_modulator"

mastersheet_info.sch <- mastersheet_info.sch %>%
  mutate(sample_date = sample_intake, .before = sample_intake)

mastersheet_info.sch <- mastersheet_info.sch %>%
  mutate(modulator = NA, .before = "age")
mastersheet_info.sch <- mastersheet_info.sch %>%
  mutate(sample_date = gsub("-0", "-", sample_date, fixed = TRUE)) %>%
  mutate(sample_date = gsub("^20", "", sample_date)) %>%
  separate(sample_date, into = c("sd_year", "sd_month", "sd_date"))
mastersheet_info.sch <- mastersheet_info.sch %>%
  mutate(sample_name = paste(sample_name, 
                             sd_date, sd_month, sd_year,
                             sep = "-")
  ) %>%
  select(-c(sd_date, sd_month, sd_year))

mastersheet_info.sch <- mastersheet_info.sch %>%
  relocate(sample_intake, .before = pre_post_modulator)

mastersheet_info.sch2 <- read_xlsx("data/ExoCF mastersheet_RPA Copenhagen SCH.xlsx",
                                  sheet = 4)
mastersheet_info.sch2 <- mastersheet_info.sch2[3:113,]

mastersheet_info.sch2 <- mastersheet_info.sch2 %>%
  mutate(Condition = gsub("CF - Pre + Post modulator", "CF_pre_post_modulator", Condition, fixed = TRUE)) %>%
  mutate(Condition = gsub("CF - non modulator", "CF_non_modulator", Condition, fixed = TRUE)) %>%
  mutate(Condition = gsub("CF BALF", "CF_BALF", Condition, fixed = TRUE)) %>%
  mutate(Condition = gsub("HC BALF", "HC_BALF", Condition, fixed = TRUE)) 

mastersheet_info.sch2[c(which(mastersheet_info.sch2$Condition == "CF_pre_post_modulator"):
                          (which(mastersheet_info.sch2$Condition == "CF_non_modulator") - 1)), 
                      "Condition"] <- "CF_pre_post_modulator"
mastersheet_info.sch2[c(which(mastersheet_info.sch2$Condition == "CF_non_modulator"):
                          (which(mastersheet_info.sch2$Condition == "HC") - 1)), 
                      "Condition"] <- "CF_non_modulator"
mastersheet_info.sch2[c(which(mastersheet_info.sch2$Condition == "HC"):
                          (which(mastersheet_info.sch2$Condition == "CF_BALF") - 1)), 
                      "Condition"] <- "HC"                     
mastersheet_info.sch2[c(which(mastersheet_info.sch2$Condition == "CF_BALF"):
                          (which(mastersheet_info.sch2$Condition == "HC_BALF") - 1)), 
                      "Condition"] <- "CF_BALF"
mastersheet_info.sch2[c(which(mastersheet_info.sch2$Condition == "HC_BALF"):
                          dim(mastersheet_info.sch2)[1]), 
                      "Condition"] <- "HC_BALF"
colnames(mastersheet_info.sch2)[2:4] <- c("sample_name", "pre_post_modulator", "modulator")
mastersheet_info.sch2 <- mastersheet_info.sch2 %>%
  filter(grepl("^11-", sample_name)) 

mastersheet_info.sch2 <- mastersheet_info.sch2 %>%
  filter(!grepl("Plasma|BALF", sample_name))

mastersheet_info.sch2 <- mastersheet_info.sch2 %>%
  mutate(sample_name = gsub(" (serum pink)", "", sample_name, fixed = TRUE))

colnames(mastersheet_info.sch2)[21] <- c("rna_extraction")
mastersheet_info.sch2 <- mastersheet_info.sch2 %>%
  select(c(1:4, 21, 46))

mastersheet_info.sch2 <- mastersheet_info.sch2 %>%
  mutate("age" = NA, .after = "modulator") %>%
  mutate("sex" = NA, .after = "rna_extraction") %>%
  mutate("FEV1" = NA, .after = "sex")
colnames(mastersheet_info.sch)[1] <- "condition"
colnames(mastersheet_info.sch2)[1] <- "condition"

mastersheet_info.sch2 <- mastersheet_info.sch2 %>%
  separate(sample_name, into = c(NA, "sample_intake1"), sep = " ", remove = FALSE) %>%
  separate(sample_name, into = c(NA, "sample_intake2"), sep = "_", remove = FALSE) %>%
  mutate(sample_intake = ifelse(is.na(sample_intake1), sample_intake2, sample_intake1)) %>%
  dplyr::select(-c(sample_intake1, sample_intake2)) %>%
  relocate(sample_intake, .after = sample_name) %>%
  mutate(sample_intake = sub("/20", "/", sample_intake, fixed = TRUE)) %>%
  mutate(sample_intake = gsub("/", ".", sample_intake, fixed = TRUE)) %>%
  mutate(sample_intake = as.Date(sample_intake, format = "%d.%m.%y"))
colnames(mastersheet_info.sch2)[10] <- "OGTT"

mastersheet_info.sch <- rbind(mastersheet_info.sch %>% mutate(OGTT = NA), mastersheet_info.sch2)
mastersheet_info.sch <- mastersheet_info.sch %>%
  mutate(sample_name = gsub("-", "_", sample_name, fixed = TRUE)) %>%
  mutate(sample_name = gsub(" ", "_", sample_name, fixed = TRUE)) %>%
  mutate(sample_name = gsub(".", "_", sample_name, fixed = TRUE)) %>%
  mutate(sample_name = gsub("/", "_", sample_name, fixed = TRUE))


multiple_same_sample_names <- mastersheet_info.sch %>%
  group_by(sample_name) %>%
  summarise(n = n()) %>%
  filter(n > 1)
#4 samples
length(unique(sample_info.sch$sample_name))
#83
length(unique(mastersheet_info.sch$sample_name))
#106

mastersheet_info.sch <- mastersheet_info.sch %>%
  filter(rna_extraction != "-")

multiple_same_sample_names <- mastersheet_info.sch %>%
  group_by(sample_name) %>%
  summarise(n = n()) %>%
  filter(n > 1)
#0 samples

length(unique(mastersheet_info.sch$sample_name))
#89

mastersheet_info.sch <- mastersheet_info.sch %>%
  unique()
#no change

sch_info <- sample_info.sch %>%
  inner_join(mastersheet_info.sch)
sum(is.na(sch_info$condition))


missing_sch <- sample_info.sch %>%
  anti_join(mastersheet_info.sch)
#5 rows

missing_sch2 <- mastersheet_info.sch %>%
  anti_join(sample_info.sch)
#11 rows

sample_info.sch <- sample_info.sch %>%
  mutate(sample_name= case_when(sample_name == "11_16_116OC_WAT8086A8" ~ "11_16_116OC_19_9_17",
                                sample_name == "11_16_162HN_WAT8086A4" ~ "11_16_162HN_15_3_17",
                                sample_name == "11_16_184EH_WAT8086A10" ~ "11_16_184EH_13_3_18",
                                sample_name == "11_16_247AF_WAT8061A9" ~ "11_16_247AF_27_2_18",
                                sample_name == "11_18_403DP_WAT8061A7" ~ "11_18_403DP",
                                TRUE ~ sample_name))

sch_info <- sample_info.sch %>%
  inner_join(mastersheet_info.sch)
sum(is.na(sch_info$condition))
#0

missing_sch <- sample_info.sch %>%
  anti_join(mastersheet_info.sch)
#0 rows

missing_sch2 <- mastersheet_info.sch %>%
  anti_join(sample_info.sch)
#6 rows
write.csv(missing_sch2, "data/missing_samples_sch2.csv", row.names = FALSE)


sch_info <- sch_info %>%
  select(-c(rna_extraction))
sch_info <- sch_info  %>%
  mutate(age = as.double(age))


#get mutation info
mutation_info.sch <- read_xlsx("data/SCH_CFTR mutations.xlsx",
                                   sheet = 1) 
colnames(mutation_info.sch) <- c("sample_id", "mutation1", "mutation2")

mutation_info.sch <- mutation_info.sch %>%
  mutate(sample_id = gsub("-", "_", sample_id, fixed = TRUE))
sch_info <- sch_info %>%
  separate(sample_name, into = c("s1", "s2", "s3", NA, NA, NA), sep = "_", remove = FALSE) %>%
  mutate(sample_id = paste(s1, s2, s3, sep = "_")) %>%
  select(-c(s1, s2, s3))

missing_mutation_info <- sch_info %>%
  anti_join(mutation_info.sch)
write.csv(missing_mutation_info, "data/missing_sch_mutation.csv", row.names = FALSE)

sch_info <- sch_info %>%
  left_join(mutation_info.sch) %>%
  select(-c(sample_id))

write.csv(format(sch_info, digits = 3), "data/formatted/sample_info_sch.csv", row.names = FALSE)


#14E
sample_info.14e <- sample_info %>%
  filter(grepl("^14E_", sample_name))
length(unique(sample_info.14e$individual_id))

mastersheet_info.14e <- mastersheet_info.rpa_sch_14e %>%
  filter(grepl("^14E-", sample_name)) %>%
  mutate(sample_name = gsub("-", "_", sample_name, fixed = TRUE)) %>%
  mutate(diabetes_status = "HC")
mastersheet_info.14e <- mastersheet_info.14e %>%
  select(c(1,2,3,6,7,10,11,12,13))
colnames(mastersheet_info.14e)[3] <- "pre_post_modulator"
mastersheet_info.14e <- mastersheet_info.14e %>%
  mutate(modulator = "", .after = pre_post_modulator)

healthy_info <- sample_info.14e %>%
  inner_join(mastersheet_info.14e)

write.csv(format(healthy_info, digits = 3), "data/formatted/sample_info_14e.csv", row.names = FALSE)

# rpa patients rename
rpa_info <- rpa_info %>%
  mutate(patient_initial = str_extract(individual_id, "[^0-9]+$"), .after = individual_id)

rpa_info_temp <- rpa_info %>%
  group_by(patient_initial) %>%
  summarise(n = n()) %>%
  filter(n > 1) %>%
  arrange(desc(n))

rpa_multiple <- rpa_info %>%
  select(sample_name, patient_initial) %>%
  filter(patient_initial %in% c(rpa_info_temp$patient_initial)) %>%
  arrange(patient_initial)

#AS and SD samples are from different patients 
rpa_multiple <- rpa_multiple %>%
  filter(!patient_initial %in% c("AS", "SD"))

rpa_info <- rpa_info %>%
  mutate(individual_id = case_when(sample_name %in% rpa_multiple$sample_name ~ patient_initial,
                                   TRUE ~ individual_id)) %>%
  arrange(individual_id) %>%
  select(-c(patient_initial))

write.csv(format(rpa_info, digits = 3), "data/formatted/sample_info_rpa.csv", row.names = FALSE)
############# combine cohorts

colnames(cph_info)
colnames(rpa_info)
colnames(sch_info)
colnames(healthy_info)

cph_info <- cph_info %>%
  relocate(age, .before = sex) %>%
  relocate(OGTT, .after = FEV1)
rpa_info <- rpa_info %>%
  relocate(sample_intake, .after = year)
sch_info <- sch_info %>%
  relocate(sample_intake, .after = year)
healthy_info <- healthy_info %>%
  relocate(sample_intake, .after = year)

rpa_info <- rpa_info %>%
  mutate(OGTT = NA_real_, .after = FEV1)
healthy_info <- healthy_info %>%
  mutate(OGTT = NA_real_, .after = FEV1)

colnames(cph_info)
colnames(rpa_info)
colnames(sch_info)
colnames(healthy_info)

colnames(cph_info) <- colnames(sch_info)
colnames(rpa_info) <- colnames(sch_info)
colnames(healthy_info) <- colnames(sch_info)

colnames(cph_info)
colnames(rpa_info)
colnames(sch_info)
colnames(healthy_info)

meta_data_au <- rbind(rpa_info %>%
                        mutate("cohort" = "RPA_NSW", .after = year) %>%
                        mutate("country" = "AU", .after = cohort),
                      sch_info %>%
                        mutate("cohort" = "SCH_NSW", .after = year) %>%
                        mutate("country" = "AU", .after = cohort), 
                      healthy_info %>%
                        mutate("cohort" = "UNSW", .after = year) %>%
                        mutate("country" = "AU", .after = cohort)) %>%
  rename(c("sample_intake_date" = "sample_intake")) %>%
  mutate(sample_intake_year = format(sample_intake_date, format = "%Y"), .after = sample_intake_date)

meta_data <- rbind(cph_info %>%
                     mutate("cohort" = "CPH", .after = year) %>%
                     mutate("country" = "DK", .after = cohort) %>%
                     rename(c("sample_intake_year" = "sample_intake")) %>%
                     mutate(sample_intake_date = as.Date(""), .before = sample_intake_year),
                   meta_data_au
                   )

meta_data <- meta_data %>%
  rename(c("patient_recruitment_year" = "year"))

levels(factor(meta_data$modulator))
levels(factor(meta_data$pre_post_modulator))

meta_data <- meta_data %>%
  mutate(modulator = case_when(modulator %in% c("", "-", "No CFTR-modulator") ~ NA_character_,
                               TRUE ~ modulator)) %>%
  mutate(pre_post_modulator = gsub(" ", "", pre_post_modulator)) %>%
  mutate(pre_post_modulator = case_when(pre_post_modulator == "-" ~ NA_character_,
                                        TRUE ~ pre_post_modulator)
         ) %>%
  mutate(pre_post_modulator = strtoi(pre_post_modulator))

levels(factor(meta_data$modulator))
levels(factor(meta_data$pre_post_modulator))

summary(meta_data$pre_post_modulator)
summary(factor(meta_data$pre_post_modulator))


summary(factor(meta_data$age))
summary(factor(meta_data$sex))
summary(factor(meta_data$FEV1))

meta_data <- meta_data %>%
  mutate(age = case_when(age == "-" ~ NA_character_,
                         TRUE ~ age)) %>%
  mutate(sex = case_when(sex == "-" ~ NA_character_,
                         TRUE ~ sex)) %>%
  mutate(FEV1 = case_when(FEV1 == "-" ~ NA_character_,
                          TRUE ~ FEV1)) %>%
  mutate(mutation1 = case_when(mutation1 == "-" ~ NA_character_,
                          TRUE ~ mutation1)) %>%
  mutate(mutation2 = case_when(mutation2 == "-" ~ NA_character_,
                               TRUE ~ mutation2))
options(digits = 3)
meta_data <- meta_data %>%
  mutate(age = as.double(age)) %>%
  mutate(FEV1 = as.double(FEV1))


### add library quality
prep_outcome <- read_xlsx("data/LibraryPrepOutcome-WAT10166-WAT10180.xlsx")[, 1:4]
colnames(prep_outcome) <- c("seq_plate", "lims_id", "sample_name", "seq_miR_library_quality")

prep_outcome <- prep_outcome %>%
  mutate(sample_name = gsub("/", "_", sample_name, fixed = TRUE)) %>%
  mutate(sample_name = gsub("-", "_", sample_name, fixed = TRUE)) %>%
  mutate(sample_name = gsub(" ", "_", sample_name, fixed = TRUE)) %>%
  mutate(sample_name = gsub(".", "_", sample_name, fixed = TRUE)) 

meta_data_with_qual <- meta_data %>%
  separate(sample_long_name, into = c("sample_name_2", NA), sep = "_S", remove = FALSE)


missing1 <- meta_data_with_qual %>%
  anti_join(prep_outcome, by = c("sample_name_2" = "sample_name"))
missing2 <-  prep_outcome %>%
  anti_join(meta_data_with_qual, by = c("sample_name" = "sample_name_2"))


meta_data_with_qual <- meta_data_with_qual %>%
  left_join(prep_outcome %>%
              select(-c(lims_id)),
            by = c("sample_name_2" = "sample_name")) %>%
  select(-c(sample_name_2))

meta_data_with_qual <- meta_data_with_qual %>%
  mutate(age_group = case_when(cohort == "SCH_NSW" ~ "child",
                               TRUE ~ "adult"), .after = age)

meta_data_with_qual <- meta_data_with_qual %>%
  mutate(modulator = case_when(pre_post_modulator == 0 ~ NA_character_,
                               TRUE ~ modulator))
levels(factor(meta_data_with_qual$modulator))

#same drugs named differently in AU and EU
meta_data_with_qual <- meta_data_with_qual %>%
  mutate(modulator = gsub("Kaftrio", "Trikafta", modulator, fixed = TRUE)) %>%
  mutate(modulator = gsub("Symkevi", "Symdeko", modulator, fixed = TRUE))

levels(factor(meta_data_with_qual$modulator))

summary(factor(meta_data_with_qual$seq_miR_library_quality))
meta_data_with_qual <- meta_data_with_qual %>%
  mutate(seq_miR_library_quality = case_when(sample_name %in% c("CPH10", "CPH22") ~ "Good",
                                             TRUE ~ seq_miR_library_quality))
summary(factor(meta_data_with_qual$seq_miR_library_quality))


summary(factor(meta_data_with_qual$mutation1))
summary(factor(meta_data_with_qual$mutation2))

meta_data_with_qual[meta_data_with_qual["individual_id"] == "115CS", "mutation1"] <- "N1303K"
meta_data_with_qual[meta_data_with_qual["individual_id"] == "115CS", "mutation2"] <- "F508del"   

meta_data_with_qual <- meta_data_with_qual %>%
  mutate(mutation1 = gsub("Andet", "unknown", mutation1)) %>%
  mutate(mutation2 = gsub("Andet", "unknown", mutation2)) %>%
  mutate(mutation = case_when(is.na(mutation1) ~ NA_character_,
                              mutation1 == "F508del" ~ paste0("F508del___", mutation2),
                              mutation2 == "F508del" ~ paste0("F508del___", mutation1),
                              mutation1 < mutation2 ~ paste(mutation1, mutation2, sep = "___"),
                              TRUE ~ paste(mutation2, mutation1, sep = "___")))
meta_data_with_qual <- meta_data_with_qual %>%
  select(-c(mutation1, mutation2))

summa <- data.frame(summary(factor(meta_data_with_qual$mutation)))
colnames(summa) <- "count"
summa <- summa %>%
  arrange(desc(count))

write.csv(summa, "data/mutation_count_summary.csv")

write.csv(format(meta_data_with_qual, digits = 3), "data/formatted/meta_data.csv", row.names = FALSE)


#check with umi counts
umi_counts <- read.csv("data/formatted/umi_counts.csv", row.names = 1)
colnames(umi_counts) <- gsub("^X", "", colnames(umi_counts))

umi_samples <- data.frame(sample_long_name = colnames(umi_counts))

check_combination <- umi_samples %>%
  inner_join(meta_data_with_qual)
#272

missing1 <- umi_samples %>%
  anti_join(meta_data_with_qual)
#0

missing2 <- meta_data_with_qual %>%
  anti_join(umi_samples)
#0


###################### including CFRD/IGT/NGT status, OGTT for some child samples from info provided later by Bernadette

results <- read_excel("prediction_pipeline/sch_pred_with_clinical_results.xlsx")
results <- results[-c(1), c(2, 6, 8, 9, 10, 11)]
colnames(results) <- c("Sample", "pre_post_mod", "prediction", "matching_dates", "OGTT", "actual")
new_label_info <- results %>%
  filter(matching_dates == "P") %>%
  filter(!is.na(actual)) %>%
  dplyr::select(c(Sample, actual, pre_post_mod, OGTT)) %>%
  rename(c("sample_long_name" = "Sample"))

meta_data_with_qual <- meta_data_with_qual %>%
  left_join(new_label_info, by = "sample_long_name")
meta_data_with_qual <- meta_data_with_qual %>%
  mutate(condition = case_when(!is.na(actual) ~ actual,
                               TRUE ~ condition))

meta_data_with_qual %>%
  filter(!is.na(pre_post_mod) & !is.na(pre_post_modulator) & pre_post_mod != pre_post_modulator)
#0 rows

meta_data_with_qual %>%
  filter(!is.na(OGTT.x) & !is.na(OGTT.y) & OGTT.x != OGTT.y)
#0 rows

meta_data_with_qual <- meta_data_with_qual %>%
  mutate(OGTT = ifelse(!is.na(OGTT.x), OGTT.x, OGTT.y)) %>%
  dplyr::select(-c(actual, pre_post_mod, OGTT.x, OGTT.y))

#OGTT values that match with sample date from sheet shared by Laura 
#only 116OC sample OGTT values are present in sheet shared by Laura thats not present in master sheet
meta_data_with_qual[meta_data_with_qual$sample_name == "11_16_116OC_16_8_16", "OGTT"] <- 5.3
meta_data_with_qual[meta_data_with_qual$sample_name == "11_16_116OC_30_4_13", "OGTT"] <- 5.3
meta_data_with_qual[meta_data_with_qual$sample_name == "11_16_116OC_8_5_12", "OGTT"] <- 5.4
meta_data_with_qual[meta_data_with_qual$sample_name == "11_16_116OC_19_9_17", "OGTT"] <- 4.1

write.csv(format(meta_data_with_qual, digits = 3), "data/formatted/meta_data.csv", row.names = FALSE)


###########

#compare condition from OGTT and condition field already present
meta_data_updated_condition <- meta_data_with_qual %>%
  mutate(OGTT = as.numeric(OGTT),
         condition_from_OGTT = case_when(is.na(OGTT) ~ NA_character_,
                                         OGTT >= 11.1 ~ "CFRD",
                                         OGTT < 7.8 ~ "NGT",
                                         TRUE ~ "IGT"),
         condition_updated = case_when(is.na(condition) ~ condition_from_OGTT,
                                       condition != "IND" & condition != "CF_pre_post_modulator" & !is.na(condition_from_OGTT) ~ condition_from_OGTT,
                                       TRUE ~ condition))
write.csv(format(meta_data_updated_condition, digits = 3), "data/formatted/meta_data.csv", row.names = FALSE)

mismatch_entries <- meta_data_updated_condition %>%
  filter(!is.na(condition_from_OGTT) & condition_from_OGTT != condition)
write.csv(format(mismatch_entries, digits = 3), "data/formatted/meta_data_condition_mismatch.csv", row.names = FALSE)


############################################

meta_data <- read.csv("data/formatted/meta_data_modified_2023June28.csv")
ogtt_bernadatette <- read.csv("data/formatted/ogtt_au.csv")
colnames(ogtt_bernadatette) <- c("sample", 
                                 "pre_post_modulator_b", 
                                 "are_dates_matching_b", 
                                 "sample_intake_date_b",
                                 "OGTT_date_b", "OGTT_2h_b",
                                 "condition_b", "comment", "condition_other")
meta_data_updated <- meta_data %>%
  left_join(ogtt_bernadatette %>% dplyr::select(-c(condition_b)), by = c("sample_long_name" = "sample"))

# subset <- meta_data_updated %>%
#   filter(!is.na(sample_intake_date_b))
# 
# subset %>%
#   filter(sample_intake_date != sample_intake_date_b)
# subset %>%
#   filter(pre_post_modulator != pre_post_modulator_b)
# s <- subset %>%
#   filter(OGTT_2h != OGTT_2h_b) %>%
#   dplyr::select(c(sample_long_name, OGTT_2h, OGTT_2h_b)) %>%
#   filter(OGTT_2h != "NA")
# sum(is.na(s$OGTT_2h))


meta_data_updated <- meta_data_updated %>%
  dplyr::select(-c(note, sample_intake_date_b, pre_post_modulator_b, ogtt_date, are_dates_matching_b)) %>%
  rename(c("OGTT_date" = "OGTT_date_b"))
write.csv(meta_data_updated, "data/formatted/meta_data_with_info_from_bernadette.csv", row.names = FALSE)



########################

#create proteomics meta-data

#copy meta-data from transcriptomics for samples available
#others use proteomics meta-data
tra_meta_data <- read.csv("data/formatted/meta_data_updated.csv") %>%
  rename(c("condition" = "condition_updated")) %>%
  rename("Sample" = "sample_long_name") %>%
  mutate(Sample = paste0("X", Sample))

prot_meta_data_from_Alex <- read.table("data/proteomics/AU_DK/proteomics_sinfo_modified.txt", 
                                           sep = "\t", header = TRUE)
prot_summary_others <- read.table("data/proteomics/AU_DK/DK_others/summary.txt", 
                                  sep = "\t", header = TRUE)[, c(1:2)]
colnames(prot_summary_others) <- c("rawfile", "label")

prot_meta_data_to_be_used <- rbind(prot_meta_data_from_Alex %>%
  dplyr::select(rawfile, label, technicalreplicate),
  prot_summary_others %>%
    filter(rawfile != "Total") %>%
    mutate(technicalreplicate = NA))

mapping <- read_excel("data/proteomics/pheno_combined_new_formatted.xlsx", sheet = "t_p_mapping")
mapping.t_avlbl <- mapping %>%
  filter(!is.na(sample_long_name_t) & !is.na(sample_long_name_p)) %>%
  dplyr::select(c(sample_long_name_p, sample_long_name_t))
mapping.t_not_avlbl <- mapping %>%
  filter(is.na(sample_long_name_t) & !is.na(sample_long_name_p))

length(unique(tra_meta_data$Sample))
# [1] 272
length(unique(mapping.t_avlbl$sample_long_name_t))
# [1] 238
length(unique(mapping.t_avlbl$sample_long_name_p))
# [1] 249
length(unique(mapping.t_not_avlbl$sample_long_name_p))
# [1] 21
length(unique(mapping$sample_long_name_p))
# [1] 271
# this is 249 + 21 + 1(NA)

prot_meta_data1 <- tra_meta_data %>%
  inner_join(mapping.t_avlbl, by = c("Sample" = "sample_long_name_t")) %>%
  inner_join(prot_meta_data_to_be_used, by = c("sample_long_name_p" = "label"))
#259
length(unique(prot_meta_data1$Sample))
# [1] 236
write.csv(prot_meta_data1, 
          "data/proteomics/prot_metadata1.csv", row.names = FALSE)
#created a new file data/proteomics/prot_metadata1_modified.csv by 
# manually entering values for missing technicalreplicate


prot_meta_data_from_Alex_sub <- read.table("data/proteomics/AU_DK/proteomics_sinfo_modified.txt", 
                                       sep = "\t", header = TRUE) %>%
  inner_join(mapping.t_not_avlbl, by = c("label" = "sample_long_name_p"))
write.csv(prot_meta_data_from_Alex_sub, "data/proteomics/prot_metadata_missing_in_tra.csv", row.names = FALSE)

#manually modifying the above file and creating a new file
# data/proteomics/prot_metadata_missing_in_tra_modified.csv



#now just rely on 2 files for proteomics metadata
#data/proteomics/prot_metadata1_modified.csv
#data/proteomics/prot_metadata_missing_in_tra_modified.csv

prot_meta_data1 <- read.csv("data/proteomics/prot_metadata1_modified.csv")
prot_meta_data_missing_tra <- read.csv("data/proteomics/prot_metadata_missing_in_tra_modified.csv") %>%
  dplyr::select(c(rawfile, label, individual_id, country, age_group, technicalreplicate,
                  biotype, sample_intake_year, FEV1, cohort, age, sample_intake_date, 
                  pre_post_modulator, modulator, sex, mutation1, mutation2,
                  OGTT_1h, OGTT_2h, condition_updated)) %>%
  mutate(mutation = case_when(is.na(mutation1) ~ NA_character_,
                              mutation1 == "F508del" ~ paste0("F508del___", mutation2),
                              mutation2 == "F508del" ~ paste0("F508del___", mutation1),
                              mutation1 < mutation2 ~ paste(mutation1, mutation2, sep = "___"),
                              TRUE ~ paste(mutation2, mutation1, sep = "___"))) %>%
  dplyr::select(-c(mutation1, mutation2)) %>%
  rename(c("condition" = "condition_updated"))

prot_meta_data1 <- prot_meta_data1 %>%
  rename(c("sample_long_name_t" = "Sample", "label" = "sample_long_name_p")) %>%
  dplyr::select(c(label, rawfile, technicalreplicate, 
                  sample_long_name_t, condition, individual_id, 
                  sample_name, cohort, country, 
                  sample_intake_date, sample_intake_year, 
                  pre_post_modulator, modulator, 
                  age, age_group, sex, FEV1, mutation, 
                  illumina_sample_number, quant_batch, biotype, 
                  patient_recruitment_year, 
                  seq_plate, seq_miR_library_quality, 
                  condition_from_OGTT, condition_initially_used, 
                  OGTT_1h, OGTT_2h,
                  OGTT_date, OGTT_2h_b, comment, condition_other))

prot_meta_data_missing_tra <- prot_meta_data_missing_tra %>%
  mutate(sample_long_name_t = NA,
           sample_name = NA,
           illumina_sample_number = NA,
           quant_batch = NA, 
           patient_recruitment_year = NA, 
           seq_plate = NA, 
           seq_miR_library_quality = NA, 
           condition_from_OGTT = NA,
           condition_initially_used = NA, 
           OGTT_date = NA,
           OGTT_2h_b = NA,
           comment = NA,
           condition_other = NA) %>%
  dplyr::select(c(label, rawfile, technicalreplicate, 
                  sample_long_name_t, condition, individual_id, 
                  sample_name, cohort, country, 
                  sample_intake_date, sample_intake_year, 
                  pre_post_modulator, modulator, 
                  age, age_group, sex, FEV1, mutation, 
                  illumina_sample_number, quant_batch, biotype, 
                  patient_recruitment_year, 
                  seq_plate, seq_miR_library_quality, 
                  condition_from_OGTT, condition_initially_used, 
                  OGTT_1h, OGTT_2h,
                  OGTT_date, OGTT_2h_b, comment, condition_other))

all.equal(colnames(prot_meta_data1), colnames(prot_meta_data_missing_tra))

prot_meta_data <- rbind(prot_meta_data1, prot_meta_data_missing_tra)

sum(prot_meta_data$rawfile %in% prot_meta_data_from_Alex$rawfile)
#262
sum(prot_meta_data$rawfile %in% prot_summary_others$rawfile)
#18

length(prot_meta_data_from_Alex$rawfile)
length(unique(prot_meta_data_from_Alex$rawfile))

sum(!prot_meta_data_from_Alex$rawfile %in% prot_meta_data$rawfile)
#16
dim(prot_meta_data_from_Alex[! prot_meta_data_from_Alex$rawfile %in% prot_meta_data$rawfile, 
               c('rawfile', 'label', 'condition')])
#16 3
missing_data <- prot_meta_data_from_Alex[! prot_meta_data_from_Alex$rawfile %in% prot_meta_data$rawfile, 
                                         c('rawfile', 'label', 'condition')]
#all glufib, QC
#and
# SW-exo-6-9-21-57-90min SW.exo.6.9.21.57.90min IGT
# SW-exo-25-8-21-13-90min SW.exo.25.8.21.13.90min NGT
#just these 2 were not prepended with 'X' in mapping file

missing_ones_required <- mapping.t_avlbl %>%
  filter(sample_long_name_p %in% c('SW.exo.6.9.21.57.90min', 'SW.exo.25.8.21.13.90min')) %>%
  mutate(sample_long_name_t = paste0('X', sample_long_name_t)) %>%
  inner_join(tra_meta_data, by = c("sample_long_name_t" = "Sample")) %>%
  inner_join(prot_meta_data_to_be_used, by = c("sample_long_name_p" = "label")) %>%
  rename(c("label" = "sample_long_name_p")) %>%
  dplyr::select(c(label, rawfile, technicalreplicate, 
                  sample_long_name_t, condition, individual_id, 
                  sample_name, cohort, country, 
                  sample_intake_date, sample_intake_year, 
                  pre_post_modulator, modulator, 
                  age, age_group, sex, FEV1, mutation, 
                  illumina_sample_number, quant_batch, biotype, 
                  patient_recruitment_year, 
                  seq_plate, seq_miR_library_quality, 
                  condition_from_OGTT, condition_initially_used, 
                  OGTT_1h, OGTT_2h,
                  OGTT_date, OGTT_2h_b, comment, condition_other))
all.equal(colnames(prot_meta_data), colnames(missing_ones_required))
prot_meta_data <- rbind(prot_meta_data, missing_ones_required)

sum(prot_meta_data$rawfile %in% prot_meta_data_from_Alex$rawfile)
#264
sum(prot_meta_data$rawfile %in% prot_summary_others$rawfile)
#18

length(prot_meta_data_from_Alex$rawfile)
length(unique(prot_meta_data_from_Alex$rawfile))

sum(!prot_meta_data_from_Alex$rawfile %in% prot_meta_data$rawfile)
#14
dim(prot_meta_data_from_Alex[! prot_meta_data_from_Alex$rawfile %in% prot_meta_data$rawfile, 
                             c('rawfile', 'label', 'condition')])
#14 3
missing_data <- prot_meta_data_from_Alex[! prot_meta_data_from_Alex$rawfile %in% prot_meta_data$rawfile, 
                                         c('rawfile', 'label', 'condition')]


dim(prot_meta_data[prot_meta_data$rawfile %in% prot_summary_others$rawfile, c('rawfile', 'label', 'condition')])
prot_meta_data[prot_meta_data$rawfile %in% prot_summary_others$rawfile, c('rawfile', 'label', 'condition')]

prot_meta_data <- prot_meta_data %>%
  mutate(mq_batch = case_when(prot_meta_data$rawfile %in% prot_summary_others$rawfile ~ 'other',
                              TRUE ~ 'main'))
summary(factor(prot_meta_data$mq_batch))

write.csv(prot_meta_data, 
          "data/proteomics/prot_metadata_all.csv", row.names = FALSE)


#identifying duplicates incorrectly created in proteomics meta-data
prot_meta_data <- read.csv("data/proteomics/prot_metadata_all.csv") 
length(unique(prot_meta_data$label))

non_unique_label <- prot_meta_data %>%
  group_by(label) %>%
  summarize(count = n()) %>%
  filter(count > 1)

duplicates <- prot_meta_data %>%
  filter(label %in% non_unique_label$label) %>%
  dplyr::select(c(label, rawfile, technicalreplicate, sample_long_name_t, 
                  individual_id, sample_intake_date, sample_intake_year)) %>%
  arrange(label)
write.csv(duplicates, "data/proteomics/incorrectly_added_duplicates.csv", row.names = FALSE)



#special case for 179AH sample 105 - corresponds to X11_16_179AH_28_03_2017_S227
tra_meta_data_179AH_2017 <- read.csv("data/formatted/meta_data_updated.csv") %>%
  rename(c("condition" = "condition_updated")) %>%
  rename("Sample" = "sample_long_name") %>%
  mutate(Sample = paste0("X", Sample)) %>%
  filter(Sample == "X11_16_179AH_28_03_2017_S227") %>%
  mutate(sample_long_name_p = "SW.exo.16.9.21.105.90min", rawfile = "SW-exo-16-9-21-105-90min",
         technicalreplicate = "t1", mq_batch = "main") %>%
  rename(c("sample_long_name_t" = "Sample", "label" = "sample_long_name_p")) %>%
  dplyr::select(c(label, rawfile, technicalreplicate, 
                  sample_long_name_t, condition, individual_id, 
                  sample_name, cohort, country, 
                  sample_intake_date, sample_intake_year, 
                  pre_post_modulator, modulator, 
                  age, age_group, sex, FEV1, mutation, 
                  illumina_sample_number, quant_batch, biotype, 
                  patient_recruitment_year, 
                  seq_plate, seq_miR_library_quality, 
                  condition_from_OGTT, condition_initially_used, 
                  OGTT_1h, OGTT_2h,
                  OGTT_date, OGTT_2h_b, comment, condition_other, mq_batch))

prot_meta_data_corrected <- read.csv("data/proteomics/prot_metadata_all_corrected.csv") 

all.equal(colnames(prot_meta_data_corrected), colnames(tra_meta_data_179AH_2017))


write.csv(tra_meta_data_179AH_2017, "data/proteomics/179AH_105_entry.csv", row.names = FALSE)
#manually copied this row into data/proteomics/prot_metadata_all_corrected.csv

#also manually corrected technicalreplicate column in data/proteomics/prot_metadata_all_corrected.csv
#   so that entries marked 't2' have corresponding 
#  't1' with same tra mapping / if no tra mapping, then same sample date


#found that for mq_batch = 'other', labels contain -
# should be updated to .
prot_meta_data_corrected <- read.csv("data/proteomics/prot_metadata_all_corrected.csv") %>%
  mutate(label = case_when(mq_batch == "other" ~ gsub("-", ".", label, fixed = TRUE),
                           TRUE ~ label)) %>%
  arrange(cohort, sample_long_name_t, label)
write.csv(prot_meta_data_corrected, "data/proteomics/prot_metadata_updated.csv", row.names = FALSE)

#main proteomics metadata file : "data/proteomics/prot_metadata_updated.csv"           



prot_meta_data <- read.csv("data/proteomics/prot_metadata_updated.csv") 
length(unique(prot_meta_data$label))
length(unique(prot_meta_data$rawfile))

summary(factor(prot_meta_data$mq_batch))
sum(is.na(prot_meta_data$sample_long_name_t))

prot_meta_data_notra <- prot_meta_data %>%
  filter(is.na(sample_long_name_t))
prot_meta_data_tra <- prot_meta_data %>%
  filter(!is.na(sample_long_name_t))

summary(factor(prot_meta_data_notra$mq_batch))
summary(factor(prot_meta_data_tra$mq_batch))

prot_meta_data <- read.csv("data/proteomics/prot_metadata_updated.csv") %>%
  dplyr::select(c(label, rawfile, sample_long_name_t, condition))
tra_meta_data <- read.csv("data/formatted/meta_data_updated.csv") %>%
  mutate(sample_long_name_t = paste0("X", sample_long_name)) %>%
  dplyr::select(c(sample_long_name_t, condition_updated)) %>%
  rename(c("condition" = "condition_updated"))
updated_mapping <- tra_meta_data %>%
  full_join(prot_meta_data)
length(unique(updated_mapping$sample_long_name_t))
length(unique(updated_mapping$label))
sum(is.na(updated_mapping$sample_long_name_t))
#21
sum(is.na(updated_mapping$label))
#33
sum(is.na(updated_mapping$rawfile))
#33
write.csv(updated_mapping, "data/formatted/updated_mapping.csv", row.names = FALSE)



#creating/including meta data of new samples
summary_info <- read.table("data/proteomics/2023_samples/2_MaxQuant_Processed/exec1_with_default_parameters/summary.txt", 
                           header = TRUE, sep = "\t") 

sample_info <- summary_info %>%
  select(Raw.file, Experiment) %>%
  filter(Raw.file != "Total")
colnames(sample_info) <- c("rawfile", "label")

summary_info <- read.table("data/proteomics/2023_samples/2_MaxQuant_Processed/exec2_with_modified_parameters/summary.txt", 
                           header = TRUE, sep = "\t") 

sample_info2 <- summary_info %>%
  select(Raw.file, Experiment) %>%
  filter(Raw.file != "Total")
colnames(sample_info2) <- c("rawfile", "label")

all.equal(sample_info, sample_info2)

library(XLConnect)
#add password before running the below line and do not commit the password
# wb <- loadWorkbook("data/CFdiabF-ShipmentDataAll2022.xlsx", password="")
data_from_bibi <- readWorksheet(wb, sheet = 1)[-c(1), ] %>%
  dplyr::select(c(record_id, age, date_visit_22, record_id_22, ogtt_0h_22, ogtt_1h_22, ogtt_2h_22,
                  diabstatus_22, early_impairment_22)) %>%
  filter(!is.na(record_id_22)) %>%
  mutate(record_id = as.integer(record_id)) %>%
  mutate(record_id_22 = as.integer(record_id_22)) %>%
  arrange(record_id_22)

data_from_bibi <- data_from_bibi  %>%
  mutate(ogtt_2h_22 = as.numeric(str_squish(ogtt_2h_22))) %>%
  mutate(diabstatus_22 = as.numeric(str_squish(diabstatus_22))) %>%
  dplyr::mutate(disease_status = case_when(is.na(ogtt_2h_22) ~ NA_character_,
                                                   ogtt_2h_22 >= 11.1 ~ "CFRD",
                                                   ogtt_2h_22 < 7.8 ~ "NGT",
                                                   TRUE ~ "IGT")) %>%
  dplyr::mutate(disease_status_directly_available = case_when(is.na(diabstatus_22) ~ NA_character_,
                                                         diabstatus_22 == 1 ~ "NGT",
                                                         diabstatus_22 == 2 ~ "Indet",
                                                         diabstatus_22 == 3 ~ "IGT",
                                                         diabstatus_22 == 4 ~ "CFRD")) %>%
  mutate(modulator_status = "postmod")

sample_info <- sample_info %>%
  separate(rawfile, into = c(NA, NA, NA, NA, "record_id_22", NA), remove = FALSE) %>%
  mutate(record_id_22 = as.integer(record_id_22))

new_sample_metadata <- sample_info %>%
  inner_join(data_from_bibi) %>%
  arrange(record_id_22)
write.csv(new_sample_metadata, 
          "data/proteomics/prot_metadata_new_samples.csv", 
          row.names = FALSE)



#check if raw files copied to external hard drive new location matches combined meta-data

meta_data1 <- read.csv("data/proteomics/prot_metadata_updated.csv") %>%
  dplyr::select(c(label, rawfile))
meta_data2 <- read.csv("data/proteomics/prot_metadata_new_samples.csv") %>%
  dplyr::select(c(label, rawfile))

combined_meta_data <- rbind(meta_data1, meta_data2) %>%
  filter(label != "21-7-23-31-90min") %>%
  mutate(rawfile = paste0(rawfile, ".raw")) %>%
  arrange(rawfile)
file_list <- read.csv("data/proteomics/exoCF_rawfile_names.csv", header = FALSE, col.names = c("rawfile")) %>%
  arrange(rawfile)

all.equal(combined_meta_data$rawfile, file_list$rawfile)
#TRUE




#create meta-data file for proteomics new mq_analysis results with all samples (333 samples)
summary_info <- read.table("data/proteomics/all/summary.txt", header = TRUE, sep = "\t") %>%
  select(Raw.file, Experiment) %>%
  filter(Raw.file != "Total")
colnames(summary_info) <- c("rawfile", "label")
# summary_info <- summary_info %>%
#   mutate(label_prev_analysis = paste0("S", label))

meta_data1 <- read.csv("data/proteomics/prot_metadata_updated.csv") %>%
  dplyr::rename(c("label_old" = "label")) 

meta_data1 <- summary_info %>%
  inner_join(meta_data1, by = "rawfile")

meta_data1 <- meta_data1 %>%
  dplyr::select(-c(label_old)) %>%
  dplyr::rename(c("batch_name" = "mq_batch"))

#get other samples - i.e. the new ones
meta_data2 <- summary_info %>%
  anti_join(meta_data1, by = "rawfile")


library(XLConnect)
#add password before running the below line and do not commit the password
# wb <- loadWorkbook("data/CFdiabF-ShipmentDataAll2022.xlsx", password="")
data_from_bibi <- readWorksheet(wb, sheet = 1)[-c(1), ]
data_from_bibi <- data_from_bibi %>%
  mutate(cf_mutation1 = case_when(cf_mutation1 == "Andet" ~ cf_mutation1other,
                                  TRUE ~ cf_mutation1),
         cf_mutation2 = case_when(cf_mutation2 == "Andet" ~ cf_mutation2other,
                                  TRUE ~ cf_mutation2)) %>%
  mutate(mutation = case_when(cf_mutation1 == "F508del" ~ paste0("F508del___", cf_mutation2),
                              cf_mutation2 == "F508del" ~ paste0("F508del___", cf_mutation1),
                              cf_mutation1 < cf_mutation2 ~ paste(cf_mutation1, cf_mutation2, sep = "___"),
                              TRUE ~ paste(cf_mutation2, cf_mutation1, sep = "___"))) %>%
  dplyr::select(c(record_id, age, date_visit_22, record_id_22, ogtt_0h_22, ogtt_1h_22, ogtt_2h_22,
                  diabstatus_22, early_impairment_22, mutation, age, sex, fev1_pp)) %>%
  mutate(MF = case_when(sex == 1 ~ "M",
                        sex == 2 ~ "F",
                        TRUE ~ NA_character_)) %>%
  mutate(sex = MF) %>%
  dplyr::select(-c(MF)) %>%
  filter(!is.na(record_id_22)) %>%
  mutate(record_id = as.integer(record_id)) %>%
  mutate(record_id_22 = as.integer(record_id_22)) %>%
  arrange(record_id_22)

data_from_bibi <- data_from_bibi  %>%
  mutate(ogtt_2h_22 = as.numeric(str_squish(ogtt_2h_22))) %>%
  mutate(diabstatus_22 = as.numeric(str_squish(diabstatus_22))) %>%
  dplyr::mutate(disease_status = case_when(is.na(ogtt_2h_22) ~ NA_character_,
                                           ogtt_2h_22 >= 11.1 ~ "CFRD",
                                           ogtt_2h_22 < 7.8 ~ "NGT",
                                           TRUE ~ "IGT")) %>%
  dplyr::mutate(disease_status_directly_available = case_when(is.na(diabstatus_22) ~ NA_character_,
                                                              diabstatus_22 == 1 ~ "NGT",
                                                              diabstatus_22 == 2 ~ "Indet",
                                                              diabstatus_22 == 3 ~ "IGT",
                                                              diabstatus_22 == 4 ~ "CFRD")) %>%
  mutate(pre_post_modulator = 1)

meta_data2 <- meta_data2 %>%
  separate(rawfile, into = c(NA, NA, NA, NA, "record_id_22", NA), remove = FALSE) %>%
  mutate(record_id_22 = as.integer(record_id_22))

new_sample_metadata <- meta_data2 %>%
  inner_join(data_from_bibi, by = "record_id_22") %>%
  arrange(record_id_22)

#taking record_id as the id that follows CPH
# i.e. record_id = 42, record_id_22 = 1   -----> CPH42
# manually checked the corresponding FEV values for few entries and they match

#The FEV value and age is from previous intake
#cant calculate current age since previous sample intake date is not available 
#to subtract from new sample intake date

new_sample_metadata <- new_sample_metadata %>%
  mutate(individual_id = paste0("CPH", record_id)) %>%
  dplyr::select(-c(record_id_22, age, early_impairment_22, 
                   fev1_pp, record_id, diabstatus_22)) %>%
  mutate(sample_intake_date = as.Date(date_visit_22, format = "%Y-%M-%d %H:%M:%S"),
         .after = "date_visit_22") %>%
  mutate(sample_intake_year = format(sample_intake_date, format = "%Y"), .after = sample_intake_date) %>%
  mutate(sample_name = paste(sample_intake_year, individual_id, sep = "_"))

