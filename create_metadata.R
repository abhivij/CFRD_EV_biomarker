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

length(unique(cph_info$individual_id))
#93

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
}

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
  select(-c(10:12, 14, 17:21))
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
colnames(mastersheet_info.sch)[4] <- "sample_date"
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
  select(c(1:4, 21))

mastersheet_info.sch2 <- mastersheet_info.sch2 %>%
  mutate("age" = NA, .after = "modulator") %>%
  mutate("sex" = NA, .after = "rna_extraction") %>%
  mutate("FEV1" = NA, .after = "sex")
colnames(mastersheet_info.sch)[1] <- "condition"
colnames(mastersheet_info.sch2)[1] <- "condition"

mastersheet_info.sch <- rbind(mastersheet_info.sch, mastersheet_info.sch2)
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
  select(c(1,2,3,7,10,11))
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
  relocate(age, .before = sex)

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

meta_data <- rbind(cph_info %>%
                     mutate("cohort" = "CPH", .after = year) %>%
                     mutate("country" = "DK", .after = cohort), 
                   rpa_info %>%
                     mutate("cohort" = "RPA_NSW", .after = year) %>%
                     mutate("country" = "AU", .after = cohort),
                   sch_info %>%
                     mutate("cohort" = "SCH_NSW", .after = year) %>%
                     mutate("country" = "AU", .after = cohort), 
                   healthy_info %>%
                     mutate("cohort" = "UNSW", .after = year) %>%
                     mutate("country" = "AU", .after = cohort))

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
                          TRUE ~ FEV1))
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