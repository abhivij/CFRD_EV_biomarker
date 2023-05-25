library(tidyverse)
library(DEP)
library(magrittr)

library(xlsx)


base_dir <- "~/UNSW/VafaeeLab/CysticFibrosisGroup/ExoCF/CFRD_EV_biomarker/"
setwd(base_dir)

data <- read.table(file = "data/proteomics/AU_DK/proteinGroups.txt", sep = "\t", header = TRUE)
summary <- read.table(file = "data/proteomics/AU_DK/summary.txt", sep = "\t", header = TRUE) %>%
  filter(Raw.file != "Total") %>%
  arrange(Raw.file) %>%
  mutate(Experiment = gsub("-", ".", Experiment, fixed = TRUE)) %>%
  mutate(Experiment = gsub("_", ".", Experiment, fixed = TRUE))

sinfo <- read.table(file = "data/proteomics/AU_DK/proteomics_sinfo.txt", sep = "\t", header = TRUE) %>%
  arrange(rawfile)
colnames(summary)[1:2] <- colnames(sinfo)[1:2]

all.equal(summary[, c(1:2)], sinfo[, c(1:2)])


data <- data[data$Reverse != "+" & data$Potential.contaminant != "+" & data$Unique.peptides >= 2,]
data <- data[data$Gene.names != "",] #removes rows with no Gene names
data <- data[complete.cases(data),] #removes NAs at the end 
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")

lfq_intensity <- data_unique %>%
  set_rownames(.$Gene.names) %>%
  select(contains("LFQ.intensity.")) %>%
  select(-contains("glufib.")) %>% #remove glufib
  select(-contains("QC")) %>% #remove QC samples
  log2()

lfq_intensity[lfq_intensity == "-Inf"] <- 0

colnames(lfq_intensity) <- sub("LFQ.intensity.","", colnames(lfq_intensity)) #subsitutes colnames with LFQ.intensity removed
mm <- match(names(lfq_intensity), sinfo$label) #match names with label
names(lfq_intensity)[!is.na(mm)] <- as.character(sinfo$sname[na.omit(mm)])



# #check number of modulator AU samples
# 
# meta_data_sub <- meta_data %>%
#   filter(country == "AU", pre_post_modulator == 1)
# 
# length(unique(meta_data_sub$individual_id))
# write.csv(meta_data_sub, "data/au_post_modulator.csv", row.names = FALSE)


#compare data from DK all
data_dkall <- read.table(file = "data/proteomics/2_UK/all_data/proteinGroups.txt", 
                         sep = "\t", header = TRUE)
summary_dkall <- read.table(file = "data/proteomics/2_UK/all_data/summary.txt", 
                            sep = "\t", header = TRUE) %>%
  filter(Raw.file != "Total") %>%
  arrange(Raw.file) %>%
  mutate(Experiment = gsub("-", ".", Experiment, fixed = TRUE)) %>%
  mutate(Experiment = gsub("_", ".", Experiment, fixed = TRUE))
colnames(summary_dkall)[1:2] <- colnames(sinfo)[1:2]
all.equal(summary_dkall[, c(1:2)], sinfo[, c(1:2)])


#compare sinfo against meta_data used for transcriptomics
meta_data <- read.csv("data/formatted/meta_data.csv")

sinfo_sub <- sinfo %>%
  dplyr::select(c(label, individualid, year, agegroup, cohort, condition)) %>%
  filter(!individualid %in% c("glufib", "QC"))

meta_data_sub <- meta_data %>%
  dplyr::select(c(sample_long_name, individual_id, patient_recruitment_year, age_group, country, condition))

colnames(sinfo_sub)[1] <- "sample_long_name"
colnames(meta_data_sub) <- colnames(sinfo_sub)

meta_data_sub <- meta_data_sub %>%
  arrange(cohort, individualid, year)
sinfo_sub <- sinfo_sub %>%
  arrange(cohort, individualid, year) %>%
  mutate(agegroup = case_when(agegroup == "Adult" ~ "adult",
                              agegroup == "Child" ~ "child")) %>%
  mutate(agegroup = ifelse(individualid %in% c("067ES", "068CG", "069WC", 
                                               "070AM", "071LO", "072CS", 
                                               "073MA", "074CS", "075RO"), "adult", agegroup)) %>%
  mutate(condition = ifelse(condition == "INDET", "IND", condition)) %>%
  mutate(individualid = sub("CPH0{1,2}", "CPH", individualid, fixed = FALSE)) 

write.csv(meta_data_sub, "data/proteomics/metadata_sub.csv", row.names = FALSE)
write.csv(sinfo_sub, "data/proteomics/sinfo_sub.csv", row.names = FALSE)

meta_data_sub.dk <- meta_data_sub %>%
  filter(cohort == "DK") %>%
  arrange(individualid) %>%
  dplyr::select(-c(cohort))
sinfo_sub.dk <- sinfo_sub %>%
  filter(cohort == "DK") %>%
  arrange(individualid) %>%
  dplyr::select(-c(cohort))

combined_metadata_sub.dk <- meta_data_sub.dk %>%
  full_join(sinfo_sub.dk, by = c("individualid", "year"), suffix = c("_t", "_p")) %>%
  mutate(ind_id = as.numeric(sub("CPH", "", individualid))) %>%
  arrange(ind_id) %>%
  dplyr::select(-c(ind_id))

#missing in tra
missing_t.dk <- combined_metadata_sub.dk %>%
  filter(is.na(sample_long_name_t)) %>%
  dplyr::select(c(individualid, year, agegroup_p, condition_p, sample_long_name_p))

#missing in prot
missing_p.dk <- combined_metadata_sub.dk %>%
  filter(is.na(sample_long_name_p)) %>%
  dplyr::select(c(individualid, year, agegroup_t, condition_t, sample_long_name_t))

write.xlsx(combined_metadata_sub.dk,
           "data/proteomics/meta_data_combined.xlsx",
           sheetName = "combined DK",
           col.names = TRUE, row.names = FALSE, append = FALSE)

combined_metadata_sub.dk.non_missing <- combined_metadata_sub.dk %>%
  filter(!is.na(sample_long_name_t) & !is.na(sample_long_name_p)) %>%
  filter(agegroup_t == agegroup_p & condition_t == condition_p) %>%
  dplyr::select(-c(agegroup_p, condition_p)) %>%
  dplyr::rename(c("agegroup" = "agegroup_t", "condition" = "condition_t"))

write.xlsx(missing_t.dk,
           "data/proteomics/meta_data_combined.xlsx",
           sheetName = "DK missing in tra",
           col.names = TRUE, row.names = FALSE, append = TRUE)
write.xlsx(missing_p.dk,
           "data/proteomics/meta_data_combined.xlsx",
           sheetName = "DK missing in prot",
           col.names = TRUE, row.names = FALSE, append = TRUE)
write.xlsx(combined_metadata_sub.dk.non_missing,
           "data/proteomics/meta_data_combined.xlsx",
           sheetName = "combined DK non missing",
           col.names = TRUE, row.names = FALSE, append = TRUE)

