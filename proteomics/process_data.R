library(tidyverse)
library(magrittr)
library(DEP)
library(missForest)


source('utils.R')

data_main <- read.table(file = "data/proteomics/AU_DK/proteinGroups.txt", sep = "\t", header = TRUE)
data_other <- read.table(file = "data/proteomics/AU_DK/DK_others/proteinGroups.txt", sep = "\t", header = TRUE)

meta_data <- read.csv("data/proteomics/prot_metadata_updated.csv") %>%
  mutate(modulator_status = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ 'postmod',
                                      TRUE ~ 'premod')) %>%
  mutate(mstat_condition_country = paste(modulator_status, condition, country, sep = "_")) %>%
  dplyr::rename(c("disease_status" = "condition"))
summary(factor(meta_data$mstat_condition_country))
#assuming that similar to tra, batches will be for country and not the cohorts within AU
# using this mstat_condition_country column to identify proteins to be filtered


data_main <- data_main[data_main$Reverse != "+" & data_main$Potential.contaminant != "+" & data_main$Unique.peptides >= 2,]
data_main <- data_main[data_main$Gene.names != "",] #removes rows with no Gene names
data_main <- data_main[complete.cases(data_main),] #removes NAs at the end 
data_main_unique <- make_unique(data_main, "Gene.names", "Protein.IDs", delim = ";") 
# data_main_unique <- data_main_unique %>%
#   select(-contains("glufib.")) %>% #remove glufib
#   select(-contains("QC"))#remove QC samples

meta_data_main <- meta_data %>%
  filter(mq_batch == "main")
meta_data_main <- add_replicate_column(meta_data_main, "mstat_condition_country")

exp_design <- meta_data_main %>%
  select(label,condition,replicate) %>%
  mutate(label=as.character(label))
exp_design$label <- paste("LFQ.intensity.",exp_design$label,sep="")

columns <- grep("LFQ.intensity.", colnames(data_main_unique))

se_main <- make_se(data_main_unique, columns, exp_design)

# plot_coverage(se_main)

dim(assays(se_main)[[1]])
# [1] 1539  252

min(summary(factor(meta_data_main$condition)) / dim(meta_data_main)[1])
#0.01190476

max(summary(factor(meta_data_main$condition)) / dim(meta_data_main)[1])
#0.1190476

summary(summary(factor(meta_data_main$condition)) / dim(meta_data_main)[1])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01190 0.04762 0.05159 0.06667 0.08929 0.11905

#proteins present in atleast number of samples equal to half of the samples of the smallest group
# this will also include proteins present in samples from one group, and few other samples from another group
# so we loose less proteins
filt_main <- filter_proteins(se_main, "fraction", min = 0.01190/2)
dim(assays(filt_main)[[1]])
# [1] 961 252




data_other <- data_other[data_other$Reverse != "+" & data_other$Potential.contaminant != "+" & data_other$Unique.peptides >= 2,]
data_other <- data_other[data_other$Gene.names != "",] #removes rows with no Gene names
data_other <- data_other[complete.cases(data_other),] #removes NAs at the end 
data_other_unique <- make_unique(data_other, "Gene.names", "Protein.IDs", delim = ";") 
# data_other_unique <- data_other_unique %>%
#   select(-contains("glufib.")) %>% #remove glufib
#   select(-contains("QC"))#remove QC samples

meta_data_other <- meta_data %>%
  filter(mq_batch == "other")
meta_data_other <- add_replicate_column(meta_data_other, "mstat_condition_country")

exp_design <- meta_data_other %>%
  select(label,condition,replicate) %>%
  mutate(label=as.character(label))
exp_design$label <- paste("LFQ.intensity.",exp_design$label,sep="")

columns <- grep("LFQ.intensity.", colnames(data_other_unique))

se_other <- make_se(data_other_unique, columns, exp_design)

plot_coverage(se_other)

dim(assays(se_other)[[1]])
# [1] 666  18

min(summary(factor(meta_data_other$condition)) / dim(meta_data_other)[1])
#0.05555556

max(summary(factor(meta_data_other$condition)) / dim(meta_data_other)[1])
#0.3888889

summary(summary(factor(meta_data_other$condition)) / dim(meta_data_other)[1])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.05556 0.06944 0.13889 0.16667 0.20833 0.38889

filt_other <- filter_proteins(se_other, "fraction", min = 0.05556/2)
dim(assays(filt_other)[[1]])
# [1] 522  18


#######create venn diagram of filtered proteins

filt_main_data <- get_df_long(filt_main) %>%
  select(name, intensity) %>%
  filter(!is.na(intensity)) %>%
  distinct(name)
filt_other_data <- get_df_long(filt_other) %>%
  select(name, intensity) %>%
  filter(!is.na(intensity)) %>%
  distinct(name)  

ggvenn::ggvenn(list("main" = filt_main_data$name,
                    "other" = filt_other_data$name),
               stroke_size = 0.1,
               set_name_size = 4,
               text_size = 3,
               fill_color = c("#F8766D", "#00BFC4")) +
  ggtitle("Proteins in the two MaxQuant Batches after filter") +
  theme(plot.title = element_text(vjust = -10, hjust = 0.5))
file_name <- paste0("plots_updated/proteomics/protein_overlap_mq_filtered.png")
ggsave(file_name)

proteins_common <- intersect(filt_main_data$name, filt_other_data$name)
proteins_specific_in_main <- setdiff(filt_main_data$name, filt_other_data$name)
proteins_specific_in_other <- setdiff(filt_other_data$name, filt_main_data$name)

#############

filt_main_data <- get_df_long(filt_main) %>%
  dplyr::select(c(label, name, intensity)) %>%
  pivot_wider(names_from = name, values_from = intensity) 
filt_main_data[, proteins_specific_in_other] <- NA

filt_other_data <- get_df_long(filt_other) %>%
  dplyr::select(c(label, name, intensity)) %>%
  pivot_wider(names_from = name, values_from = intensity)
filt_other_data[, proteins_specific_in_main] <- NA

filt_main_data <- filt_main_data %>%
  dplyr::select(label, all_of(c(proteins_common, proteins_specific_in_main, proteins_specific_in_other)))
filt_other_data <- filt_other_data %>%
  dplyr::select(label, all_of(c(proteins_common, proteins_specific_in_main, proteins_specific_in_other)))
all.equal(colnames(filt_main_data), colnames(filt_other_data))

data.combined <- rbind(filt_main_data, filt_other_data) %>%
  mutate(label = sub("LFQ.intensity.", "lfq_", label)) %>%
  column_to_rownames("label")
data.main <- filt_main_data %>%
  mutate(label = sub("LFQ.intensity.", "lfq_", label)) %>%
  column_to_rownames("label") 
data.other <- filt_other_data %>%
  mutate(label = sub("LFQ.intensity.", "lfq_", label)) %>%
  column_to_rownames("label") 

sum(!is.na(data.combined)) / (nrow(data.combined) * ncol(data.combined))
#0.2346025

sum(!is.na(data.main)) / (nrow(data.main) * ncol(data.main))
#0.2352363

sum(!is.na(data.other)) / (nrow(data.other) * ncol(data.other))
#0.2257292

#to be noted : percentage of non-NAs are really low

# all data in sample x proteins

data.imputed.combined <- t(missForest(data.combined, verbose = TRUE)$ximp)
data.imputed.main <- t(missForest(data.main, verbose = TRUE)$ximp)

data.imputed.combined <- as.data.frame(data.imputed.combined)
data.imputed.main <- as.data.frame(data.imputed.main)

sum(is.na(data.imputed.combined))
#0
sum(is.na(data.imputed.main))
#0

data.zeroimputed.combined <- data.combined
data.zeroimputed.combined[is.na(data.zeroimputed.combined)] <- 0
data.zeroimputed.combined <- as.data.frame(t(data.zeroimputed.combined))

data.zeroimputed.main <- data.main[, c(proteins_common, proteins_specific_in_main)]
data.zeroimputed.main[is.na(data.zeroimputed.main)] <- 0
data.zeroimputed.main <- as.data.frame(t(data.zeroimputed.main))

write.csv(data.imputed.combined, "data/proteomics/imputed_combined.csv")
data.imputed.combined2 <- read.csv("data/proteomics/imputed_combined.csv", row.names = 1)
all.equal(data.imputed.combined, data.imputed.combined2)

write.csv(data.zeroimputed.combined, "data/proteomics/zeroimputed_combined.csv")
data.zeroimputed.combined2 <- read.csv("data/proteomics/zeroimputed_combined.csv", row.names = 1)
all.equal(data.zeroimputed.combined, data.zeroimputed.combined2)

write.csv(data.imputed.main, "data/proteomics/imputed_main.csv")
data.imputed.main2 <- read.csv("data/proteomics/imputed_main.csv", row.names = 1)
all.equal(data.imputed.main, data.imputed.main2)

write.csv(data.zeroimputed.main, "data/proteomics/zeroimputed_main.csv")
data.zeroimputed.main2 <- read.csv("data/proteomics/zeroimputed_main.csv", row.names = 1)
all.equal(data.zeroimputed.main, data.zeroimputed.main2)




####################
################################################

#processing the MQ result with all samples i.e. 333 samples
#         this includes the DK samples provided in 2023

data <- read.table(file = "data/proteomics/all/proteinGroups.txt", sep = "\t", header = TRUE)

meta_data <- read.csv("data/proteomics/prot_metadata_all_2023Aug.csv") %>%
  mutate(modulator_status = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ 'postmod',
                                      TRUE ~ 'premod')) %>%
  mutate(mstat_condition_country = paste(modulator_status, condition, country, sep = "_")) %>%
  dplyr::rename(c("disease_status" = "condition"))
summary(factor(meta_data$mstat_condition_country))
#assuming that similar to tra, batches will be for country and not the cohorts within AU
# using this mstat_condition_country column to identify proteins to be filtered


data <- data[data$Reverse != "+" & data$Potential.contaminant != "+" & data$Unique.peptides >= 2,]
data <- data[data$Gene.names != "",] #removes rows with no Gene names
data <- data[complete.cases(data),] #removes NAs at the end 
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";") 
# data_unique <- data_unique %>%
#   select(-contains("glufib.")) %>% #remove glufib
#   select(-contains("QC"))#remove QC samples


######create protein names file

protein_names <- data_unique %>% dplyr::select("Protein.IDs", "Protein.names", "Gene.names", "name", "ID")
colnames(protein_names) <- c("protein_id", "protein_name", "gene_name", "name", "ID")
write.csv(protein_names, file="data/proteomics/protein_names.csv", row.names = F)

######################


meta_data <- add_replicate_column(meta_data, "mstat_condition_country")

exp_design <- meta_data %>%
  select(label,condition,replicate) %>%
  mutate(label=as.character(label))
exp_design$label <- paste("LFQ.intensity.",exp_design$label,sep="")

columns <- grep("LFQ.intensity.", colnames(data_unique))

se <- make_se(data_unique, columns, exp_design)

#plot_coverage(se)

dim(assays(se)[[1]])

summary(summary(factor(meta_data$condition)) / dim(meta_data)[1])
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.009009 0.037538 0.066066 0.066667 0.084084 0.186186

filt <- filter_proteins(se, "fraction", min = 0.009009/2)
dim(assays(filt)[[1]])

filt_data <- get_df_long(filt) %>%
  dplyr::select(c(label, name, intensity)) %>%
  pivot_wider(names_from = name, values_from = intensity) %>%
  mutate(label = sub("LFQ.intensity.", "lfq_", label)) %>%
  column_to_rownames("label") 


sum(!is.na(filt_data)) / (nrow(filt_data) * ncol(filt_data))
#0.2377119

#to be noted : percentage of non-NAs are really low

# data in sample x proteins

data.imputed.mf <- t(missForest(filt_data, verbose = TRUE)$ximp)
# missForest iteration 1 in progress...done!
#   estimated error(s): 0.3976678 
# difference(s): 0.00098186 
# time: 259.585 seconds
# 
# missForest iteration 2 in progress...done!
#   estimated error(s): 0.3959015 
# difference(s): 6.932778e-05 
# time: 297.052 seconds
# 
# missForest iteration 3 in progress...done!
#   estimated error(s): 0.3953103 
# difference(s): 4.393199e-05 
# time: 292.296 seconds
# 
# missForest iteration 4 in progress...done!
#   estimated error(s): 0.3965392 
# difference(s): 4.068473e-05 
# time: 304.098 seconds
# 
# missForest iteration 5 in progress...done!
#   estimated error(s): 0.395222 
# difference(s): 3.894847e-05 
# time: 294.523 seconds
# 
# missForest iteration 6 in progress...done!
#   estimated error(s): 0.3968719 
# difference(s): 4.013596e-05 
# time: 287.009 seconds


data.imputed.mf <- as.data.frame(data.imputed.mf)

sum(is.na(data.imputed.mf))
#0

data.imputed.zero <- filt_data
sum(is.na(data.imputed.zero))
#393455
data.imputed.zero[is.na(data.imputed.zero)] <- 0
data.imputed.zero <- as.data.frame(t(data.imputed.zero))
sum(is.na(data.imputed.zero))

write.csv(data.imputed.mf, "data/proteomics/data_333samples_imputed_mf.csv")

write.csv(data.imputed.zero, "data/proteomics/data_333samples_imputed_zero.csv")

#writing once again to get column names in format converted to by read.csv
data.imputed.mf <- read.csv("data/proteomics/data_333samples_imputed_mf.csv", row.names = 1)
data.imputed.zero <- read.csv("data/proteomics/data_333samples_imputed_zero.csv", row.names = 1)

write.csv(data.imputed.mf, "data/proteomics/data_333samples_imputed_mf.csv")
write.csv(data.imputed.zero, "data/proteomics/data_333samples_imputed_zero.csv")
##########################################
#ignoring the 'other' batch samples

data <- read.table(file = "data/proteomics/all/proteinGroups.txt", sep = "\t", header = TRUE)

meta_data_no_other <- read.csv("data/proteomics/prot_metadata_all_2023Aug.csv") %>%
  mutate(modulator_status = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ 'postmod',
                                      TRUE ~ 'premod')) %>%
  mutate(mstat_condition_country = paste(modulator_status, condition, country, sep = "_")) %>%
  dplyr::rename(c("disease_status" = "condition")) %>%
  filter(batch_name != "other")
summary(factor(meta_data_no_other$mstat_condition_country))


data <- data[data$Reverse != "+" & data$Potential.contaminant != "+" & data$Unique.peptides >= 2,]
data <- data[data$Gene.names != "",] #removes rows with no Gene names
data <- data[complete.cases(data),] #removes NAs at the end 
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";") 
# data_unique <- data_unique %>%
#   select(-contains("glufib.")) %>% #remove glufib
#   select(-contains("QC"))#remove QC samples


meta_data_no_other <- add_replicate_column(meta_data_no_other, "mstat_condition_country")

exp_design_no_other <- meta_data_no_other %>%
  select(label,condition,replicate) %>%
  mutate(label=as.character(label))
exp_design_no_other$label <- paste("LFQ.intensity.",exp_design_no_other$label,sep="")

columns <- grep("LFQ.intensity.", colnames(data_unique))

se_no_other <- make_se(data_unique, columns, exp_design_no_other)

#plot_coverage(se_no_other)

dim(assays(se_no_other)[[1]])
# [1] 1658  315

summary(summary(factor(meta_data_no_other$condition)) / dim(meta_data_no_other)[1])
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.009524 0.038095 0.060317 0.066667 0.084127 0.184127 

filt_no_other <- filter_proteins(se_no_other, "fraction", min = 0.009524/2)
dim(assays(filt_no_other)[[1]])
# [1] 1546  315

filt_data_no_other <- get_df_long(filt_no_other) %>%
  dplyr::select(c(label, name, intensity)) %>%
  pivot_wider(names_from = name, values_from = intensity) %>%
  mutate(label = sub("LFQ.intensity.", "lfq_", label)) %>%
  column_to_rownames("label") 


sum(!is.na(filt_data_no_other)) / (nrow(filt_data_no_other) * ncol(filt_data_no_other))
#0.2388673

#to be noted : percentage of non-NAs still low 
# 0.2377119 became 0.2388673


# data in sample x proteins
data_no_other.imputed.mf <- t(missForest(filt_data_no_other, verbose = TRUE)$ximp)
# missForest iteration 1 in progress...done!
#   estimated error(s): 0.396851 
# difference(s): 0.0009781329 
# time: 229.65 seconds
# 
# missForest iteration 2 in progress...done!
#   estimated error(s): 0.3956837 
# difference(s): 6.474067e-05 
# time: 268.447 seconds
# 
# missForest iteration 3 in progress...done!
#   estimated error(s): 0.3959579 
# difference(s): 4.357137e-05 
# time: 272.157 seconds
# 
# missForest iteration 4 in progress...done!
#   estimated error(s): 0.3953812 
# difference(s): 4.03907e-05 
# time: 288.288 seconds
# 
# missForest iteration 5 in progress...done!
#   estimated error(s): 0.3946 
# difference(s): 3.813403e-05 
# time: 288.526 seconds
# 
# missForest iteration 6 in progress...done!
#   estimated error(s): 0.3959512 
# difference(s): 3.794324e-05 
# time: 288.897 seconds
# 
# missForest iteration 7 in progress...done!
#   estimated error(s): 0.3956366 
# difference(s): 3.889752e-05 
# time: 289.967 seconds
data_no_other.imputed.mf <- as.data.frame(data_no_other.imputed.mf)

sum(is.na(data_no_other.imputed.mf))
#0

data_no_other.imputed.zero <- filt_data_no_other
sum(is.na(data_no_other.imputed.zero))
#370664
data_no_other.imputed.zero[is.na(data_no_other.imputed.zero)] <- 0
data_no_other.imputed.zero <- as.data.frame(t(data_no_other.imputed.zero))
sum(is.na(data_no_other.imputed.zero))

write.csv(data_no_other.imputed.mf, "data/proteomics/data_no_other_315samples_imputed_mf.csv")

write.csv(data_no_other.imputed.zero, "data/proteomics/data_no_other_315samples_imputed_zero.csv")

#writing once again to get column names in format converted to by read.csv
data_no_other.imputed.mf <- read.csv("data/proteomics/data_no_other_315samples_imputed_mf.csv", row.names = 1)
data_no_other.imputed.zero <- read.csv("data/proteomics/data_no_other_315samples_imputed_zero.csv", row.names = 1)

write.csv(data_no_other.imputed.mf, "data/proteomics/data_no_other_315samples_imputed_mf.csv")
write.csv(data_no_other.imputed.zero, "data/proteomics/data_no_other_315samples_imputed_zero.csv")
