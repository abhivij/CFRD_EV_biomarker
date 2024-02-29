library(sva)

base_dir <- "~/UNSW/VafaeeLab/CysticFibrosisGroup/ExoCF/CFRD_EV_biomarker/"
setwd(base_dir)

source("utils_diff.R")
source("utils.R")


data <- read.csv("data/formatted/rna_all/umi_counts_filter90.csv", row.names = 1)

# filter out HC and unknown
# quantile norm
# combat

phenotype <- read.table("data/formatted/tra_phenotype_2024Jan.txt", header = TRUE, sep = "\t") %>%
  mutate(modstatus_condition_country = paste(PreModulatorVsPostModulator, condition, country, sep = "_")) %>%
  mutate(modstatus_condition = paste(PreModulatorVsPostModulator, condition, sep = "_"))

summary(factor(phenotype$condition))
# CFRD                 HC                IGT                NGT UNKNOWN TO PREDICT 
# 82                  6                 80                138                 28


# %>%
#   dplyr::rename(c("disease_status" = "condition"))

phenotype <- phenotype %>%
  filter(condition != "UNKNOWN TO PREDICT" & condition != "HC")
summary(factor(phenotype$modstatus_condition_country))
summary(factor(phenotype$modstatus_condition))

data <- data[, phenotype$Sample]

#normalize
data <- edgeR::cpm(data, log = TRUE) 

all.equal(colnames(data), phenotype$Sample)
#TRUE

phenotype <- phenotype %>%
  mutate(country_batch_name = paste(country, batch_name, sep = "_"))
summary(factor(phenotype$country_batch_name))

# AU_initial DK_initial     DK_new 
# 105        133         62 