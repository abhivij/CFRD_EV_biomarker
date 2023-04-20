library(xlsx)
library(ggvenn)
library(ComplexHeatmap)
library(readxl)
library(tidyverse)

results <- read_excel("prediction_pipeline/results_sch.xlsx")
avatar_master_sheet_info <- xlsx::read.xlsx("data/Avatar Master Database_1MAR2022.xlsx", 
                                            password = "CF", sheetName = "Children_CF_Master") %>%
  dplyr::select(c(11,15))
colnames(avatar_master_sheet_info) <- avatar_master_sheet_info[1,]
avatar_master_sheet_info <- avatar_master_sheet_info[-c(1),]

rowSums(is.na(avatar_master_sheet_info)) != ncol(avatar_master_sheet_info)

data_of_interest <- avatar_master_sheet_info[rowSums(is.na(avatar_master_sheet_info)) != ncol(avatar_master_sheet_info), ]

results_with_MRN <- results %>% inner_join(data_of_interest, by = c("Patient ID" = "Patient ID "))
results_with_MRN <- as.data.frame(results_with_MRN)

write.xlsx(results_with_MRN,
           "prediction_pipeline/results_sch_with_MRN.xlsx",
           col.names = TRUE, row.names = FALSE)
