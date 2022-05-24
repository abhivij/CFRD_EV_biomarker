library(multiMiR)
library(tidyverse)

#obtain mapping with IGFBP7

validated_mirnas <- get_multimir(
  org     = "hsa",
  target  = "IGFBP7",
  table   = "validated",
  summary = TRUE
)
validated_mirnas@data
validated_mirnas@summary

validated_mirnas.data <- validated_mirnas@data
validated_mirnas.summary <- validated_mirnas@summary

length(unique(validated_mirnas.data$mature_mirna_id))
#38 same as summary

# predicted_mirnas <- get_multimir(org     = "hsa",
#                          target  = "IGFBP7",
#                          table   = "validated",
#                          summary = TRUE,
#                          predicted.cutoff      = 35,
#                          predicted.cutoff.type = "p",
#                          predicted.site        = "all")

predicted_mirnas <- get_multimir(
  org     = "hsa",
  target  = "IGFBP7",
  table   = "predicted",
  summary = TRUE
)
predicted_mirnas.data <- predicted_mirnas@data
predicted_mirnas.summary <- predicted_mirnas@summary

length(unique(predicted_mirnas.data$mature_mirna_id))
#58 same as summary



umi_counts <- read.csv("data/formatted/umi_counts.csv", row.names = 1)
colnames(umi_counts) <- gsub("^X", "", colnames(umi_counts))

umi_counts <- data.frame(t(umi_counts))


validated_mirnas.data <- validated_mirnas@data %>%
  arrange(mature_mirna_id)
validated_mirnas.data <- validated_mirnas.data %>%
  select(mature_mirna_id, target_symbol, database, support_type) %>%
  mutate(mature_mirna_id = gsub("-", ".", mature_mirna_id, fixed = TRUE)) %>%
  arrange(mature_mirna_id)
#some mirnas have both positive and negative support type
mirnas_of_interest <- unique(validated_mirnas.data$mature_mirna_id)

length(colnames(umi_counts) %in% mirnas_of_interest)
sum(colnames(umi_counts) %in% mirnas_of_interest)

umi_counts.moi <- umi_counts[, colnames(umi_counts) %in% mirnas_of_interest]
sum(colnames(umi_counts.moi) %in% mirnas_of_interest)
