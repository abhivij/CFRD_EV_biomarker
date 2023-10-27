library(tidyverse)
library(readxl)
library(ROCR)

library(ggvenn)
library(ComplexHeatmap)
library(viridis)
library(RColorBrewer)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/CysticFibrosisGroup/ExoCF/CFRD_EV_biomarker/"
setwd(base_dir)

#create heatmap of average prediction probabilities

# best_features <- read.csv("Data/selected_features/best_features_with_add_col.csv") %>%
#   filter(is_best == 1, 
#          grepl("_combined_", dataset_id),
#          grepl('_common_combat_', dataset_id)) %>%
#   mutate(comparison = gsub("GBM_combined_transcriptomic_common_combat_|GBM_combined_proteomic_common_combat_",
#                            "", dataset_id)) %>%
#   mutate(dataset_id = paste(dataset_id, description, 
#                             min_iter_feature_presence, comparison,
#                             sep = "_")) %>%
#   dplyr::select(c(dataset_id, comparison))
# 
# sample_wise_results.proteomics <- read.csv("fem_pipeline_results_combined_pr_common_combat_subset/all_samplewise_result_df.csv") %>%
#   inner_join(best_features %>% filter(grepl("proteomic", dataset_id)), 
#              by = c("DataSetId" = "dataset_id")) %>%
#   filter(FSM == "all") %>%
#   dplyr::select(-c(FSM))
# sample_wise_results.transcriptomics <- read.csv("fem_pipeline_results_combined_tr_common_combat_subset/all_samplewise_result_df.csv") %>%
#   inner_join(best_features %>% filter(grepl("transcriptomic", dataset_id)), 
#              by = c("DataSetId" = "dataset_id")) %>%
#   filter(FSM == "all") %>%
#   dplyr::select(-c(FSM))
# # 
# sample_type <- "test"
# comparison_of_interest <- "POSTOPE_TPVsREC_TP"
# classes <- c("REC_TP", "POSTOPE_TP")
# 
# write.csv(head(sample_wise_results.transcriptomics),
#           "Data/tra_results.subset.csv", 
#           row.names = FALSE)


results_file_path <- "integration_prediction_result/CFRDVsIGT.csv"
comparison_of_interest = "CFRDVsIGT"
classes = c("IGT", "CFRD")
sample_type = "train"

#classes assumed to be given as (neg_class, pos_class)
create_mean_prob_heatmap <- function(results_file_path,
                                     comparison_of_interest,
                                     classes,
                                     sample_type){
  
  results <- read.csv(results_file_path) %>%
    dplyr::filter(Type == sample_type)
  
  results <- results %>%
    group_by(omics_type, sample, TrueLabel, model) %>%
    summarize(mean_prob = mean(Pred_prob)) %>% 
    ungroup()
  
  sample_label <- results %>%
    dplyr::select(c(sample, TrueLabel)) %>%
    distinct()
  
  # results.tra <- sample_wise_results.transcriptomics %>%
  #   filter(Type == sample_type, comparison == comparison_of_interest) %>%
  #   dplyr::select(-c(DataSetId, comparison, Type)) %>%
  #   group_by(Sample, TrueLabel, Model) %>%
  #   summarize(mean_prob = mean(PredProb))
  # results.prot <- sample_wise_results.proteomics %>%
  #   filter(Type == sample_type, comparison == comparison_of_interest) %>%
  #   dplyr::select(-c(DataSetId, comparison, Type)) %>%
  #   group_by(Sample, TrueLabel, Model) %>%
  #   summarize(mean_prob = mean(PredProb))
  # 
  # sample_label.tra <- results.tra %>%
  #   dplyr::select(c(Sample, TrueLabel)) %>%
  #   distinct() %>%
  #   arrange(Sample)
  # sample_label.prot <- results.prot %>%
  #   dplyr::select(c(Sample, TrueLabel)) %>%
  #   distinct() %>%
  #   mutate(Sample = sub("HB0", "HB", Sample)) %>%
  #   arrange(Sample)
  # 
  # setdiff(sample_label.prot$Sample, sample_label.tra$Sample)
  # setdiff(sample_label.tra$Sample, sample_label.prot$Sample)
  # 
  # results.prot <- results.prot %>%
  #   mutate(Sample = sub("HB0", "HB", Sample)) 
  
  results <- results %>%
    mutate(omics_type = ifelse(omics_type == "prot", "proteomics", "transcriptomics")) %>%
    mutate(model_omics_type = paste(model, omics_type, sep = "_"))
  
  #show replicate samples together
  all_replicate_samples <- read.table("data/formatted/replicate_prot_samples.txt", sep = "\t",
                                  header = TRUE)
  all_sample_name_pt_mapping <- read.table("data/formatted/pt_mapping_and_common_sample_name.txt", sep = "\t",
                                       header = TRUE) 
  all_replicate_sample_name <- all_replicate_samples %>%
    inner_join(all_sample_name_pt_mapping, by = c("tra_sample_name", "prot_sample_name")) %>%
    arrange(desc(n), tra_sample_name, prot_sample_name)
  
  data_to_plot <- results %>%
    dplyr::select(c(sample, model_omics_type, mean_prob)) %>%
    pivot_wider(names_from = model_omics_type, values_from = mean_prob)
  
  #START : doing the below to specifically show samples with replicates at the start
  sample_name_not_replicate <- sort(data_to_plot$sample[!data_to_plot$sample %in% all_replicate_sample_name$sample_name])
  #using all_replicate_sample_name$sample_name so that the same order as in all_replicate_sample_name is used
  sample_name_replicate <- all_replicate_sample_name$sample_name[all_replicate_sample_name$sample_name %in% data_to_plot$sample]
  sample_name_ordered <- c(sample_name_replicate, sample_name_not_replicate)
  
  data_to_plot <- data_to_plot %>%
    column_to_rownames("sample")
  data_to_plot <- data_to_plot[sample_name_ordered, ]
  #END : doing the below to specifically show samples with replicates at the start
  
  data_to_plot <- data.matrix(data_to_plot)
  
  model_names <- c("Simple logistic regression",
                   "Radial Kernel SVM",
                   "Sigmoid Kernel SVM",
                   "Random Forest",
                   "L1 Regularized logistic regression",
                   "L2 Regularized logistic regression",
                   "Elastic net logistic regression")
  
  meta_data.row <- data.frame(sample = rownames(data_to_plot)) %>% 
    inner_join(sample_label) %>%
    mutate(TrueLabel = factor(TrueLabel, levels = rev(classes))) %>%
    mutate(replicates = case_when(sample %in% sample_name_replicate ~ "replicates present",
                                  TRUE ~ "no replicates"))
  meta_data.col <- data.frame(model_omics_type = colnames(data_to_plot)) %>%
    separate(model_omics_type, into = c("model", "omics_type"), sep = "_", remove = FALSE) %>%
    mutate(model = factor(model, levels = model_names))
  
  row_col <- list()
  row_col[["TrueLabel"]] <- c("#440154", "#fde725")
  names(row_col[["TrueLabel"]]) <- classes
  
  row_col[["Replicates"]] <- c("darkgreen", NA)
  names(row_col[["Replicates"]]) <- c("replicates present", "no replicates")
  
  column_col <- list("Omics type" = c("proteomics" = "skyblue1",
                                      "transcriptomics" = "indianred1"))
  column_col[["Model"]] <- brewer.pal(n = 7, name = "Paired")
  names(column_col[["Model"]]) <- model_names
  
  ht <- Heatmap(data_to_plot, name = "Mean Prediction probability",
                col = viridis(5),
                rect_gp = gpar(col = "white", lwd = 1),
                cluster_columns = TRUE,
                cluster_rows = FALSE,
                row_title = "Samples",
                row_names_side = "left",
                row_split = meta_data.row$TrueLabel,
                column_split = meta_data.col$omics_type,
                column_title = NULL,
                show_column_names = FALSE,
                row_names_gp = gpar(fontsize = 4.5),
                bottom_annotation = HeatmapAnnotation(
                  "Omics type" = meta_data.col$omics_type,
                  "Model" = meta_data.col$model,
                  annotation_name_side = "left",
                  col = column_col
                )) + 
    HeatmapAnnotation("TrueLabel" = meta_data.row$TrueLabel, 
                      "Replicates" = meta_data.row$replicates,
                      which = "row", 
                      col = row_col)
  plot_dir_path <- "plots_updated/integration/base_model/"
  if(!dir.exists(plot_dir_path)){
    dir.create(plot_dir_path, recursive = TRUE)
  }
  png(paste0(plot_dir_path, comparison_of_interest, 
             sample_type,
             ".png"), units = "cm", width = 20, height = 15, res = 1200)  
  draw(ht, column_title = paste(sample_type, sub("Vs", " Vs ", comparison_of_interest)))    
  dev.off() 
  
}

create_mean_prob_heatmap(results_file_path = "integration_prediction_result/CFRDVsIGT.csv",
                         comparison_of_interest = "CFRDVsIGT",
                         classes = c("IGT", "CFRD"),
                         sample_type = "train")
create_mean_prob_heatmap(results_file_path = "integration_prediction_result/CFRDVsNGT.csv",
                         comparison_of_interest = "CFRDVsNGT",
                         classes = c("NGT", "CFRD"),
                         sample_type = "train")
create_mean_prob_heatmap(results_file_path = "integration_prediction_result/IGTVsNGT.csv",
                         comparison_of_interest = "IGTVsNGT",
                         classes = c("NGT", "IGT"),
                         sample_type = "train")

create_mean_prob_heatmap(results_file_path = "integration_prediction_result/CFRDVsIGT.csv",
                         comparison_of_interest = "CFRDVsIGT",
                         classes = c("IGT", "CFRD"),
                         sample_type = "test")
create_mean_prob_heatmap(results_file_path = "integration_prediction_result/CFRDVsNGT.csv",
                         comparison_of_interest = "CFRDVsNGT",
                         classes = c("NGT", "CFRD"),
                         sample_type = "test")
create_mean_prob_heatmap(results_file_path = "integration_prediction_result/IGTVsNGT.csv",
                         comparison_of_interest = "IGTVsNGT",
                         classes = c("NGT", "IGT"),
                         sample_type = "test")                     
