library(tidyverse)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/CysticFibrosisGroup/ExoCF/CFRD_EV_biomarker/"
setwd(base_dir)

source('FEMPipeline/custom_plots_create.R')

best_features <- read.csv("data/selected_features/best_features_with_is_best.csv")   %>%
  filter(is_best == 1, grepl("CF_EV_adult_filtered_seurat3_norm_find_var_none_", dataset_id)) %>%
  separate(dataset_id, into = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, "comparison"), remove = FALSE) %>%
  mutate(dataset_id = paste(dataset_id, description, 
                            min_iter_feature_presence, comparison,
                            sep = "_")) %>%
  mutate(cm = c("Random Forest", "L2 Regularized logistic regression", "Sigmoid Kernel SVM")) %>%
  dplyr::select(c(dataset_id, cm, comparison))

sample_wise_results <- read.csv("fem_pipeline_results_adult_filtered_then_seurat3_norm_find_var_none_subset/all_samplewise_result_df.csv") %>%
  inner_join(best_features, by = c("DataSetId" = "dataset_id", "Model" = "cm")) %>%
  mutate(Model = sub("L2 Regularized logistic regression", "L2 Regularized\nlogistic regression", Model)) %>%
  mutate(Model = sub("Sigmoid Kernel SVM", "Sigmoid\nKernel SVM", Model))


sample_wise_results.train <- sample_wise_results %>%
  filter(Type == "train")
sample_wise_results.test <- sample_wise_results %>%
  filter(Type == "test")

current_comparison <- "CFRDVsNGT"
classes <- c("NGT", "CFRD")

sample_type <- "test"
dir_path <- "prediction_pipeline/plots/ROC/"


plot_roc(sample_wise_results, sample_type = "test", 
         current_comparison = "CFRDVsIGT",
         classes = c("IGT", "CFRD"))
plot_roc(sample_wise_results, sample_type = "test", 
         current_comparison = "CFRDVsNGT",
         classes = c("NGT", "CFRD"))
plot_roc(sample_wise_results, sample_type = "test", 
         current_comparison = "IGTVsNGT",
         classes = c("NGT", "IGT"))
plot_roc(sample_wise_results, sample_type = "train", 
         current_comparison = "CFRDVsIGT",
         classes = c("IGT", "CFRD"))
plot_roc(sample_wise_results, sample_type = "train", 
         current_comparison = "CFRDVsNGT",
         classes = c("NGT", "CFRD"))
plot_roc(sample_wise_results, sample_type = "train", 
         current_comparison = "IGTVsNGT",
         classes = c("NGT", "IGT"))






# plot(perf,
#      lty = 3, lwd = 0.5,
#      colorize = TRUE,
#      main='ROC curves from 6 repeats of 5-fold cross-validation')
# 
# plot(perf,
#      avg = 'horizontal',
#      spread.estimate='boxplot',
#      lwd=2, add = TRUE)
# 
# plot(perf,
#      avg = 'vertical',
#      spread.estimate='boxplot',
#      lwd=2, add = TRUE)
# 
# plot(perf,
#      colorize = TRUE,
#      avg = 'threshold',
#      lwd = 2, add = TRUE)
# 
# 
# for(i in c(1:30)){
#   i <- 1
#   results_subset <- sample_wise_results.test %>%
#     filter(Iter == i, comparison == current_comparison)
#   
#   pr <- ROCR::prediction(results_subset$PredProb, 
#                          results_subset$TrueLabel, label.ordering = classes)
#   #compute ROC curve, and AUC  
#   prf <- ROCR::performance(pr, measure = "tpr", x.measure = "fpr")
#   if(i == 1){
#     plot(prf, col = i)  
#   } else{
#     # plot(prf, col = i, add = TRUE)
#   }
#   
#   print(ROCR::performance(pr, measure = "auc")@y.values[[1]])
# }
# 
# 
# library(ROCR)
# data(ROCR.simple)
# pred <- prediction(ROCR.simple$predictions,ROCR.simple$labels)
# pred
# prf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
# plot(prf, col = i)  
# 
# 
# data(ROCR.xval)
# 
# predictions <- ROCR.xval$predictions
# labels <- ROCR.xval$labels
# length(predictions)
# 
# pred <- prediction(predictions, labels)
# perf <- performance(pred,'tpr','fpr')
# 
# plot(perf,
#      colorize=TRUE,
#      lwd=2,
#      main='ROC curves from 10-fold cross-validation')
# 
# 
# 
# df <- data.frame(x = 1:2)
# df$y <- matrix(1:4, ncol = 2)


data(aSAH)
roc.s100b <- roc(aSAH$outcome, aSAH$s100b)
roc.wfns <- roc(aSAH$outcome, aSAH$wfns)
roc.ndka <- roc(aSAH$outcome, aSAH$wfns)

# Simple example:
plot(roc.s100b)




create_mean_prob_heatmap(sample_wise_results,
                         comparison_of_interest = "CFRDVsIGT",
                         classes = c("IGT", "CFRD"),
                         sample_type = "test")
create_mean_prob_heatmap(sample_wise_results,
                         comparison_of_interest = "CFRDVsIGT",
                         classes = c("IGT", "CFRD"),
                         sample_type = "train")




create_all_iter_pred_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsIGT",
                             classes = c("IGT", "CFRD"),
                             sample_type = "test", phenotype,
                             plot_dir_path = "plots/custom_heatmap/per_iter_pred/")
create_all_iter_pred_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsNGT",
                             classes = c("NGT", "CFRD"),
                             sample_type = "test", phenotype,
                             plot_dir_path = "plots/custom_heatmap/per_iter_pred/")
create_all_iter_pred_heatmap(sample_wise_results,
                             comparison_of_interest = "IGTVsNGT",
                             classes = c("NGT", "IGT"),
                             sample_type = "test", phenotype,
                             plot_dir_path = "plots/custom_heatmap/per_iter_pred/")

create_all_iter_prob_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsIGT",
                             classes = c("IGT", "CFRD"),
                             sample_type = "test", phenotype,
                             plot_dir_path = "plots/custom_heatmap/per_iter_prob/")
create_all_iter_prob_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsNGT",
                             classes = c("NGT", "CFRD"),
                             sample_type = "test", phenotype,
                             plot_dir_path = "plots/custom_heatmap/per_iter_prob/")
create_all_iter_prob_heatmap(sample_wise_results,
                             comparison_of_interest = "IGTVsNGT",
                             classes = c("NGT", "IGT"),
                             sample_type = "test", phenotype,
                             plot_dir_path = "plots/custom_heatmap/per_iter_prob/")


create_all_iter_pred_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsIGT",
                             classes = c("IGT", "CFRD"),
                             sample_type = "train", phenotype,
                             plot_dir_path = "plots/custom_heatmap/per_iter_pred/")
create_all_iter_pred_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsNGT",
                             classes = c("NGT", "CFRD"),
                             sample_type = "train", phenotype,
                             plot_dir_path = "plots/custom_heatmap/per_iter_pred/")
create_all_iter_pred_heatmap(sample_wise_results,
                             comparison_of_interest = "IGTVsNGT",
                             classes = c("NGT", "IGT"),
                             sample_type = "train", phenotype,
                             plot_dir_path = "plots/custom_heatmap/per_iter_pred/")

create_all_iter_prob_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsIGT",
                             classes = c("IGT", "CFRD"),
                             sample_type = "train", phenotype,
                             plot_dir_path = "plots/custom_heatmap/per_iter_prob/")
create_all_iter_prob_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsNGT",
                             classes = c("NGT", "CFRD"),
                             sample_type = "train", phenotype,
                             plot_dir_path = "plots/custom_heatmap/per_iter_prob/")
create_all_iter_prob_heatmap(sample_wise_results,
                             comparison_of_interest = "IGTVsNGT",
                             classes = c("NGT", "IGT"),
                             sample_type = "train", phenotype,
                             plot_dir_path = "plots/custom_heatmap/per_iter_prob/")


########################

#heatmaps for transcriptomics - created using latest set of predictions (after looking into proteomics data)

best_features <- read.csv("data/selected_features/best_features_with_is_best.csv") %>%
  filter(is_best == 1, grepl("CF_EV_tra_combat_", dataset_id)) %>%
  separate(dataset_id, into = c(NA, NA, NA, NA, "comparison"), remove = FALSE) %>%
  mutate(dataset_id = paste(dataset_id, description, 
                            min_iter_feature_presence, comparison,
                            sep = "_")) %>%
  mutate(cm = c("Sigmoid Kernel SVM", "L2 Regularized logistic regression", "Random Forest")) %>%
  dplyr::select(c(dataset_id, cm, comparison))

sample_wise_results <- read.csv("fem_pipeline_results_tra_combat_subset/all_samplewise_result_df.csv") %>%
  inner_join(best_features, by = c("DataSetId" = "dataset_id", "Model" = "cm")) %>%
  mutate(Model = sub("L2 Regularized logistic regression", "L2 Regularized\nlogistic regression", Model)) %>%
  mutate(Model = sub("Sigmoid Kernel SVM", "Sigmoid\nKernel SVM", Model))

phenotype <- read.table("data/formatted/phenotype.txt", sep = "\t", header = TRUE)


levels(factor(sample_wise_results$comparison))



create_all_iter_pred_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsIGT",
                             classes = c("IGT", "CFRD"),
                             sample_type = "test", phenotype,
                             plot_dir_path = "plots_updated/custom_heatmap/per_iter_pred/")
create_all_iter_pred_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsNGT",
                             classes = c("NGT", "CFRD"),
                             sample_type = "test", phenotype,
                             plot_dir_path = "plots_updated/custom_heatmap/per_iter_pred/")
create_all_iter_pred_heatmap(sample_wise_results,
                             comparison_of_interest = "IGTVsNGT",
                             classes = c("NGT", "IGT"),
                             sample_type = "test", phenotype,
                             plot_dir_path = "plots_updated/custom_heatmap/per_iter_pred/")

create_all_iter_prob_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsIGT",
                             classes = c("IGT", "CFRD"),
                             sample_type = "test", phenotype,
                             plot_dir_path = "plots_updated/custom_heatmap/per_iter_prob/")
create_all_iter_prob_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsNGT",
                             classes = c("NGT", "CFRD"),
                             sample_type = "test", phenotype,
                             plot_dir_path = "plots_updated/custom_heatmap/per_iter_prob/")
create_all_iter_prob_heatmap(sample_wise_results,
                             comparison_of_interest = "IGTVsNGT",
                             classes = c("NGT", "IGT"),
                             sample_type = "test", phenotype,
                             plot_dir_path = "plots_updated/custom_heatmap/per_iter_prob/")


create_all_iter_pred_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsIGT",
                             classes = c("IGT", "CFRD"),
                             sample_type = "train", phenotype,
                             plot_dir_path = "plots_updated/custom_heatmap/per_iter_pred/")
create_all_iter_pred_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsNGT",
                             classes = c("NGT", "CFRD"),
                             sample_type = "train", phenotype,
                             plot_dir_path = "plots_updated/custom_heatmap/per_iter_pred/")
create_all_iter_pred_heatmap(sample_wise_results,
                             comparison_of_interest = "IGTVsNGT",
                             classes = c("NGT", "IGT"),
                             sample_type = "train", phenotype,
                             plot_dir_path = "plots_updated/custom_heatmap/per_iter_pred/")

create_all_iter_prob_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsIGT",
                             classes = c("IGT", "CFRD"),
                             sample_type = "train", phenotype,
                             plot_dir_path = "plots_updated/custom_heatmap/per_iter_prob/")
create_all_iter_prob_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsNGT",
                             classes = c("NGT", "CFRD"),
                             sample_type = "train", phenotype,
                             plot_dir_path = "plots_updated/custom_heatmap/per_iter_prob/")
create_all_iter_prob_heatmap(sample_wise_results,
                             comparison_of_interest = "IGTVsNGT",
                             classes = c("NGT", "IGT"),
                             sample_type = "train", phenotype,
                             plot_dir_path = "plots_updated/custom_heatmap/per_iter_prob/")






########################

#heatmaps for proteomics

best_features <- read.csv("data/selected_features/best_features_with_is_best.csv") %>%
  filter(is_best == 1, grepl("CF_EV_prot_mf_quantile_combat_", dataset_id)) %>%
  separate(dataset_id, into = c(NA, NA, NA, NA, NA, NA, "comparison"), remove = FALSE) %>%
  mutate(dataset_id = paste(dataset_id, description, 
                            min_iter_feature_presence, comparison,
                            sep = "_")) %>%
  mutate(cm = c("Radial Kernel SVM", "Random Forest", "Random Forest")) %>%
  dplyr::select(c(dataset_id, cm, comparison))

sample_wise_results <- read.csv("fem_pipeline_results_prot_mf_quantile_combat_subset/all_samplewise_result_df.csv") %>%
  inner_join(best_features, by = c("DataSetId" = "dataset_id", "Model" = "cm")) %>%
  mutate(Model = sub("Radial Kernel SVM", "Radial\nKernel SVM", Model))

phenotype <- read.table("data/formatted/prot_phenotype_333.txt", sep = "\t", header = TRUE)

levels(factor(sample_wise_results$comparison))
levels(factor(sample_wise_results$DataSetId))

create_all_iter_pred_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsIGT",
                             classes = c("IGT", "CFRD"),
                             sample_type = "test", phenotype,
                             plot_dir_path = "plots_updated/custom_heatmap_prot/per_iter_pred/")
create_all_iter_pred_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsNGT",
                             classes = c("NGT", "CFRD"),
                             sample_type = "test", phenotype,
                             plot_dir_path = "plots_updated/custom_heatmap_prot/per_iter_pred/")
create_all_iter_pred_heatmap(sample_wise_results,
                             comparison_of_interest = "IGTVsNGT",
                             classes = c("NGT", "IGT"),
                             sample_type = "test", phenotype,
                             plot_dir_path = "plots_updated/custom_heatmap_prot/per_iter_pred/")

create_all_iter_prob_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsIGT",
                             classes = c("IGT", "CFRD"),
                             sample_type = "test", phenotype,
                             plot_dir_path = "plots_updated/custom_heatmap_prot/per_iter_prob/")
create_all_iter_prob_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsNGT",
                             classes = c("NGT", "CFRD"),
                             sample_type = "test", phenotype,
                             plot_dir_path = "plots_updated/custom_heatmap_prot/per_iter_prob/")
create_all_iter_prob_heatmap(sample_wise_results,
                             comparison_of_interest = "IGTVsNGT",
                             classes = c("NGT", "IGT"),
                             sample_type = "test", phenotype,
                             plot_dir_path = "plots_updated/custom_heatmap_prot/per_iter_prob/")


create_all_iter_pred_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsIGT",
                             classes = c("IGT", "CFRD"),
                             sample_type = "train", phenotype,
                             plot_dir_path = "plots_updated/custom_heatmap_prot/per_iter_pred/")
create_all_iter_pred_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsNGT",
                             classes = c("NGT", "CFRD"),
                             sample_type = "train", phenotype,
                             plot_dir_path = "plots_updated/custom_heatmap_prot/per_iter_pred/")
create_all_iter_pred_heatmap(sample_wise_results,
                             comparison_of_interest = "IGTVsNGT",
                             classes = c("NGT", "IGT"),
                             sample_type = "train", phenotype,
                             plot_dir_path = "plots_updated/custom_heatmap_prot/per_iter_pred/")

create_all_iter_prob_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsIGT",
                             classes = c("IGT", "CFRD"),
                             sample_type = "train", phenotype,
                             plot_dir_path = "plots_updated/custom_heatmap_prot/per_iter_prob/")
create_all_iter_prob_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsNGT",
                             classes = c("NGT", "CFRD"),
                             sample_type = "train", phenotype,
                             plot_dir_path = "plots_updated/custom_heatmap_prot/per_iter_prob/")
create_all_iter_prob_heatmap(sample_wise_results,
                             comparison_of_interest = "IGTVsNGT",
                             classes = c("NGT", "IGT"),
                             sample_type = "train", phenotype,
                             plot_dir_path = "plots_updated/custom_heatmap_prot/per_iter_prob/")


########################

#2024 Sept - create samplewise prediction heatmaps

#for prot same as previous set of calls - but creating in a different directory now

best_features <- read.csv("data/selected_features/best_features_with_is_best.csv") %>%
  filter(is_best == 1, grepl("CF_EV_prot_mf_quantile_combat_", dataset_id)) %>%
  separate(dataset_id, into = c(NA, NA, NA, NA, NA, NA, "comparison"), remove = FALSE) %>%
  mutate(dataset_id = paste(dataset_id, description, 
                            min_iter_feature_presence, comparison,
                            sep = "_")) %>%
  mutate(cm = c("Radial Kernel SVM", "Random Forest", "Random Forest")) %>%
  dplyr::select(c(dataset_id, cm, comparison))

sample_wise_results <- read.csv("fem_pipeline_results_prot_mf_quantile_combat_subset/all_samplewise_result_df.csv") %>%
  inner_join(best_features, by = c("DataSetId" = "dataset_id", "Model" = "cm")) %>%
  mutate(Model = sub("Radial Kernel SVM", "Radial\nKernel SVM", Model))

phenotype <- read.table("data/formatted/prot_phenotype_333.txt", sep = "\t", header = TRUE)

levels(factor(sample_wise_results$comparison))
levels(factor(sample_wise_results$DataSetId))

create_all_iter_pred_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsIGT",
                             classes = c("IGT", "CFRD"),
                             sample_type = "test", phenotype,
                             plot_dir_path = "plots_Sep2024/per_iter_pred_prot/",
                             pdf_plot = TRUE)
create_all_iter_pred_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsNGT",
                             classes = c("NGT", "CFRD"),
                             sample_type = "test", phenotype,
                             plot_dir_path = "plots_Sep2024/per_iter_pred_prot/",
                             pdf_plot = TRUE)
create_all_iter_pred_heatmap(sample_wise_results,
                             comparison_of_interest = "IGTVsNGT",
                             classes = c("NGT", "IGT"),
                             sample_type = "test", phenotype,
                             plot_dir_path = "plots_Sep2024/per_iter_pred_prot/",
                             pdf_plot = TRUE)
create_all_iter_pred_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsIGT",
                             classes = c("IGT", "CFRD"),
                             sample_type = "train", phenotype,
                             plot_dir_path = "plots_Sep2024/per_iter_pred_prot/",
                             pdf_plot = TRUE)
create_all_iter_pred_heatmap(sample_wise_results,
                             comparison_of_interest = "CFRDVsNGT",
                             classes = c("NGT", "CFRD"),
                             sample_type = "train", phenotype,
                             plot_dir_path = "plots_Sep2024/per_iter_pred_prot/",
                             pdf_plot = TRUE)
create_all_iter_pred_heatmap(sample_wise_results,
                             comparison_of_interest = "IGTVsNGT",
                             classes = c("NGT", "IGT"),
                             sample_type = "train", phenotype,
                             plot_dir_path = "plots_Sep2024/per_iter_pred_prot/",
                             pdf_plot = TRUE)
