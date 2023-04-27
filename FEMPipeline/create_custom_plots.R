library(tidyverse)
library(readxl)
library(ROCR)

library(ggvenn)


base_dir <- "/home/abhivij/UNSW/VafaeeLab/CysticFibrosisGroup/ExoCF/CFRD_EV_biomarker/"
setwd(base_dir)


best_features <- read.csv("data/selected_features/best_features_with_is_best.csv")   %>%
  filter(is_best == 1, grepl("CF_EV_adult_filtered_seurat3_norm_find_var_none_", dataset_id)) %>%
  separate(dataset_id, into = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, "comparison"), remove = FALSE) %>%
  mutate(dataset_id = paste(dataset_id, description, 
                            min_iter_feature_presence, comparison,
                            sep = "_")) %>%
  mutate(cm = c("Random Forest", "L2 Regularized logistic regression", "Sigmoid Kernel SVM")) %>%
  dplyr::select(c(dataset_id, cm, comparison))

sample_wise_results <- read.csv("fem_pipeline_results_adult_filtered_then_seurat3_norm_find_var_none_subset/all_samplewise_result_df.csv") %>%
  inner_join(best_features, by = c("DataSetId" = "dataset_id", "Model" = "cm"))


# sample_wise_results.train <- sample_wise_results %>%
#   filter(Type == "train")
# sample_wise_results.test <- sample_wise_results %>%
#   filter(Type == "test")

current_comparison <- "CFRDVsNGT"
classes <- c("NGT", "CFRD")

sample_type <- "test"
dir_path <- "prediction_pipeline/plots/ROC/"

plot_roc <- function(sample_wise_results, sample_type, current_comparison,
                     classes,
                     dir_path = "prediction_pipeline/plots/ROC/"){
  results <- sample_wise_results %>%
    filter(Type == sample_type, comparison == current_comparison)
  pred_prob <- list()
  true_label <- list()
  for(i in c(1:30)){
    pred_prob[[i]] <- results[results$Iter == i, "PredProb"]
    true_label[[i]] <- results[results$Iter == i, "TrueLabel"]
  }
  
  pred <- prediction(pred_prob, true_label, label.ordering = classes)
  perf <- performance(pred, "tpr", "fpr")
  
  title <- paste(sub("Vs", " Vs ", current_comparison), 
                 ": ROC curves for the 30 sets of", sample_type, "samples")
  if(!dir.exists(dir_path)){
    dir.create(path = dir_path, recursive = TRUE)
  }
  file_name <- paste0(current_comparison, "_", sample_type, ".png")
  png(paste0(dir_path, file_name))
  
  plot(perf,
       lty = 3, lwd = 0.5,
       colorize = TRUE,
       main = title)
  # plot(perf,
  #      avg = 'horizontal', spread.estimate = 'stderr',
  #      lwd = 2, col = 1, add = TRUE)
  # plot(perf,
  #      avg = 'vertical', spread.estimate = 'boxplot',
  #      lwd = 2, col = 2, add = TRUE)
  
  pred_prob <- results[, "PredProb"]
  true_label <- results[, "TrueLabel"]
  pred <- prediction(pred_prob, true_label, label.ordering = classes)
  perf <- performance(pred, "tpr", "fpr")
  
  plot(perf,
       lwd = 2, add = TRUE)
  dev.off()
}

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
