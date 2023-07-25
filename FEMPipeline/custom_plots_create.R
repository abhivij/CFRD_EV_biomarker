library(tidyverse)
library(readxl)
library(ROCR)
library(pROC)

library(ggvenn)
library(ComplexHeatmap)
library(viridis)
library(RColorBrewer)




plot_roc <- function(sample_wise_results, sample_type, current_comparison,
                     classes,
                     dir_path = "prediction_pipeline/plots/ROC/"){
  results <- sample_wise_results %>%
    filter(Type == sample_type, comparison == current_comparison)
  # pred_prob <- list()
  # true_label <- list()
  # for(i in c(1:30)){
  #   pred_prob[[i]] <- results[results$Iter == i, "PredProb"]
  #   true_label[[i]] <- results[results$Iter == i, "TrueLabel"]
  # }
  # 
  # pred <- prediction(pred_prob, true_label, label.ordering = classes)
  # perf <- performance(pred, "tpr", "fpr")
  
  title <- paste(sub("Vs", " Vs ", current_comparison), 
                 ": ROC curve for the 30 sets of", sample_type, "samples")
  if(!dir.exists(dir_path)){
    dir.create(path = dir_path, recursive = TRUE)
  }
  file_name <- paste0(current_comparison, "_", sample_type, ".png")
  png(paste0(dir_path, file_name))
  # 
  # plot(perf,
  #      lty = 3, lwd = 0.5,
  #      colorize = TRUE,
  #      main = title)
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
       lwd = 2, main = title)
  dev.off()
}


plot_roc <- function(sample_wise_results, sample_type, current_comparison,
                     classes,
                     dir_path = "prediction_pipeline/plots/ROC/"){
  results <- sample_wise_results %>%
    filter(Type == sample_type, comparison == current_comparison)

  title <- paste(sub("Vs", " Vs ", current_comparison), 
                 ": ROC curves for all iter combined", sample_type, "samples")
  if(!dir.exists(dir_path)){
    dir.create(path = dir_path, recursive = TRUE)
  }
  file_name <- paste0(current_comparison, "_", sample_type, ".png")
  png(paste0(dir_path, file_name))
  
  roc_curve <- roc(results[, "TrueLabel"],
                   1-results[, "PredProb"], levels = classes, direction = '<')
  print(results[, "PredProb"][1:10])
  print(1-results[, "PredProb"][1:10])
  plot(roc_curve, main = title) 
  
  # for(i in c(1:30)){
  #   if(i == 1){
  #     print(results[results$Iter == i, "PredProb"])
  #     print(1-results[results$Iter == i, "PredProb"])
  #     roc_curve <- roc(results[results$Iter == i, "TrueLabel"],
  #                      1-results[results$Iter == i, "PredProb"], levels = rev(classes))
  #     plot(roc_curve, col = i, main = title)      
  #   } else{
  #     #plot(roc_curve, col = i, add = TRUE)
  #   }
  # }
  dev.off()
}


# sample_wise_results
# sample_type = "test"
# comparison_of_interest = "CFRDVsIGT"
# classes = c("IGT", "CFRD")
# phenotype <- read.table("data/formatted/phenotype.txt", sep = "\t", header = TRUE)





#classes assumed to be given as (neg_class, pos_class)
create_mean_prob_heatmap <- function(sample_wise_results,
                                     comparison_of_interest,
                                     classes,
                                     sample_type){
  results <- sample_wise_results %>%
    filter(Type == sample_type, comparison == comparison_of_interest) %>%
    dplyr::select(-c(DataSetId, comparison, Type)) %>%
    group_by(Sample, TrueLabel, Model) %>%
    summarize(mean_prob = mean(PredProb))
  
  sample_label <- results %>%
    dplyr::select(c(Sample, TrueLabel)) %>%
    distinct() %>%
    arrange(Sample)
  
  results <- results %>%
    ungroup()
  
  data_to_plot <- results %>%
    dplyr::select(c(Sample, Model, mean_prob)) %>%
    pivot_wider(names_from = Model, values_from = mean_prob) %>%
    arrange(Sample) %>%
    column_to_rownames("Sample")
  data_to_plot <- data.matrix(data_to_plot)
  # 
  # model_names <- c("Simple logistic regression",
  #                  "Radial Kernel SVM",
  #                  "Sigmoid Kernel SVM",
  #                  "Random Forest",
  #                  "L1 Regularized logistic regression",
  #                  "L2 Regularized logistic regression",
  #                  "Elastic net logistic regression")
  
  meta_data.row <- data.frame(Sample = rownames(data_to_plot)) %>% 
    inner_join(sample_label) %>%
    mutate(TrueLabel = factor(TrueLabel, levels = rev(classes)))
  meta_data.col <- data.frame(model = colnames(data_to_plot))
  
  row_col <- list()
  row_col[["TrueLabel"]] <- c("#0d0887", "#f0f921")
  names(row_col[["TrueLabel"]]) <- classes
  # column_col <- list("Omics type" = c("proteomics" = "skyblue1",
  #                                     "transcriptomics" = "indianred1"))
  # column_col[["Model"]] <- brewer.pal(n = 7, name = "Paired")
  # names(column_col[["Model"]]) <- model_names
  # 
  Heatmap(data_to_plot, name = "Mean Prediction probability",
                col = plasma(5),
                rect_gp = gpar(col = "white", lwd = 1),
                cluster_columns = FALSE,
                cluster_rows = FALSE,
                row_title = "Samples",
                row_names_side = "left",
                row_split = meta_data.row$TrueLabel,
                column_split = meta_data.col$model,
                column_title = NULL,
                show_column_names = FALSE,
                row_names_gp = gpar(fontsize = 5),
                bottom_annotation = HeatmapAnnotation(
                  "Model" = meta_data.col$model,
                  annotation_name_side = "left",
                  col = list("Model" = c("Sigmoid Kernel SVM" = "skyblue1"))
                )) + 
    HeatmapAnnotation("TrueLabel" = meta_data.row$TrueLabel, 
                      which = "row", 
                      col = row_col)
  
  png(paste0("plots/custom_heatmap/", comparison_of_interest, 
             sample_type,
             ".png"), units = "cm", width = 20, height = 15, res = 1200)  
  draw(ht, column_title = paste(sample_type, sub("Vs", " Vs ", comparison_of_interest)))    
  dev.off() 
  
}


# plot_dir_path = "plots/custom_heatmap/per_iter_prob/"

#classes assumed to be given as (neg_class, pos_class)
create_all_iter_prob_heatmap <- function(sample_wise_results,
                                    comparison_of_interest,
                                    classes,
                                    sample_type, phenotype,
                                    plot_dir_path = "plots/custom_heatmap/per_iter_prob/"){
  results <- sample_wise_results %>%
    filter(Type == sample_type, comparison == comparison_of_interest)%>%
    dplyr::select(c(Sample, Iter, PredProb, PredictedLabel, TrueLabel, Model))

  data_to_plot <- results %>%
    dplyr::select(c(Sample, Iter, PredProb)) %>%
    pivot_wider(names_from = Iter, values_from = PredProb) %>%
    arrange(Sample) %>%
    column_to_rownames("Sample")
  data_to_plot <- data.matrix(data_to_plot)
  # 
  # model_names <- c("Simple logistic regression",
  #                  "Radial Kernel SVM",
  #                  "Sigmoid Kernel SVM",
  #                  "Random Forest",
  #                  "L1 Regularized logistic regression",
  #                  "L2 Regularized logistic regression",
  #                  "Elastic net logistic regression")
  
  meta_data.row <- data.frame(Sample = rownames(data_to_plot)) %>% 
    inner_join(phenotype %>% 
                 dplyr::select(Sample, cohort, age_group, sex, condition) %>%
                 dplyr::rename("TrueLabel" = "condition")) %>%
    mutate(TrueLabel = factor(TrueLabel, levels = rev(classes)))
  meta_data.col <- data.frame(Iter = colnames(data_to_plot)) %>%
    mutate(Model = unique(results$Model)) %>%
    mutate(Repeat = paste("Repeat", ceiling(as.numeric(Iter)/5)))
  
  row_col <- list()
  row_col[["TrueLabel"]] <- c("#0d0887", "#f0f921")
  names(row_col[["TrueLabel"]]) <- classes
  
  row_col[["Sex"]] <- c("M" = "steelblue1", "F" = "plum1")
  
  cohorts <- unique(meta_data.row$cohort)
  row_col[["Cohort"]] <- brewer.pal(n = 8, name = "Paired")[1:length(cohorts)]
  names(row_col[["Cohort"]]) <- cohorts
  
  age_groups <- unique(meta_data.row$age_group)
  row_col[["Age Group"]] <- brewer.pal(n = 8, name = "Set2")[1:length(age_groups)]
  names(row_col[["Age Group"]]) <- age_groups
  
  column_col <- list()
  column_col[["Model"]] <- c("skyblue1")
  names(column_col[["Model"]]) <- unique(results$Model)
  
  column_col[["CV repeats"]] <- brewer.pal(n = 6, name = "Paired")
  names(column_col[["CV repeats"]]) <- paste("Repeat", c(1:6))
   
  ht <- Heatmap(data_to_plot, name = "Prediction probability\nin each iteration",
          col = plasma(5),
          rect_gp = gpar(col = "white", lwd = 1),
          cluster_columns = FALSE,
          cluster_rows = FALSE,
          row_title = "Samples",
          row_names_side = "left",
          row_split = meta_data.row$TrueLabel,
          column_split = meta_data.col$Repeat,
          row_names_gp = gpar(fontsize = 5),
          row_title_gp = gpar(fontsize = 9),
          column_names_gp = gpar(fontsize = 4),
          column_names_rot = 45,
          column_title = "Repeated CV iterations",
          column_title_side = "top",
          column_title_gp = gpar(fontsize = 9),
          column_names_side = "top",
          bottom_annotation = HeatmapAnnotation(
            "CV repeats" = meta_data.col$Repeat,
            "Model" = meta_data.col$Model,
            annotation_name_side = "left",
            col = column_col
          )) + 
    HeatmapAnnotation("TrueLabel" = meta_data.row$TrueLabel, 
                      "Cohort" = meta_data.row$cohort,
                      "Age Group" = meta_data.row$age_group,
                      "Sex" = meta_data.row$sex,
                      which = "row", 
                      col = row_col)
  
  if(!dir.exists(plot_dir_path)){
    dir.create(plot_dir_path, recursive = TRUE)
  }
  png(paste0(plot_dir_path, comparison_of_interest, 
             sample_type,
             ".png"), units = "cm", width = 20, height = 15, res = 1200)  
  draw(ht, 
       column_title = paste("Probabilities : ", sample_type, "subset samples of", sub("Vs", " Vs ", comparison_of_interest)))    
  dev.off() 
  
}




create_all_iter_pred_heatmap <- function(sample_wise_results,
                                         comparison_of_interest,
                                         classes,
                                         sample_type, phenotype,
                                         plot_dir_path = "plots/custom_heatmap/per_iter_pred/"){
  results <- sample_wise_results %>%
    filter(Type == sample_type, comparison == comparison_of_interest)%>%
    dplyr::select(c(Sample, Iter, PredProb, PredictedLabel, TrueLabel, Model))
  
  data_to_plot <- results %>%
    dplyr::select(c(Sample, Iter, PredictedLabel)) %>%
    pivot_wider(names_from = Iter, values_from = PredictedLabel) %>%
    arrange(Sample) %>%
    column_to_rownames("Sample")
  data_to_plot <- as.matrix(data_to_plot, na.strings = "NA")
  # 
  # model_names <- c("Simple logistic regression",
  #                  "Radial Kernel SVM",
  #                  "Sigmoid Kernel SVM",
  #                  "Random Forest",
  #                  "L1 Regularized logistic regression",
  #                  "L2 Regularized logistic regression",
  #                  "Elastic net logistic regression")
  
  meta_data.row <- data.frame(Sample = rownames(data_to_plot)) %>% 
    inner_join(phenotype %>% 
                 dplyr::select(Sample, cohort, age_group, sex, condition) %>%
                 dplyr::rename("TrueLabel" = "condition"))
  meta_data.col <- data.frame(Iter = colnames(data_to_plot)) %>%
    mutate(Model = unique(results$Model)) %>%
    mutate(Repeat = paste("Repeat", ceiling(as.numeric(Iter)/5)))
  
  row_col <- list()
  row_col[["True Label"]] <- c("#0d0887", "#f0f921")
  names(row_col[["True Label"]]) <- classes
  
  row_col[["Sex"]] <- c("M" = "steelblue1", "F" = "plum1")
  
  cohorts <- unique(meta_data.row$cohort)
  row_col[["Cohort"]] <- brewer.pal(n = 8, name = "Paired")[1:length(cohorts)]
  names(row_col[["Cohort"]]) <- cohorts
  
  age_groups <- unique(meta_data.row$age_group)
  row_col[["Age Group"]] <- brewer.pal(n = 8, name = "Set2")[1:length(age_groups)]
  names(row_col[["Age Group"]]) <- age_groups
  
  column_col <- list()
  column_col[["Model"]] <- c("skyblue1")
  names(column_col[["Model"]]) <- unique(results$Model)
  
  column_col[["CV repeats"]] <- brewer.pal(n = 6, name = "Paired")
  names(column_col[["CV repeats"]]) <- paste("Repeat", c(1:6))
  
  ht <- Heatmap(data_to_plot, name = "Predicted Label",
                col = rev(row_col[["True Label"]]),
                rect_gp = gpar(col = "white", lwd = 1),
                cluster_columns = FALSE,
                cluster_rows = FALSE,
                row_title = "Samples",
                row_names_side = "left",
                row_split = meta_data.row$TrueLabel,
                column_split = meta_data.col$Repeat,
                row_names_gp = gpar(fontsize = 5),
                row_title_gp = gpar(fontsize = 9),
                column_names_gp = gpar(fontsize = 4),
                column_names_rot = 45,
                column_title = "Repeated CV iterations",
                column_title_side = "top",
                column_title_gp = gpar(fontsize = 9),
                column_names_side = "top",
                bottom_annotation = HeatmapAnnotation(
                  "CV repeats" = meta_data.col$Repeat,
                  "Model" = meta_data.col$Model,
                  annotation_name_side = "left",
                  col = column_col
                )) + 
    HeatmapAnnotation("True Label" = factor(meta_data.row$TrueLabel, levels = rev(classes)), 
                      "Cohort" = meta_data.row$cohort,
                      "Age Group" = meta_data.row$age_group,
                      "Sex" = meta_data.row$sex,
                      which = "row", 
                      col = row_col)
  
  if(!dir.exists(plot_dir_path)){
    dir.create(plot_dir_path, recursive = TRUE)
  }
  png(paste0(plot_dir_path, comparison_of_interest, 
             sample_type,
             ".png"), units = "cm", width = 20, height = 15, res = 1200)  
  draw(ht, 
       column_title = paste("Predictions : ", sample_type, "subset samples of", sub("Vs", " Vs ", comparison_of_interest)))    
  dev.off() 
  
}
# 
# sample_wise_results
# sample_type = "test" 
# comparison_of_interest = "CFRDVsIGT"
# classes = c("IGT", "CFRD")
# phenotype <- read.table("data/formatted/phenotype.txt", sep = "\t", header = TRUE)
