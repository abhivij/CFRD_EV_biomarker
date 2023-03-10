#custom pipeline to test best biomarkers and best models on child pre-post modulator samples

library(tidyverse)
library(umap)

base_dir <- "~/UNSW/VafaeeLab/CysticFibrosisGroup/ExoCF/CFRD_EV_biomarker/"
setwd(base_dir)


source("prediction_pipeline/cm_logistic_regression.R")
source("prediction_pipeline/cm_svm.R")
source("prediction_pipeline/cm_rf.R")

######################################################################


dataset_replace_str = "CF_EV_adult_filtered_seurat3_norm_find_var_none_"
train_data_file_path <- "data/formatted/umi_counts_filtered_seurat3_with_norm_and_find_var_feat.csv"
test_data_file_path <- "data/formatted/umi_counts_prepost_filtered_seurat3.csv"

comparison = "CFRDVsIGT"
classes = c("IGT", "CFRD")
best_features_file_path = "data/selected_features/best_features_with_is_best.csv"

comparison = "CFRDVsNGT"
classes = c("NGT", "CFRD")

result_file_dir_path = "data/prediction_result_pre_post/"
result_file_name = "CFRDVsIGT.csv"

comparison = "IGTVsNGT"
classes = c("NGT", "IGT")

result_file_name = "IGTVsNGT.csv"

data <- read.table("data/formatted/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                   nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")



custom_pipeline <- function(comparison, classes, 
                            best_features_file_path,
                            class_colours,
                            dataset_replace_str = "CF_EV_adult_filtered_seurat3_norm_find_var_none_",
                            result_file_dir_path = "data/prediction_result_pre_post/",
                            result_file_name = "CFRDVsIGT.csv",
                            model = "rf",
                            train_data_file_path = "data/formatted/umi_counts_filtered_seurat3_with_norm_and_find_var_feat.csv",
                            test_data_file_path = "data/formatted/umi_counts_prepost_filtered_seurat3.csv"){
  
  best_features <- read.csv(best_features_file_path)  
  best_features_sub <- best_features %>%
    mutate(dataset_id = gsub(dataset_replace_str, "", dataset_id)) %>%
    filter(is_best == 1, dataset_id == comparison)
  
  biomarkers <- strsplit(best_features_sub$biomarkers, split = "|", fixed = TRUE)[[1]]  
  
  train.data <- read.table(train_data_file_path, header=TRUE, sep=",", row.names=1, skip=0,
                           nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
  test.data <- read.table(test_data_file_path, header=TRUE, sep=",", row.names=1, skip=0,
                          nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")

  train.phenotype <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t") %>%
    filter(age_group == "adult") %>%
    filter(!is.na(CFRDVsIGT) | !is.na(CFRDVsNGT) | !is.na(IGTVsNGT))
  
  test.phenotype <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t") %>%
    filter(age_group == "child" & condition == "CF_pre_post_modulator")

  
  train.output_labels <- train.phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes) %>%
    dplyr::select(Sample, Label) %>%
    dplyr::mutate(Label = factor(Label)) %>%
    arrange(Label, Sample)
  train.data <- train.data[, train.output_labels$Sample]

  
  test.output_labels <- test.phenotype %>%
    rename("Label" = condition) %>%
    dplyr::select(Sample, Label, pre_post_modulator, modulator) %>%
    dplyr::mutate(Label = factor(Label)) %>%    
    arrange(pre_post_modulator, Sample)
  test.data <- test.data[, test.output_labels$Sample]

  #currently train.data, test.data format : (transcripts x samples)
  train.data <- as.data.frame(t(as.matrix(train.data)))
  test.data <- as.data.frame(t(as.matrix(test.data)))  
  #now train.data, test.data format : (samples x transcripts)
  
  create_dim_red_plot_adult_and_child_data(train.data, test.data, 
                                           train.output_labels, test.output_labels,
                                           paste0(classes[2], ", ",
                                                  classes[1], 
                                                  " and pre-post modulator samples with",
                                                  " all filtered transcripts"),
                                           plot_dir_path = "prediction_pipeline/plots/dim_red_with_prepost/",
                                           plot_file_name = paste0(classes[2], "_", classes[1],
                                                                   "_all_filtered.jpg"),
                                           class_colours)
  
  #get best biomarkers only
  features_with_slash <- colnames(train.data)[grepl("/", colnames(train.data), fixed = TRUE)] 
  for(f in features_with_slash){
    f_replaced <- gsub("/|-", ".", f) 
    if(f_replaced %in% biomarkers){
      biomarkers[biomarkers == f_replaced] = f
    }
  }
  biomarkers <- gsub(".", "-", biomarkers, fixed = TRUE)
  
  train.data <- train.data[, biomarkers]
  test.data <- test.data[, biomarkers]
  
  create_dim_red_plot_adult_and_child_data(train.data, test.data, 
                                           train.output_labels, test.output_labels,
                                           paste0(classes[2], ", ",
                                                  classes[1], 
                                                  " and pre-post modulator samples with",
                                                  " all selected biomarkers"),
                                           plot_dir_path = "prediction_pipeline/plots/dim_red_with_prepost/",
                                           plot_file_name = paste0(classes[2], "_", classes[1],
                                                                   "_selected_biomarkers.jpg"),
                                           class_colours)
  
  if(model == "rf"){
    result_df <- rf_model(train.data, train.output_labels, test.data, test.output_labels, classes)
  } else if(model == "radial_svm"){
    result_df <- svm_model(train.data, train.output_labels,
                           test.data, test.output_labels,
                           classes, kernel = "radial")
  } else if(model == "sigmoid_svm"){
    result_df <- svm_model(train.data, train.output_labels,
                           test.data, test.output_labels,
                           classes, kernel = "sigmoid")
  } else if(model == "l2_log_reg"){
    result_df <- log_reg_model(train.data, train.output_labels, test.data, test.output_labels,
                               classes, regularize = "l2")
  } else if(model == "l1_log_reg"){
    result_df <- log_reg_model(train.data, train.output_labels, test.data, test.output_labels,
                               classes, regularize = "l1")
  } else if(model == "log_reg"){
    result_df <- log_reg_model(train.data, train.output_labels, test.data, test.output_labels,
                               classes)
  }
  
  result_df <- result_df %>%
    left_join(test.output_labels, by = c("sample" = "Sample"))
  
  if(!dir.exists(result_file_dir_path)){
    dir.create(result_file_dir_path, recursive = TRUE)
  }
  write.csv(format(result_df, digits = 3), paste0(result_file_dir_path, result_file_name), row.names = FALSE)
}  

# class_colours = c("red", "orange", "green")
# title <- "CFRD, IGT and pre-post modulator samples with all filtered transcripts"

create_dim_red_plot_adult_and_child_data <- function(train.data, test.data,
                                                     train.output_labels, test.output_labels,
                                                     title,
                                                     plot_dir_path,
                                                     plot_file_name,
                                                     class_colours){
  
  data <- rbind(train.data, test.data)
  output_labels <- rbind(train.output_labels %>% 
                           mutate("pre_post_modulator" = "N/A", "modulator"= "N/A"),
                         test.output_labels) %>%
    mutate(modulator = ifelse(is.na(modulator), "N/A", modulator)) %>%
    mutate(pre_post_modulator = factor(case_when(
                                        pre_post_modulator == "0" ~ "pre-modulator",
                                        pre_post_modulator == "1" ~ "post-modulator",
                                        TRUE ~ pre_post_modulator),
                                levels = c("pre-modulator", "post-modulator", "N/A"))
           )
  
  data <- data[output_labels$Sample, ]
  
  all.equal(rownames(data), output_labels$Sample)
  #TRUE
  
  result <- umap(data)
  dim_red_df <- data.frame(x = result$layout[,1], y = result$layout[,2])  
  xlab <- "UMAP 1"
  ylab <- "UMAP 2"
  
  #note overall point diam is size + stroke
  #so total : 3+1 = 4
  #hence in legend, with 0 stroke, 4 size is given
  ggplot2::ggplot(dim_red_df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(ggplot2::aes(fill = output_labels$Label,
                                     shape = output_labels$pre_post_modulator,
                                     colour = output_labels$modulator), 
                        size = 3, stroke = 1) +
    ggplot2::guides(fill = guide_legend(override.aes = list(shape = 21, stroke = 0, size = 4))) +
    ggplot2::scale_fill_manual(name = "Condition", values = class_colours) +
    ggplot2::guides(shape = guide_legend(override.aes = list(fill = c("black", "black", "black")))) +
    ggplot2::scale_shape_manual(name = "pre/post modulator", values = c(21, 22, 23)) +
    ggplot2::guides(colour = guide_legend(override.aes = list(shape = 1))) +
    ggplot2::scale_colour_manual(name = "Modulator", values = c("black", "purple", "royalblue1", "cyan")) +
    ggplot2::labs(title = title,
                  caption = paste("Data dimension :", paste(dim(data), collapse = "x"))) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab)
  
  if(!dir.exists(plot_dir_path)){
    dir.create(plot_dir_path, recursive = TRUE)
  }
  file_path <- paste(plot_dir_path, plot_file_name, sep = "/")
  ggplot2::ggsave(file_path, units = "cm", width = 30)
}


show_metrics <- function(comparison, classes, result_file_path,
                         metric_output_file_path = "data/prediction_result/metrics.csv"){
  print(comparison)
  result_df <- read.csv(result_file_path)  
  
  results.train <- result_df %>%
    filter(Type == "train")
  results.test <- result_df %>%
    filter(Type == "test")
  
  acc.train <- mean(results.train$TrueLabel == results.train$PredictedLabel)
  acc.test <- mean(results.test$TrueLabel == results.test$PredictedLabel)
  
  pr <- ROCR::prediction(results.train$Pred_prob, results.train$TrueLabel, label.ordering = classes)
  auc.train <- ROCR::performance(pr, measure = "auc")@y.values[[1]]
  
  pr <- ROCR::prediction(results.test$Pred_prob, results.test$TrueLabel, label.ordering = classes)
  auc.test <- ROCR::performance(pr, measure = "auc")@y.values[[1]]
  
  #classes[2] assigned to class1 deliberately
  #classes are provided as c(neg_class, pos_class)
  class1 <- classes[2]
  class2 <- classes[1]
  class1_count <- results.train %>% dplyr::filter(TrueLabel == class1) %>% nrow()
  class2_count <- results.train %>% dplyr::filter(TrueLabel == class2) %>% nrow() 
  train_count <- paste0(class1, "=", class1_count, " | ", class2, "=", class2_count)
  class1_count <- results.test %>% dplyr::filter(TrueLabel == class1) %>% nrow()
  class2_count <- results.test %>% dplyr::filter(TrueLabel == class2) %>% nrow() 
  test_count <- paste0(class1, "=", class1_count, " | ", class2, "=", class2_count)
  metrics <- data.frame(Comparison = comparison,
                        TrainCount = train_count,
                        TestCount = test_count,
                        TrainAccuracy = acc.train,
                        TrainAUC = auc.train,
                        TestAccuracy = acc.test,
                        TestAUC = auc.test)
  
  write.table(x = metrics, file = metric_output_file_path, append = TRUE, 
              col.names = !file.exists(metric_output_file_path), sep = ",",
              row.names = FALSE)
}

custom_pipeline(
  comparison = "CFRDVsIGT",
  classes = c("IGT", "CFRD"),
  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
  class_colours = c("red", "orange", "green"),
  dataset_replace_str = "CF_EV_adult_filtered_seurat3_norm_find_var_none_",
  result_file_dir_path = "data/prediction_result_pre_post/",
  result_file_name = "CFRDVsIGT.csv",
  model = "rf"
)

custom_pipeline(
  comparison = "CFRDVsNGT",
  classes = c("NGT", "CFRD"),
  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
  class_colours = c("red", "yellow", "green"),
  dataset_replace_str = "CF_EV_adult_filtered_seurat3_norm_find_var_none_",
  result_file_dir_path = "data/prediction_result_pre_post/",
  result_file_name = "CFRDVsNGT.csv",
  model = "l2_log_reg"
)

custom_pipeline(
  comparison = "IGTVsNGT",
  classes = c("NGT", "IGT"),
  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
  class_colours = c("orange", "yellow", "green"),
  dataset_replace_str = "CF_EV_adult_filtered_seurat3_norm_find_var_none_",
  result_file_dir_path = "data/prediction_result_pre_post/",
  result_file_name = "IGTVsNGT.csv",
  model = "sigmoid_svm"
)

