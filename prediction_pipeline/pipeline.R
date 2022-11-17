library(tidyverse)

base_dir <- "~/UNSW/VafaeeLab/CysticFibrosisGroup/ExoCF/CFRD_EV_biomarker/"
setwd(base_dir)


source("prediction_pipeline/cm_logistic_regression.R")
source("prediction_pipeline/cm_svm.R")
source("prediction_pipeline/cm_rf.R")

######################################################################


comparison = "CFRDVsIGT"
classes = c("IGT", "CFRD")
best_features_file_path = "data/selected_features/best_features_with_is_best.csv"
only_adults = TRUE
dataset_replace_str = "CF_EV_AU_adult_logtmm_"
result_file_dir = "data/prediction_result_adult/"
result_file_name = "CFRDVsIGT.csv"
perform_filter = TRUE
norm = "log_tmm"

pipeline <- function(comparison, classes, 
                     best_features_file_path,
                     only_adults = FALSE,
                     dataset_replace_str = "CF_EV_AU_zlogtmm_",
                     result_file_dir = "data/prediction_result/",
                     result_file_name = "CFRDVsIGT.csv",
                     perform_filter = TRUE,
                     norm = "norm_log_tmm"){
  
  best_features <- read.csv(best_features_file_path)  
  best_features_sub <- best_features %>%
    mutate(dataset_id = gsub(dataset_replace_str, "", dataset_id)) %>%
    filter(is_best == 1, dataset_id == comparison)
  
  biomarkers <- strsplit(best_features_sub$biomarkers, split = "|", fixed = TRUE)[[1]]  
  
  data <- read.table("data/formatted/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                     nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
  phenotype <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t")
  
  if(only_adults){
    phenotype <- phenotype %>% filter(age_group == "adult")
  }
  
  output_labels.train <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes, country == "AU") %>%
    dplyr::select(Sample, Label, age, age_group, sex, FEV1) %>%
    dplyr::mutate(Label = factor(Label), age = as.numeric(age), 
                  age_group = factor(age_group), sex = factor(sex),
                  FEV1 = as.numeric(FEV1)) %>%
    arrange(Label, Sample)
  data.train <- data[, output_labels.train$Sample]
  print(summary(output_labels.train))

  output_labels.test <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes, country == "DK") %>%
    dplyr::select(Sample, Label, age, age_group, sex, FEV1) %>%
    dplyr::mutate(Label = factor(Label), age = as.numeric(age), 
                  age_group = factor(age_group), sex = factor(sex),
                  FEV1 = as.numeric(FEV1)) %>%    
    arrange(Label, Sample)
  data.test <- data[, output_labels.test$Sample]
  print(summary(output_labels.test))
  
  #currently data.train, data.test format : (transcripts x samples)
  
  if(perform_filter){
    keep <- edgeR::filterByExpr(data.train, group = output_labels.train$Label)
    data.train <- data.train[keep, ]
    data.test <- data.test[keep, ]  
  }
  
  
  if(norm == "norm_log_tmm"){
    #calculating norm log tmm
    dge <- edgeR::DGEList(counts = data.train, group = output_labels.train$Label)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    tmm <- edgeR::cpm(dge, log = TRUE)
    data.train <- tmm
    
    dge <- edgeR::DGEList(counts = data.test, group = output_labels.test$Label)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    tmm <- edgeR::cpm(dge, log = TRUE)
    data.test <- tmm
    
    data.train <- as.data.frame(t(as.matrix(data.train)))
    data.test <- as.data.frame(t(as.matrix(data.test)))  
    
    #normalizing the data
    normparam <- caret::preProcess(data.train) 
    data.train <- predict(normparam, data.train)
    data.test <- predict(normparam, data.test) #normalizing test data using params from train data   
    
  } else if(norm == "log_tmm"){
    #calculating norm log tmm
    dge <- edgeR::DGEList(counts = data.train, group = output_labels.train$Label)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    tmm <- edgeR::cpm(dge, log = TRUE)
    data.train <- tmm
    
    dge <- edgeR::DGEList(counts = data.test, group = output_labels.test$Label)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    tmm <- edgeR::cpm(dge, log = TRUE)
    data.test <- tmm
    
    data.train <- as.data.frame(t(as.matrix(data.train)))
    data.test <- as.data.frame(t(as.matrix(data.test)))  
  }
  #now data.train, data.test format : (samples x transcripts)
  
  #get best biomarkers only

  
  features_with_slash <- colnames(data.train)[grepl("/", colnames(data.train), fixed = TRUE)] 
  for(f in features_with_slash){
    f_replaced <- gsub("/|-", ".", f) 
    if(f_replaced %in% biomarkers){
      biomarkers[biomarkers == f_replaced] = f
    }
  }
  biomarkers <- gsub(".", "-", biomarkers, fixed = TRUE)
  
  
  data.train <- data.train[, biomarkers]
  data.test <- data.test[, biomarkers]

  #for old biomarkers identified - i.e. on AU adult+child with norm_log_tmm
  # if(comparison == "CFRDVsIGT"){
  #   #l2 log reg
  #   result_df <- log_reg_model(data.train, output_labels.train, data.test, output_labels.test, classes, regularize = "l2") 
  # } else if(comparison == "CFRDVsNGT"){
  #   #sigmoid kernel svm
  #   result_df <- svm_model(data.train, output_labels.train, data.test, output_labels.test, classes, kernel = "sigmoid")
  # } else if(comparison == "IGTVsNGT"){
  #   #random forest
  #   result_df <- rf_model(data.train, output_labels.train, data.test, output_labels.test, classes)
  # }
  
  if(comparison == "CFRDVsIGT"){
    #random forest
    result_df <- rf_model(data.train, output_labels.train, data.test, output_labels.test, classes)
  } else if(comparison == "CFRDVsNGT"){
    #radial kernel svm
    result_df <- svm_model(data.train, output_labels.train, data.test, output_labels.test, classes, kernel = "radial")
  } else if(comparison == "IGTVsNGT"){
    #radial kernel svm
    result_df <- svm_model(data.train, output_labels.train, data.test, output_labels.test, classes, kernel = "radial")
  }

  if(!dir.exists(result_file_dir)){
    dir.create(result_file_dir, recursive = TRUE)
  }
  write.csv(format(result_df, digits = 3), paste0(result_file_dir, result_file_name), row.names = FALSE)
  
}  



pipeline(
  comparison = "CFRDVsIGT",
  classes = c("IGT", "CFRD"),
  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
  dataset_replace_str = "CF_EV_AU_zlogtmm_",
  result_file_dir = "data/prediction_result/",
  result_file_name = "CFRDVsIGT.csv",
  perform_filter = TRUE,
  norm = "norm_log_tmm"
)


pipeline(
  comparison = "CFRDVsNGT",
  classes = c("NGT", "CFRD"),
  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
  dataset_replace_str = "CF_EV_AU_zlogtmm_",
  result_file_dir = "data/prediction_result/",
  result_file_name = "CFRDVsNGT.csv",
  perform_filter = TRUE,
  norm = "norm_log_tmm"
)


pipeline(
  comparison = "IGTVsNGT",
  classes = c("NGT", "IGT"),
  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
  dataset_replace_str = "CF_EV_AU_zlogtmm_",
  result_file_dir = "data/prediction_result/",
  result_file_name = "IGTVsNGT.csv",
  perform_filter = TRUE,
  norm = "norm_log_tmm"
)



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

  metrics <- data.frame(Comparison = comparison,
                        TrainAccuracy = acc.train,
                        TrainAUC = auc.train,
                        TestAccuracy = acc.test,
                        TestAUC = auc.test)
  
  write.table(x = metrics, file = metric_output_file_path, append = TRUE, 
              col.names = !file.exists(metric_output_file_path), sep = ",",
              row.names = FALSE)
    
}



#CFRDVsIGT

# comparison = "CFRDVsIGT"
# result_file_path <- "data/prediction_result/CFRDVsIGT.csv"


show_metrics(comparison = "CFRDVsIGT",
             classes = c("IGT", "CFRD"),
             result_file_path = "data/prediction_result/CFRDVsIGT.csv")


show_metrics(comparison = "CFRDVsNGT",
             classes = c("NGT", "CFRD"),
             result_file_path = "data/prediction_result/CFRDVsNGT.csv")


show_metrics(comparison = "IGTVsNGT",
             classes = c("NGT", "IGT"),
             result_file_path = "data/prediction_result/IGTVsNGT.csv")




#############################
#only adults with logtmm

pipeline(
  comparison = "CFRDVsIGT",
  classes = c("IGT", "CFRD"),
  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
  only_adults = TRUE,
  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
  result_file_dir = "data/prediction_result_adult/",
  result_file_name = "CFRDVsIGT.csv",
  perform_filter = TRUE,
  norm = "log_tmm"
)


pipeline(
  comparison = "CFRDVsNGT",
  classes = c("NGT", "CFRD"),
  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
  only_adults = TRUE,
  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
  result_file_dir = "data/prediction_result_adult/",
  result_file_name = "CFRDVsNGT.csv",
  perform_filter = TRUE,
  norm = "log_tmm"
)


pipeline(
  comparison = "IGTVsNGT",
  classes = c("NGT", "IGT"),
  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
  only_adults = TRUE,
  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
  result_file_dir = "data/prediction_result_adult/",
  result_file_name = "IGTVsNGT.csv",
  perform_filter = TRUE,
  norm = "log_tmm"
)

show_metrics(comparison = "CFRDVsIGT",
             classes = c("IGT", "CFRD"),
             result_file_path = "data/prediction_result_adult/CFRDVsIGT.csv",
             metric_output_file_path = "data/prediction_result_adult/metrics.csv")
show_metrics(comparison = "CFRDVsNGT",
             classes = c("NGT", "CFRD"),
             result_file_path = "data/prediction_result_adult/CFRDVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult/metrics.csv")
show_metrics(comparison = "IGTVsNGT",
             classes = c("NGT", "IGT"),
             result_file_path = "data/prediction_result_adult/IGTVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult/metrics.csv")

