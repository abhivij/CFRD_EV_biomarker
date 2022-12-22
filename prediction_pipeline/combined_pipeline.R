library(tidyverse)

base_dir <- "~/UNSW/VafaeeLab/CysticFibrosisGroup/ExoCF/CFRD_EV_biomarker/"
setwd(base_dir)


source("prediction_pipeline/cm_logistic_regression.R")
source("prediction_pipeline/cm_svm.R")
source("prediction_pipeline/cm_rf.R")

######################################################################


comparison = "CFRDVsNGT"
classes = c("NGT", "CFRD")
best_features_file_path  = "data/selected_features/best_features_with_is_best.csv"
only_adults = TRUE
train_cohort_country = "AU"
use_phenotype_info = FALSE
dataset_replace_str = "CF_EV_AU_adult_logcpm_"
result_file_dir = "data/prediction_result_adult_AU/"
result_file_name = "CFRDVsNGT.csv"
perform_filter = TRUE
norm = "log_cpm"
upsample = FALSE

combined_pipeline <- function(comparison, classes, 
                     best_features_file_path,
                     only_adults = FALSE,
                     train_cohort_country = "AU",
                     use_phenotype_info = TRUE,
                     dataset_replace_str = "CF_EV_AU_zlogtmm_",
                     result_file_dir = "data/prediction_result/",
                     result_file_name = "CFRDVsIGT.csv",
                     perform_filter = TRUE,
                     norm = "norm_log_tmm",
                     model = "rf",
                     upsample = FALSE){
  
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
  
  test_cohort_country <- ifelse(train_cohort_country == "AU", "DK", "AU") 
  
  train.output_labels <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes, country == train_cohort_country) %>%
    dplyr::select(Sample, Label, age, age_group, sex, FEV1) %>%
    dplyr::mutate(Label = factor(Label), age = as.numeric(age), 
                  age_group = factor(age_group), sex = factor(sex),
                  FEV1 = as.numeric(FEV1)) %>%
    arrange(Label, Sample)
  train.tra_data <- data[, train.output_labels$Sample]
  
  #if upsample is true, re-include (randomly) samples present in class with lower count 
  if(upsample == TRUE){
    #classes[2] assigned to class1 deliberately
    #classes are provided as c(neg_class, pos_class)
    class1 <- classes[2]
    class2 <- classes[1]
    class1_count <- train.output_labels %>% dplyr::filter(Label == class1) %>% nrow()
    class2_count <- train.output_labels %>% dplyr::filter(Label == class2) %>% nrow() 
    
    if(class1_count > class2_count){
      num_samples_to_add <- class1_count - class2_count
      class_to_add <- class2
    } else if(class1_count < class2_count){
      num_samples_to_add <- class2_count - class1_count 
      class_to_add <- class1
    }
    
    set.seed(1000)
    samples_in_class_to_add <- train.output_labels %>% dplyr::filter(Label == class_to_add)  
    
    #generate random numbers from 1 to number of samples in lower count class
    #generate num_samples_to_add number of random numbers
    samples_to_choose <- sample(c(1:nrow(samples_in_class_to_add)), num_samples_to_add)
    repeated_samples <- samples_in_class_to_add[samples_to_choose, ]  
    
    train.tra_data_repeat <- data[, repeated_samples$Sample, drop = FALSE]
    colnames(train.tra_data_repeat) <- paste0(colnames(train.tra_data_repeat), "_rep")
    train.tra_data <- cbind(train.tra_data, train.tra_data_repeat)
    
    repeated_samples <- repeated_samples %>%
      dplyr::mutate(Sample = paste0(Sample, "_rep"))
    train.output_labels <- rbind(train.output_labels,
                                 repeated_samples)
    rownames(train.output_labels) <- NULL
  }
  print(summary(train.output_labels))
  
  test.output_labels <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes, country == test_cohort_country) %>%
    dplyr::select(Sample, Label, age, age_group, sex, FEV1) %>%
    dplyr::mutate(Label = factor(Label), age = as.numeric(age), 
                  age_group = factor(age_group), sex = factor(sex),
                  FEV1 = as.numeric(FEV1)) %>%    
    arrange(Label, Sample)
  test.tra_data <- data[, test.output_labels$Sample]
  print(summary(test.output_labels))

  if(use_phenotype_info){
    train.pheno_data <- train.output_labels %>%
      dplyr::select(Sample, age, sex, FEV1) %>%
      mutate(age = as.numeric(str_trim(age)), FEV1 = as.numeric(str_trim(FEV1))) %>%
      mutate(is_M = ifelse(sex == "M", 1, 0),
             is_F = ifelse(sex == "F", 1, 0)) %>%
      dplyr::select(-c(sex)) %>%
      column_to_rownames("Sample")
    
    test.pheno_data <- test.output_labels %>%
      dplyr::select(Sample, age, sex, FEV1) %>%
      mutate(age = as.numeric(str_trim(age)), FEV1 = as.numeric(str_trim(FEV1))) %>%
      mutate(is_M = ifelse(sex == "M", 1, 0),
             is_F = ifelse(sex == "F", 1, 0)) %>%
      dplyr::select(-c(sex)) %>%
      column_to_rownames("Sample")
  }
  train.output_labels <- train.output_labels %>%
    dplyr::select(Sample, Label)
  test.output_labels <- test.output_labels %>%
    dplyr::select(Sample, Label)
  
  #currently train.tra_data, test.tra_data format : (transcripts x samples)
  if(perform_filter){
    keep <- edgeR::filterByExpr(train.tra_data, group = train.output_labels$Label)
    train.tra_data <- train.tra_data[keep, ]
    test.tra_data <- test.tra_data[keep, ]  
  }
  
  if(norm == "norm_log_tmm"){
    #calculating norm log tmm
    dge <- edgeR::DGEList(counts = train.tra_data, group = train.output_labels$Label)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    tmm <- edgeR::cpm(dge, log = TRUE)
    train.tra_data <- tmm
    
    dge <- edgeR::DGEList(counts = test.tra_data, group = test.output_labels$Label)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    tmm <- edgeR::cpm(dge, log = TRUE)
    test.tra_data <- tmm
    
    train.tra_data <- as.data.frame(t(as.matrix(train.tra_data)))
    test.tra_data <- as.data.frame(t(as.matrix(test.tra_data)))  
    
    #normalizing the data
    normparam <- caret::preProcess(train.tra_data) 
    train.tra_data <- predict(normparam, train.tra_data)
    test.tra_data <- predict(normparam, test.tra_data) #normalizing test data using params from train data   
    
  } else if(norm == "log_tmm"){
    #calculating norm log tmm
    dge <- edgeR::DGEList(counts = train.tra_data, group = train.output_labels$Label)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    tmm <- edgeR::cpm(dge, log = TRUE)
    train.tra_data <- tmm
    
    dge <- edgeR::DGEList(counts = test.tra_data, group = test.output_labels$Label)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    tmm <- edgeR::cpm(dge, log = TRUE)
    test.tra_data <- tmm
    
    train.tra_data <- as.data.frame(t(as.matrix(train.tra_data)))
    test.tra_data <- as.data.frame(t(as.matrix(test.tra_data)))  
  
  } else if(norm == "log_cpm"){
    
    train.tra_data <- edgeR::cpm(train.tra_data, log = TRUE)
    test.tra_data <- edgeR::cpm(test.tra_data, log = TRUE)
    train.tra_data <- as.data.frame(t(as.matrix(train.tra_data)))
    test.tra_data <- as.data.frame(t(as.matrix(test.tra_data)))  
    
  } else if(norm == "log"){
    train.tra_data[train.tra_data == 0] <- 2^-30
    train.tra_data <- log2(train.tra_data)
    
    test.tra_data[test.tra_data == 0] <- 2^-30
    test.tra_data <- log2(test.tra_data)
    
    train.tra_data <- as.data.frame(t(as.matrix(train.tra_data)))
    test.tra_data <- as.data.frame(t(as.matrix(test.tra_data)))  
  }
  
  #now train.tra_data, test.tra_data format : (samples x transcripts)
  
  #get best biomarkers only
  features_with_slash <- colnames(train.tra_data)[grepl("/", colnames(train.tra_data), fixed = TRUE)] 
  for(f in features_with_slash){
    f_replaced <- gsub("/|-", ".", f) 
    if(f_replaced %in% biomarkers){
      biomarkers[biomarkers == f_replaced] = f
    }
  }
  biomarkers <- gsub(".", "-", biomarkers, fixed = TRUE)
  
  train.tra_data <- train.tra_data[, biomarkers]
  test.tra_data <- test.tra_data[, biomarkers]
  
  if(use_phenotype_info){
    #normalize train.pheno_data, test.pheno_data
    normparam <- caret::preProcess(train.pheno_data) 
    train.pheno_data <- predict(normparam, train.pheno_data)
    test.pheno_data <- predict(normparam, test.pheno_data)
    
    #combine transcriptomic and phenotypic data
    train.data <- cbind(train.pheno_data, train.tra_data)
    test.data <- cbind(test.pheno_data, test.tra_data)    
  } else{
    train.data <- train.tra_data
    test.data <- test.tra_data
  }

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
  }
  
  # 


  
  if(!dir.exists(result_file_dir)){
    dir.create(result_file_dir, recursive = TRUE)
  }
  write.csv(format(result_df, digits = 3), paste0(result_file_dir, result_file_name), row.names = FALSE)
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

combined_pipeline(
  comparison = "CFRDVsIGT",
  classes = c("IGT", "CFRD"),
  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
  only_adults = TRUE,
  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
  result_file_dir = "data/prediction_result_tra_and_pheno_AUtrain/",
  result_file_name = "CFRDVsIGT.csv",
  perform_filter = TRUE,
  norm = "log_tmm"
)


combined_pipeline(
  comparison = "CFRDVsNGT",
  classes = c("NGT", "CFRD"),
  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
  only_adults = TRUE,
  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
  result_file_dir = "data/prediction_result_tra_and_pheno_AUtrain/",
  result_file_name = "CFRDVsNGT.csv",
  perform_filter = TRUE,
  norm = "log_tmm"
)


combined_pipeline(
  comparison = "IGTVsNGT",
  classes = c("NGT", "IGT"),
  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
  only_adults = TRUE,
  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
  result_file_dir = "data/prediction_result_tra_and_pheno_AUtrain/",
  result_file_name = "IGTVsNGT.csv",
  perform_filter = TRUE,
  norm = "log_tmm"
)

show_metrics(comparison = "CFRDVsIGT",
             classes = c("IGT", "CFRD"),
             result_file_path = "data/prediction_result_tra_and_pheno_AUtrain/CFRDVsIGT.csv",
             metric_output_file_path = "data/prediction_result_tra_and_pheno_AUtrain/metrics.csv")
show_metrics(comparison = "CFRDVsNGT",
             classes = c("NGT", "CFRD"),
             result_file_path = "data/prediction_result_tra_and_pheno_AUtrain/CFRDVsNGT.csv",
             metric_output_file_path = "data/prediction_result_tra_and_pheno_AUtrain/metrics.csv")
show_metrics(comparison = "IGTVsNGT",
             classes = c("NGT", "IGT"),
             result_file_path = "data/prediction_result_tra_and_pheno_AUtrain/IGTVsNGT.csv",
             metric_output_file_path = "data/prediction_result_tra_and_pheno_AUtrain/metrics.csv")




combined_pipeline(
  comparison = "CFRDVsIGT",
  classes = c("IGT", "CFRD"),
  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
  only_adults = TRUE,
  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
  result_file_dir = "data/prediction_result_tra_and_pheno_AUtrain/",
  result_file_name = "CFRDVsIGT_with_upsample.csv",
  perform_filter = TRUE,
  norm = "log_tmm",
  upsample = TRUE
)
combined_pipeline(
  comparison = "CFRDVsNGT",
  classes = c("NGT", "CFRD"),
  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
  only_adults = TRUE,
  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
  result_file_dir = "data/prediction_result_tra_and_pheno_AUtrain/",
  result_file_name = "CFRDVsNGT_with_upsample.csv",
  perform_filter = TRUE,
  norm = "log_tmm",
  upsample = TRUE
)
combined_pipeline(
  comparison = "IGTVsNGT",
  classes = c("NGT", "IGT"),
  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
  only_adults = TRUE,
  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
  result_file_dir = "data/prediction_result_tra_and_pheno_AUtrain/",
  result_file_name = "IGTVsNGT_with_upsample.csv",
  perform_filter = TRUE,
  norm = "log_tmm",
  upsample = TRUE
)

show_metrics(comparison = "CFRDVsIGT",
             classes = c("IGT", "CFRD"),
             result_file_path = "data/prediction_result_tra_and_pheno_AUtrain/CFRDVsIGT_with_upsample.csv",
             metric_output_file_path = "data/prediction_result_tra_and_pheno_AUtrain/metrics.csv")
show_metrics(comparison = "CFRDVsNGT",
             classes = c("NGT", "CFRD"),
             result_file_path = "data/prediction_result_tra_and_pheno_AUtrain/CFRDVsNGT_with_upsample.csv",
             metric_output_file_path = "data/prediction_result_tra_and_pheno_AUtrain/metrics.csv")
show_metrics(comparison = "IGTVsNGT",
             classes = c("NGT", "IGT"),
             result_file_path = "data/prediction_result_tra_and_pheno_AUtrain/IGTVsNGT_with_upsample.csv",
             metric_output_file_path = "data/prediction_result_tra_and_pheno_AUtrain/metrics.csv")




####################




combined_pipeline(
  comparison = "CFRDVsIGT",
  classes = c("IGT", "CFRD"),
  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
  only_adults = TRUE,
  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
  result_file_dir = "data/prediction_result_tra_and_pheno_AUtrain/",
  result_file_name = "CFRDVsIGT_radsvm.csv",
  perform_filter = TRUE,
  norm = "log_tmm"
)


combined_pipeline(
  comparison = "CFRDVsNGT",
  classes = c("NGT", "CFRD"),
  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
  only_adults = TRUE,
  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
  result_file_dir = "data/prediction_result_tra_and_pheno_AUtrain/",
  result_file_name = "CFRDVsNGT_radsvm.csv",
  perform_filter = TRUE,
  norm = "log_tmm"
)


combined_pipeline(
  comparison = "IGTVsNGT",
  classes = c("NGT", "IGT"),
  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
  only_adults = TRUE,
  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
  result_file_dir = "data/prediction_result_tra_and_pheno_AUtrain/",
  result_file_name = "IGTVsNGT_radsvm.csv",
  perform_filter = TRUE,
  norm = "log_tmm"
)

show_metrics(comparison = "CFRDVsIGT",
             classes = c("IGT", "CFRD"),
             result_file_path = "data/prediction_result_tra_and_pheno_AUtrain/CFRDVsIGT_radsvm.csv",
             metric_output_file_path = "data/prediction_result_tra_and_pheno_AUtrain/metrics_radsvm.csv")
show_metrics(comparison = "CFRDVsNGT",
             classes = c("NGT", "CFRD"),
             result_file_path = "data/prediction_result_tra_and_pheno_AUtrain/CFRDVsNGT_radsvm.csv",
             metric_output_file_path = "data/prediction_result_tra_and_pheno_AUtrain/metrics_radsvm.csv")
show_metrics(comparison = "IGTVsNGT",
             classes = c("NGT", "IGT"),
             result_file_path = "data/prediction_result_tra_and_pheno_AUtrain/IGTVsNGT_radsvm.csv",
             metric_output_file_path = "data/prediction_result_tra_and_pheno_AUtrain/metrics_radsvm.csv")




combined_pipeline(
  comparison = "CFRDVsIGT",
  classes = c("IGT", "CFRD"),
  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
  only_adults = TRUE,
  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
  result_file_dir = "data/prediction_result_tra_and_pheno_AUtrain/",
  result_file_name = "CFRDVsIGT_with_upsample_radsvm.csv",
  perform_filter = TRUE,
  norm = "log_tmm",
  upsample = TRUE
)
combined_pipeline(
  comparison = "CFRDVsNGT",
  classes = c("NGT", "CFRD"),
  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
  only_adults = TRUE,
  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
  result_file_dir = "data/prediction_result_tra_and_pheno_AUtrain/",
  result_file_name = "CFRDVsNGT_with_upsample_radsvm.csv",
  perform_filter = TRUE,
  norm = "log_tmm",
  upsample = TRUE
)
combined_pipeline(
  comparison = "IGTVsNGT",
  classes = c("NGT", "IGT"),
  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
  only_adults = TRUE,
  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
  result_file_dir = "data/prediction_result_tra_and_pheno_AUtrain/",
  result_file_name = "IGTVsNGT_with_upsample_radsvm.csv",
  perform_filter = TRUE,
  norm = "log_tmm",
  upsample = TRUE
)

show_metrics(comparison = "CFRDVsIGT",
             classes = c("IGT", "CFRD"),
             result_file_path = "data/prediction_result_tra_and_pheno_AUtrain/CFRDVsIGT_with_upsample_radsvm.csv",
             metric_output_file_path = "data/prediction_result_tra_and_pheno_AUtrain/metrics_radsvm.csv")
show_metrics(comparison = "CFRDVsNGT",
             classes = c("NGT", "CFRD"),
             result_file_path = "data/prediction_result_tra_and_pheno_AUtrain/CFRDVsNGT_with_upsample_radsvm.csv",
             metric_output_file_path = "data/prediction_result_tra_and_pheno_AUtrain/metrics_radsvm.csv")
show_metrics(comparison = "IGTVsNGT",
             classes = c("NGT", "IGT"),
             result_file_path = "data/prediction_result_tra_and_pheno_AUtrain/IGTVsNGT_with_upsample_radsvm.csv",
             metric_output_file_path = "data/prediction_result_tra_and_pheno_AUtrain/metrics_radsvm.csv")


##############
#train on DK, test on AU

combined_pipeline(comparison = "CFRDVsIGT", classes = c("IGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "DK",
                  use_phenotype_info = FALSE,
                  dataset_replace_str = "CF_EV_DK_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_DK/",
                  result_file_name = "CFRDVsIGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  upsample = FALSE)
combined_pipeline(comparison = "CFRDVsNGT", classes = c("NGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "DK",
                  use_phenotype_info = FALSE,
                  dataset_replace_str = "CF_EV_DK_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_DK/",
                  result_file_name = "CFRDVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  upsample = FALSE)
combined_pipeline(comparison = "IGTVsNGT", classes = c("NGT", "IGT"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "DK",
                  use_phenotype_info = FALSE,
                  dataset_replace_str = "CF_EV_DK_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_DK/",
                  result_file_name = "IGTVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  upsample = FALSE)

show_metrics(comparison = "CFRDVsIGT",
             classes = c("IGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_DK/CFRDVsIGT.csv",
             metric_output_file_path = "data/prediction_result_adult_DK/metrics.csv")
show_metrics(comparison = "CFRDVsNGT",
             classes = c("NGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_DK/CFRDVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_DK/metrics.csv")
show_metrics(comparison = "IGTVsNGT",
             classes = c("NGT", "IGT"),
             result_file_path = "data/prediction_result_adult_DK/IGTVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_DK/metrics.csv")


###############

combined_pipeline(comparison = "CFRDVsIGT", classes = c("IGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "DK",
                  use_phenotype_info = FALSE,
                  dataset_replace_str = "CF_EV_DK_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_DK_upsample/",
                  result_file_name = "CFRDVsIGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  upsample = TRUE)
combined_pipeline(comparison = "CFRDVsNGT", classes = c("NGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "DK",
                  use_phenotype_info = FALSE,
                  dataset_replace_str = "CF_EV_DK_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_DK_upsample/",
                  result_file_name = "CFRDVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  upsample = TRUE)
combined_pipeline(comparison = "IGTVsNGT", classes = c("NGT", "IGT"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "DK",
                  use_phenotype_info = FALSE,
                  dataset_replace_str = "CF_EV_DK_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_DK_upsample/",
                  result_file_name = "IGTVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  upsample = TRUE)

show_metrics(comparison = "CFRDVsIGT",
             classes = c("IGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_DK_upsample/CFRDVsIGT.csv",
             metric_output_file_path = "data/prediction_result_adult_DK_upsample/metrics.csv")
show_metrics(comparison = "CFRDVsNGT",
             classes = c("NGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_DK_upsample/CFRDVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_DK_upsample/metrics.csv")
show_metrics(comparison = "IGTVsNGT",
             classes = c("NGT", "IGT"),
             result_file_path = "data/prediction_result_adult_DK_upsample/IGTVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_DK_upsample/metrics.csv")


#############################

combined_pipeline(comparison = "CFRDVsIGT", classes = c("IGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "DK",
                  use_phenotype_info = FALSE,
                  dataset_replace_str = "CF_EV_DK_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_DK_l2logreg/",
                  result_file_name = "CFRDVsIGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm", model = "l2_log_reg",
                  upsample = FALSE)
combined_pipeline(comparison = "CFRDVsNGT", classes = c("NGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "DK",
                  use_phenotype_info = FALSE,
                  dataset_replace_str = "CF_EV_DK_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_DK_l2logreg/",
                  result_file_name = "CFRDVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm", model = "l2_log_reg",
                  upsample = FALSE)
combined_pipeline(comparison = "IGTVsNGT", classes = c("NGT", "IGT"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "DK",
                  use_phenotype_info = FALSE,
                  dataset_replace_str = "CF_EV_DK_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_DK_l2logreg/",
                  result_file_name = "IGTVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm", model = "l2_log_reg",
                  upsample = FALSE)

show_metrics(comparison = "CFRDVsIGT",
             classes = c("IGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_DK_l2logreg/CFRDVsIGT.csv",
             metric_output_file_path = "data/prediction_result_adult_DK_l2logreg/metrics.csv")
show_metrics(comparison = "CFRDVsNGT",
             classes = c("NGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_DK_l2logreg/CFRDVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_DK_l2logreg/metrics.csv")
show_metrics(comparison = "IGTVsNGT",
             classes = c("NGT", "IGT"),
             result_file_path = "data/prediction_result_adult_DK_l2logreg/IGTVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_DK_l2logreg/metrics.csv")



#############################

combined_pipeline(comparison = "CFRDVsIGT", classes = c("IGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "DK",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_DK_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_DK_with_pheno/",
                  result_file_name = "CFRDVsIGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  upsample = FALSE)
combined_pipeline(comparison = "CFRDVsNGT", classes = c("NGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "DK",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_DK_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_DK_with_pheno/",
                  result_file_name = "CFRDVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  upsample = FALSE)
combined_pipeline(comparison = "IGTVsNGT", classes = c("NGT", "IGT"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "DK",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_DK_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_DK_with_pheno/",
                  result_file_name = "IGTVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  upsample = FALSE)

show_metrics(comparison = "CFRDVsIGT",
             classes = c("IGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_DK_with_pheno/CFRDVsIGT.csv",
             metric_output_file_path = "data/prediction_result_adult_DK_with_pheno/metrics.csv")
show_metrics(comparison = "CFRDVsNGT",
             classes = c("NGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_DK_with_pheno/CFRDVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_DK_with_pheno/metrics.csv")
show_metrics(comparison = "IGTVsNGT",
             classes = c("NGT", "IGT"),
             result_file_path = "data/prediction_result_adult_DK_with_pheno/IGTVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_DK_with_pheno/metrics.csv")

#################


combined_pipeline(comparison = "CFRDVsIGT", classes = c("IGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "DK",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_DK_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_DK_with_pheno_with_std_scale/",
                  result_file_name = "CFRDVsIGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  upsample = FALSE)
combined_pipeline(comparison = "CFRDVsNGT", classes = c("NGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "DK",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_DK_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_DK_with_pheno_with_std_scale/",
                  result_file_name = "CFRDVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  upsample = FALSE)
combined_pipeline(comparison = "IGTVsNGT", classes = c("NGT", "IGT"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "DK",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_DK_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_DK_with_pheno_with_std_scale/",
                  result_file_name = "IGTVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  upsample = FALSE)

show_metrics(comparison = "CFRDVsIGT",
             classes = c("IGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_DK_with_pheno_with_std_scale/CFRDVsIGT.csv",
             metric_output_file_path = "data/prediction_result_adult_DK_with_pheno_with_std_scale/metrics.csv")
show_metrics(comparison = "CFRDVsNGT",
             classes = c("NGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_DK_with_pheno_with_std_scale/CFRDVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_DK_with_pheno_with_std_scale/metrics.csv")
show_metrics(comparison = "IGTVsNGT",
             classes = c("NGT", "IGT"),
             result_file_path = "data/prediction_result_adult_DK_with_pheno_with_std_scale/IGTVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_DK_with_pheno_with_std_scale/metrics.csv")


combined_pipeline(comparison = "CFRDVsIGT", classes = c("IGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "DK",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_DK_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_DK_l2logreg_with_pheno_std_scale/",
                  result_file_name = "CFRDVsIGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm", model = "l2_log_reg",
                  upsample = FALSE)
combined_pipeline(comparison = "CFRDVsNGT", classes = c("NGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "DK",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_DK_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_DK_l2logreg_with_pheno_std_scale/",
                  result_file_name = "CFRDVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm", model = "l2_log_reg",
                  upsample = FALSE)
combined_pipeline(comparison = "IGTVsNGT", classes = c("NGT", "IGT"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "DK",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_DK_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_DK_l2logreg_with_pheno_std_scale/",
                  result_file_name = "IGTVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm", model = "l2_log_reg",
                  upsample = FALSE)

show_metrics(comparison = "CFRDVsIGT",
             classes = c("IGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_DK_l2logreg_with_pheno_std_scale/CFRDVsIGT.csv",
             metric_output_file_path = "data/prediction_result_adult_DK_l2logreg_with_pheno_std_scale/metrics.csv")
show_metrics(comparison = "CFRDVsNGT",
             classes = c("NGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_DK_l2logreg_with_pheno_std_scale/CFRDVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_DK_l2logreg_with_pheno_std_scale/metrics.csv")
show_metrics(comparison = "IGTVsNGT",
             classes = c("NGT", "IGT"),
             result_file_path = "data/prediction_result_adult_DK_l2logreg_with_pheno_std_scale/IGTVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_DK_l2logreg_with_pheno_std_scale/metrics.csv")

#####################

combined_pipeline(comparison = "CFRDVsIGT", classes = c("IGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "DK",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_DK_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_DK_l2logreg_with_pheno/",
                  result_file_name = "CFRDVsIGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm", model = "l2_log_reg",
                  upsample = FALSE)
combined_pipeline(comparison = "CFRDVsNGT", classes = c("NGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "DK",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_DK_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_DK_l2logreg_with_pheno/",
                  result_file_name = "CFRDVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm", model = "l2_log_reg",
                  upsample = FALSE)
combined_pipeline(comparison = "IGTVsNGT", classes = c("NGT", "IGT"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "DK",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_DK_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_DK_l2logreg_with_pheno/",
                  result_file_name = "IGTVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm", model = "l2_log_reg",
                  upsample = FALSE)

show_metrics(comparison = "CFRDVsIGT",
             classes = c("IGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_DK_l2logreg_with_pheno/CFRDVsIGT.csv",
             metric_output_file_path = "data/prediction_result_adult_DK_l2logreg_with_pheno/metrics.csv")
show_metrics(comparison = "CFRDVsNGT",
             classes = c("NGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_DK_l2logreg_with_pheno/CFRDVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_DK_l2logreg_with_pheno/metrics.csv")
show_metrics(comparison = "IGTVsNGT",
             classes = c("NGT", "IGT"),
             result_file_path = "data/prediction_result_adult_DK_l2logreg_with_pheno/IGTVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_DK_l2logreg_with_pheno/metrics.csv")


#####################

#rerun AU since RF model cutoffs were changed

combined_pipeline(comparison = "CFRDVsIGT", classes = c("IGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = FALSE,
                  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_AU/",
                  result_file_name = "CFRDVsIGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  upsample = FALSE)
combined_pipeline(comparison = "CFRDVsNGT", classes = c("NGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = FALSE,
                  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_AU/",
                  result_file_name = "CFRDVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  upsample = FALSE)
combined_pipeline(comparison = "IGTVsNGT", classes = c("NGT", "IGT"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = FALSE,
                  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_AU/",
                  result_file_name = "IGTVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  upsample = FALSE)

show_metrics(comparison = "CFRDVsIGT",
             classes = c("IGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_AU/CFRDVsIGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU/metrics.csv")
show_metrics(comparison = "CFRDVsNGT",
             classes = c("NGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_AU/CFRDVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU/metrics.csv")
show_metrics(comparison = "IGTVsNGT",
             classes = c("NGT", "IGT"),
             result_file_path = "data/prediction_result_adult_AU/IGTVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU/metrics.csv")

#########################################

combined_pipeline(comparison = "CFRDVsIGT", classes = c("IGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = FALSE,
                  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_AU_with_upsample/",
                  result_file_name = "CFRDVsIGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  upsample = TRUE)
combined_pipeline(comparison = "CFRDVsNGT", classes = c("NGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = FALSE,
                  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_AU_with_upsample/",
                  result_file_name = "CFRDVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  upsample = TRUE)
combined_pipeline(comparison = "IGTVsNGT", classes = c("NGT", "IGT"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = FALSE,
                  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_AU_with_upsample/",
                  result_file_name = "IGTVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  upsample = TRUE)

show_metrics(comparison = "CFRDVsIGT",
             classes = c("IGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_AU_with_upsample/CFRDVsIGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_with_upsample/metrics.csv")
show_metrics(comparison = "CFRDVsNGT",
             classes = c("NGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_AU_with_upsample/CFRDVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_with_upsample/metrics.csv")
show_metrics(comparison = "IGTVsNGT",
             classes = c("NGT", "IGT"),
             result_file_path = "data/prediction_result_adult_AU_with_upsample/IGTVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_with_upsample/metrics.csv")


##############

combined_pipeline(comparison = "CFRDVsIGT", classes = c("IGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_AU_with_pheno/",
                  result_file_name = "CFRDVsIGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  upsample = FALSE)
combined_pipeline(comparison = "CFRDVsNGT", classes = c("NGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_AU_with_pheno/",
                  result_file_name = "CFRDVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  upsample = FALSE)
combined_pipeline(comparison = "IGTVsNGT", classes = c("NGT", "IGT"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_AU_with_pheno/",
                  result_file_name = "IGTVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  upsample = FALSE)

show_metrics(comparison = "CFRDVsIGT",
             classes = c("IGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_AU_with_pheno/CFRDVsIGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_with_pheno/metrics.csv")
show_metrics(comparison = "CFRDVsNGT",
             classes = c("NGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_AU_with_pheno/CFRDVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_with_pheno/metrics.csv")
show_metrics(comparison = "IGTVsNGT",
             classes = c("NGT", "IGT"),
             result_file_path = "data/prediction_result_adult_AU_with_pheno/IGTVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_with_pheno/metrics.csv")

##################

combined_pipeline(comparison = "CFRDVsIGT", classes = c("IGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_AU_with_pheno_upsample/",
                  result_file_name = "CFRDVsIGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  upsample = TRUE)
combined_pipeline(comparison = "CFRDVsNGT", classes = c("NGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_AU_with_pheno_upsample/",
                  result_file_name = "CFRDVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  upsample = TRUE)
combined_pipeline(comparison = "IGTVsNGT", classes = c("NGT", "IGT"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_AU_with_pheno_upsample/",
                  result_file_name = "IGTVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  upsample = TRUE)

show_metrics(comparison = "CFRDVsIGT",
             classes = c("IGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_AU_with_pheno_upsample/CFRDVsIGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_with_pheno_upsample/metrics.csv")
show_metrics(comparison = "CFRDVsNGT",
             classes = c("NGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_AU_with_pheno_upsample/CFRDVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_with_pheno_upsample/metrics.csv")
show_metrics(comparison = "IGTVsNGT",
             classes = c("NGT", "IGT"),
             result_file_path = "data/prediction_result_adult_AU_with_pheno_upsample/IGTVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_with_pheno_upsample/metrics.csv")


##############

combined_pipeline(comparison = "CFRDVsIGT", classes = c("IGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = FALSE,
                  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_AU_radsvm/",
                  result_file_name = "CFRDVsIGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  model = "radial_svm",
                  upsample = FALSE)
combined_pipeline(comparison = "CFRDVsNGT", classes = c("NGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = FALSE,
                  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_AU_radsvm/",
                  result_file_name = "CFRDVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  model = "radial_svm",
                  upsample = FALSE)
combined_pipeline(comparison = "IGTVsNGT", classes = c("NGT", "IGT"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = FALSE,
                  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_AU_radsvm/",
                  result_file_name = "IGTVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  model = "radial_svm",
                  upsample = FALSE)

show_metrics(comparison = "CFRDVsIGT",
             classes = c("IGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_AU_radsvm/CFRDVsIGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_radsvm/metrics.csv")
show_metrics(comparison = "CFRDVsNGT",
             classes = c("NGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_AU_radsvm/CFRDVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_radsvm/metrics.csv")
show_metrics(comparison = "IGTVsNGT",
             classes = c("NGT", "IGT"),
             result_file_path = "data/prediction_result_adult_AU_radsvm/IGTVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_radsvm/metrics.csv")

##############

combined_pipeline(comparison = "CFRDVsIGT", classes = c("IGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = FALSE,
                  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_AU_radsvm_with_upsample/",
                  result_file_name = "CFRDVsIGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  model = "radial_svm",
                  upsample = TRUE)
combined_pipeline(comparison = "CFRDVsNGT", classes = c("NGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = FALSE,
                  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_AU_radsvm_with_upsample/",
                  result_file_name = "CFRDVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  model = "radial_svm",
                  upsample = TRUE)
combined_pipeline(comparison = "IGTVsNGT", classes = c("NGT", "IGT"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = FALSE,
                  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_AU_radsvm_with_upsample/",
                  result_file_name = "IGTVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  model = "radial_svm",
                  upsample = TRUE)

show_metrics(comparison = "CFRDVsIGT",
             classes = c("IGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_AU_radsvm_with_upsample/CFRDVsIGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_radsvm_with_upsample/metrics.csv")
show_metrics(comparison = "CFRDVsNGT",
             classes = c("NGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_AU_radsvm_with_upsample/CFRDVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_radsvm_with_upsample/metrics.csv")
show_metrics(comparison = "IGTVsNGT",
             classes = c("NGT", "IGT"),
             result_file_path = "data/prediction_result_adult_AU_radsvm_with_upsample/IGTVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_radsvm_with_upsample/metrics.csv")

##############

combined_pipeline(comparison = "CFRDVsIGT", classes = c("IGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_AU_radsvm_with_pheno/",
                  result_file_name = "CFRDVsIGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  model = "radial_svm",
                  upsample = FALSE)
combined_pipeline(comparison = "CFRDVsNGT", classes = c("NGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_AU_radsvm_with_pheno/",
                  result_file_name = "CFRDVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  model = "radial_svm",
                  upsample = FALSE)
combined_pipeline(comparison = "IGTVsNGT", classes = c("NGT", "IGT"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_AU_radsvm_with_pheno/",
                  result_file_name = "IGTVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  model = "radial_svm",
                  upsample = FALSE)

show_metrics(comparison = "CFRDVsIGT",
             classes = c("IGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_AU_radsvm_with_pheno/CFRDVsIGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_radsvm_with_pheno/metrics.csv")
show_metrics(comparison = "CFRDVsNGT",
             classes = c("NGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_AU_radsvm_with_pheno/CFRDVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_radsvm_with_pheno/metrics.csv")
show_metrics(comparison = "IGTVsNGT",
             classes = c("NGT", "IGT"),
             result_file_path = "data/prediction_result_adult_AU_radsvm_with_pheno/IGTVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_radsvm_with_pheno/metrics.csv")

##############

combined_pipeline(comparison = "CFRDVsIGT", classes = c("IGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_AU_radsvm_with_pheno_with_upsample/",
                  result_file_name = "CFRDVsIGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  model = "radial_svm",
                  upsample = TRUE)
combined_pipeline(comparison = "CFRDVsNGT", classes = c("NGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_AU_radsvm_with_pheno_with_upsample/",
                  result_file_name = "CFRDVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  model = "radial_svm",
                  upsample = TRUE)
combined_pipeline(comparison = "IGTVsNGT", classes = c("NGT", "IGT"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
                  result_file_dir = "data/prediction_result_adult_AU_radsvm_with_pheno_with_upsample/",
                  result_file_name = "IGTVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_tmm",
                  model = "radial_svm",
                  upsample = TRUE)

show_metrics(comparison = "CFRDVsIGT",
             classes = c("IGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_AU_radsvm_with_pheno_with_upsample/CFRDVsIGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_radsvm_with_pheno_with_upsample/metrics.csv")
show_metrics(comparison = "CFRDVsNGT",
             classes = c("NGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_AU_radsvm_with_pheno_with_upsample/CFRDVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_radsvm_with_pheno_with_upsample/metrics.csv")
show_metrics(comparison = "IGTVsNGT",
             classes = c("NGT", "IGT"),
             result_file_path = "data/prediction_result_adult_AU_radsvm_with_pheno_with_upsample/IGTVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_radsvm_with_pheno_with_upsample/metrics.csv")

##############

combined_pipeline(comparison = "CFRDVsIGT", classes = c("IGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = FALSE,
                  dataset_replace_str = "CF_EV_AU_adult_logcpm_",
                  result_file_dir = "data/prediction_result_adult_AU_adult_logcpm/",
                  result_file_name = "CFRDVsIGT.csv",
                  perform_filter = TRUE,
                  norm = "log_cpm",
                  model = "rf",
                  upsample = FALSE)
combined_pipeline(comparison = "CFRDVsNGT", classes = c("NGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = FALSE,
                  dataset_replace_str = "CF_EV_AU_adult_logcpm_",
                  result_file_dir = "data/prediction_result_adult_AU_adult_logcpm/",
                  result_file_name = "CFRDVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_cpm",
                  model = "rf",
                  upsample = FALSE)
combined_pipeline(comparison = "IGTVsNGT", classes = c("NGT", "IGT"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = FALSE,
                  dataset_replace_str = "CF_EV_AU_adult_logcpm_",
                  result_file_dir = "data/prediction_result_adult_AU_adult_logcpm/",
                  result_file_name = "IGTVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_cpm",
                  model = "rf",
                  upsample = FALSE)

show_metrics(comparison = "CFRDVsIGT",
             classes = c("IGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_AU_adult_logcpm/CFRDVsIGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_adult_logcpm/metrics.csv")
show_metrics(comparison = "CFRDVsNGT",
             classes = c("NGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_AU_adult_logcpm/CFRDVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_adult_logcpm/metrics.csv")
show_metrics(comparison = "IGTVsNGT",
             classes = c("NGT", "IGT"),
             result_file_path = "data/prediction_result_adult_AU_adult_logcpm/IGTVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_adult_logcpm/metrics.csv")

##############

combined_pipeline(comparison = "CFRDVsIGT", classes = c("IGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_AU_adult_logcpm_",
                  result_file_dir = "data/prediction_result_adult_AU_adult_logcpm_with_pheno/",
                  result_file_name = "CFRDVsIGT.csv",
                  perform_filter = TRUE,
                  norm = "log_cpm",
                  model = "rf",
                  upsample = FALSE)
combined_pipeline(comparison = "CFRDVsNGT", classes = c("NGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_AU_adult_logcpm_",
                  result_file_dir = "data/prediction_result_adult_AU_adult_logcpm_with_pheno/",
                  result_file_name = "CFRDVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_cpm",
                  model = "rf",
                  upsample = FALSE)
combined_pipeline(comparison = "IGTVsNGT", classes = c("NGT", "IGT"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_AU_adult_logcpm_",
                  result_file_dir = "data/prediction_result_adult_AU_adult_logcpm_with_pheno/",
                  result_file_name = "IGTVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log_cpm",
                  model = "rf",
                  upsample = FALSE)

show_metrics(comparison = "CFRDVsIGT",
             classes = c("IGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_AU_adult_logcpm_with_pheno/CFRDVsIGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_adult_logcpm_with_pheno/metrics.csv")
show_metrics(comparison = "CFRDVsNGT",
             classes = c("NGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_AU_adult_logcpm_with_pheno/CFRDVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_adult_logcpm_with_pheno/metrics.csv")
show_metrics(comparison = "IGTVsNGT",
             classes = c("NGT", "IGT"),
             result_file_path = "data/prediction_result_adult_AU_adult_logcpm_with_pheno/IGTVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_adult_logcpm_with_pheno/metrics.csv")

##############

combined_pipeline(comparison = "CFRDVsIGT", classes = c("IGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = FALSE,
                  dataset_replace_str = "CF_EV_AU_adult_log_",
                  result_file_dir = "data/prediction_result_adult_AU_adult_log/",
                  result_file_name = "CFRDVsIGT.csv",
                  perform_filter = TRUE,
                  norm = "log",
                  model = "rf",
                  upsample = FALSE)
combined_pipeline(comparison = "CFRDVsNGT", classes = c("NGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = FALSE,
                  dataset_replace_str = "CF_EV_AU_adult_log_",
                  result_file_dir = "data/prediction_result_adult_AU_adult_log/",
                  result_file_name = "CFRDVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log",
                  model = "rf",
                  upsample = FALSE)
combined_pipeline(comparison = "IGTVsNGT", classes = c("NGT", "IGT"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = FALSE,
                  dataset_replace_str = "CF_EV_AU_adult_log_",
                  result_file_dir = "data/prediction_result_adult_AU_adult_log/",
                  result_file_name = "IGTVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log",
                  model = "rf",
                  upsample = FALSE)

show_metrics(comparison = "CFRDVsIGT",
             classes = c("IGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_AU_adult_log/CFRDVsIGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_adult_log/metrics.csv")
show_metrics(comparison = "CFRDVsNGT",
             classes = c("NGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_AU_adult_log/CFRDVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_adult_log/metrics.csv")
show_metrics(comparison = "IGTVsNGT",
             classes = c("NGT", "IGT"),
             result_file_path = "data/prediction_result_adult_AU_adult_log/IGTVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_adult_log/metrics.csv")


#######################


combined_pipeline(comparison = "CFRDVsIGT", classes = c("IGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_AU_adult_log_",
                  result_file_dir = "data/prediction_result_adult_AU_adult_log_with_pheno/",
                  result_file_name = "CFRDVsIGT.csv",
                  perform_filter = TRUE,
                  norm = "log",
                  model = "rf",
                  upsample = FALSE)
combined_pipeline(comparison = "CFRDVsNGT", classes = c("NGT", "CFRD"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_AU_adult_log_",
                  result_file_dir = "data/prediction_result_adult_AU_adult_log_with_pheno/",
                  result_file_name = "CFRDVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log",
                  model = "rf",
                  upsample = FALSE)
combined_pipeline(comparison = "IGTVsNGT", classes = c("NGT", "IGT"), 
                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv",
                  only_adults = TRUE,
                  train_cohort_country = "AU",
                  use_phenotype_info = TRUE,
                  dataset_replace_str = "CF_EV_AU_adult_log_",
                  result_file_dir = "data/prediction_result_adult_AU_adult_log_with_pheno/",
                  result_file_name = "IGTVsNGT.csv",
                  perform_filter = TRUE,
                  norm = "log",
                  model = "rf",
                  upsample = FALSE)

show_metrics(comparison = "CFRDVsIGT",
             classes = c("IGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_AU_adult_log_with_pheno/CFRDVsIGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_adult_log_with_pheno/metrics.csv")
show_metrics(comparison = "CFRDVsNGT",
             classes = c("NGT", "CFRD"),
             result_file_path = "data/prediction_result_adult_AU_adult_log_with_pheno/CFRDVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_adult_log_with_pheno/metrics.csv")
show_metrics(comparison = "IGTVsNGT",
             classes = c("NGT", "IGT"),
             result_file_path = "data/prediction_result_adult_AU_adult_log_with_pheno/IGTVsNGT.csv",
             metric_output_file_path = "data/prediction_result_adult_AU_adult_log_with_pheno/metrics.csv")

