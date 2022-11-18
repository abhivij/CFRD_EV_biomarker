library(tidyverse)

base_dir <- "~/UNSW/VafaeeLab/CysticFibrosisGroup/ExoCF/CFRD_EV_biomarker/"
setwd(base_dir)


source("prediction_pipeline/cm_logistic_regression.R")
source("prediction_pipeline/cm_svm.R")
source("prediction_pipeline/cm_rf.R")

######################################################################

# 
# comparison = "CFRDVsIGT"
# classes = c("IGT", "CFRD")
# best_features_file_path = "data/selected_features/best_features_with_is_best.csv"
# only_adults = TRUE
# dataset_replace_str = "CF_EV_AU_adult_logtmm_"
# result_file_dir = "data/prediction_result_adult/"
# result_file_name = "CFRDVsIGT.csv"
# perform_filter = TRUE
# norm = "log_tmm"

#using RF - since RF performs well in transcriptomics and exceptionally well on phenotypic data
combined_pipeline <- function(comparison, classes, 
                     best_features_file_path,
                     only_adults = FALSE,
                     dataset_replace_str = "CF_EV_AU_zlogtmm_",
                     result_file_dir = "data/prediction_result/",
                     result_file_name = "CFRDVsIGT.csv",
                     perform_filter = TRUE,
                     norm = "norm_log_tmm",
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
  
  train.output_labels <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes, country == "AU") %>%
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
    
    train.tra_data_repeat <- data[, repeated_samples$Sample]
    colnames(train.tra_data_repeat) <- paste0(colnames(train.tra_data_repeat), "_rep")
    train.tra_data <- cbind(train.tra_data, train.tra_data_repeat)
    
    repeated_samples <- repeated_samples %>%
      dplyr::mutate(Sample = paste0(Sample, "_rep"))
    train.output_labels <- rbind(train.output_labels,
                                 repeated_samples)
    rownames(train.output_labels) <- NULL
  }
  
  print(summary(train.output_labels))
  train.pheno_data <- train.output_labels %>%
    dplyr::select(Sample, age, sex, FEV1) %>%
    mutate(age = as.numeric(str_trim(age)), FEV1 = as.numeric(str_trim(FEV1))) %>%
    mutate(is_M = ifelse(sex == "M", 1, 0),
           is_F = ifelse(sex == "F", 1, 0)) %>%
    dplyr::select(-c(sex)) %>%
    column_to_rownames("Sample")
  train.output_labels <- train.output_labels %>%
    dplyr::select(Sample, Label)
  
  test.output_labels <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes, country == "DK") %>%
    dplyr::select(Sample, Label, age, age_group, sex, FEV1) %>%
    dplyr::mutate(Label = factor(Label), age = as.numeric(age), 
                  age_group = factor(age_group), sex = factor(sex),
                  FEV1 = as.numeric(FEV1)) %>%    
    arrange(Label, Sample)
  test.tra_data <- data[, test.output_labels$Sample]
  print(summary(test.output_labels))
  test.pheno_data <- test.output_labels %>%
    dplyr::select(Sample, age, sex, FEV1) %>%
    mutate(age = as.numeric(str_trim(age)), FEV1 = as.numeric(str_trim(FEV1))) %>%
    mutate(is_M = ifelse(sex == "M", 1, 0),
           is_F = ifelse(sex == "F", 1, 0)) %>%
    dplyr::select(-c(sex)) %>%
    column_to_rownames("Sample")
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
  
  #normalize train.pheno_data, test.pheno_data
  normparam <- caret::preProcess(train.pheno_data) 
  train.pheno_data <- predict(normparam, train.pheno_data)
  test.pheno_data <- predict(normparam, test.pheno_data)
  
  #combine transcriptomic and phenotypic data
  train.data <- cbind(train.pheno_data, train.tra_data)
  test.data <- cbind(test.pheno_data, test.tra_data)
  
  # result_df <- rf_model(train.data, train.output_labels, test.data, test.output_labels, classes)
  result_df <- svm_model(train.data, train.output_labels, test.data, test.output_labels, classes, kernel = "radial")
  
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
