library(tidyverse)

base_dir <- "~/UNSW/VafaeeLab/CysticFibrosisGroup/ExoCF/CFRD_EV_biomarker/"
setwd(base_dir)

source("prediction_pipeline/cm_logistic_regression.R")
source("prediction_pipeline/cm_svm.R")
source("prediction_pipeline/cm_rf.R")

comparison = "CFRDVsIGT"
classes = c("IGT", "CFRD")
result_file_dir = "data/prediction_result_phenotype/"

phenotype_pipeline <- function(comparison, classes,
                               result_file_dir = "data/prediction_result_phenotype/"){
  
  phenotype <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t")
  
  data <- phenotype %>%
    filter(is.na(pre_post_modulator) | pre_post_modulator == 0) %>%
    filter(condition %in% c("CFRD", "IGT", "NGT"), age_group == "adult") %>%
    select(Sample, condition, country, age, sex, FEV1) %>%
    mutate(age = as.numeric(str_trim(age)), FEV1 = as.numeric(str_trim(FEV1))) %>%
    mutate(is_M = ifelse(sex == "M", 1, 0),
           is_F = ifelse(sex == "F", 1, 0)) %>%
    column_to_rownames("Sample")
  data.train <- data %>%
    filter(country == "AU")
  data.test <- data %>%
    filter(country == "DK")
  
  #currently child samples have been filtered out, which results in no NAs for FEV1
  # #replace missing values in FEV1 with mean
  # data.train <- data.train %>%
  #   inner_join(data.train %>%
  #                group_by(condition) %>%
  #                summarize(mean_FEV1 = mean(FEV1, na.rm = TRUE)))
  # data.train <- data.train %>%
  #   mutate(FEV1 = case_when(is.na(FEV1) ~ mean_FEV1,
  #                           TRUE ~ FEV1)) %>%
  #   select(-c(mean_FEV1))
  # 
  
  #select conditions of interest to perform classification
  #doing this later so that replace FEV1 by mean of groups from whole data for train can be done if required
  #                                               i.e. for the case when child samples are not filtered out
  data.train <- data.train %>%
    filter(condition %in% classes)
  data.test <- data.test %>%
    filter(condition %in% classes) 
  
  x.train <- data.train %>%
    select(-c(condition, country, sex))
  x.test <- data.test %>%
    select(-c(condition, country, sex))
  
  y.train <- data.train %>%
    select(condition) %>%
    rename("Label" = condition)
  y.test <- data.test %>%
    select(condition) %>%
    rename("Label" = condition)
  
  #normalize 
  
  #uncomment this to not normalize one-hot encoded features i.e. col 3, 4)
  # x.train.norm_part <- x.train[, c(1:2)]
  # x.train.rest <- x.train[, c(3:4)]
  # 
  # x.test.norm_part <- x.test[, c(1:2)]
  # x.test.rest <- x.test[, c(3:4)]
  # 
  # normparam <- caret::preProcess(x.train.norm_part) 
  # x.train.norm_part <- predict(normparam, x.train.norm_part)
  # normparam <- caret::preProcess(x.test.norm_part) 
  # x.test.norm_part <- predict(normparam, x.test.norm_part) 
  # 
  # x.train <- cbind(x.train.norm_part, x.train.rest)
  # x.test <- cbind(x.test.norm_part, x.test.rest)
  
  normparam <- caret::preProcess(x.train) 
  x.train <- predict(normparam, x.train)
  x.test <- predict(normparam, x.test)
  
  result_df <- rf_model(x.train, y.train, x.test, y.test, classes)
  
  # radial kernel SVM doesn't give good results 
  # result_df <- svm_model(x.train, y.train, x.test, y.test, classes, kernel = "radial")
  
  if(!dir.exists(result_file_dir)){
    dir.create(result_file_dir, recursive = TRUE)
  }
  result_file_name <- paste0("phenotype_based_", comparison, ".csv")
  write.csv(format(result_df, digits = 3), paste0(result_file_dir, result_file_name), row.names = FALSE)
}  


phenotype_pipeline(comparison = "CFRDVsIGT",
                   classes = c("IGT", "CFRD"))
phenotype_pipeline(comparison = "CFRDVsNGT",
                   classes = c("NGT", "CFRD"))
phenotype_pipeline(comparison = "IGTVsNGT",
                   classes = c("NGT", "IGT"))


show_metrics <- function(comparison, classes, result_file_path){
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
  
  file_path = "data/prediction_result_phenotype/metrics.csv"
  write.table(x = metrics, file = file_path, append = TRUE, 
              col.names = !file.exists(file_path), sep = ",",
              row.names = FALSE)
  
}


show_metrics(comparison = "CFRDVsIGT",
             classes = c("IGT", "CFRD"),
             result_file_path = "data/prediction_result_phenotype/phenotype_based_CFRDVsIGT.csv")
show_metrics(comparison = "CFRDVsNGT",
             classes = c("NGT", "CFRD"),
             result_file_path = "data/prediction_result_phenotype/phenotype_based_CFRDVsNGT.csv")
show_metrics(comparison = "IGTVsNGT",
             classes = c("NGT", "IGT"),
             result_file_path = "data/prediction_result_phenotype/phenotype_based_IGTVsNGT.csv")

