#to combine predictions from multiple models into 1 using a specified classification model
library(tidyverse)

comparison = "CFRDVsNGT"
conditions = c("NGT", "CFRD")
data_file_path = "integration_prediction_result/CFRDVsNGT.csv"
o_type = "tra"
stack_model = "rf"
base_model_metrics_file_path = "integration_prediction_result/metrics_base_models.csv"
output_dir_path = "integration_prediction_result/stacked_top2/"

level2_stack <- function(comparison,
                         conditions,
                         data_file_path,
                         o_type,
                         base_model_metrics_file_path = NA,
                         output_dir_path = "integration_prediction_result/stacked/",
                         stack_model = "rf",
                         k_prod_times = 30,
                         prot_n = 2,
                         tra_n = 2){
  data <- read.csv(data_file_path) %>%
    dplyr::select(-c(PredictedLabel, cutoff))
  if(o_type != "both"){
    data <- data %>% 
      filter(omics_type == o_type)
  }
  if(!is.na(base_model_metrics_file_path)){
    #choose top 2 models only for stacking
    metrics <- read.csv(base_model_metrics_file_path) %>%
      dplyr::filter(Comparison == comparison, Type == "test") %>%
      dplyr::select(c(Omics_type, model, MeanAUC))
    metrics.prot <- metrics %>%
      filter(Omics_type == "prot") %>%
      slice_max(order_by = MeanAUC, n = prot_n)
    metrics.tra <- metrics %>%
      filter(Omics_type == "tra") %>%
      slice_max(order_by = MeanAUC, n = tra_n)
    
    data <- data %>%
      dplyr::filter((omics_type == "prot" & model %in% metrics.prot$model) | (omics_type == "tra" & model %in% metrics.tra$model))
  }
  
  data <- data %>%
    pivot_wider(names_from = c(omics_type, model), values_from = Pred_prob)
  data[is.na(data)] <- 0
  
  for(i in c(1:k_prod_times)){
    # i <- 1
    data_sub <- data %>%
      filter(iter == i) %>%
      dplyr::select(-c(iter))
    data.train <- data_sub %>%
      filter(Type == "train") %>%
      dplyr::select(-c(Type, TrueLabel)) %>%
      column_to_rownames("sample")
    label.train <- data_sub %>%
      filter(Type == "train") %>%
      dplyr::select(sample, TrueLabel) %>%
      rename("Label" = "TrueLabel")
    
    data.test <- data_sub %>%
      filter(Type == "test") %>%
      dplyr::select(-c(Type, TrueLabel)) %>%
      column_to_rownames("sample")
    label.test <- data_sub %>%
      filter(Type == "test") %>%
      dplyr::select(sample, TrueLabel) %>%
      rename("Label" = "TrueLabel")
    
    # result_df <- cbind(iter = i,
    #                    omics_type = paste(o_type, "stacked", sep = "_"),
    #                    log_reg_model(data.train, label.train, data.test, label.test, 
    #                                  classes = conditions, regularize = 'l2'))
    
    if(stack_model == "rf"){
      model_result <- rf_model(data.train, label.train, data.test, label.test,
                               classes = conditions)
    } else if(stack_model == "l2_log_reg"){
      model_result <- log_reg_model(data.train, label.train, data.test, label.test,
                                    classes = conditions, regularize = 'l2')
    } else if(stack_model == "l1_log_reg"){
      model_result <- log_reg_model(data.train, label.train, data.test, label.test,
                                    classes = conditions, regularize = 'l1')
    } else if(stack_model == "el_log_reg"){
      model_result <- log_reg_model(data.train, label.train, data.test, label.test,
                                    classes = conditions, regularize = 'el')
    } else if(stack_model == "sig_svm"){
      model_result <- svm_model(data.train, label.train, data.test, label.test,
                                    classes = conditions, kernel = "sigmoid")
    } else if(stack_model == "rad_svm"){
      model_result <- svm_model(data.train, label.train, data.test, label.test,
                                classes = conditions, kernel = "radial")
    }
    
    result_df <- cbind(iter = i,
                       omics_type = paste(o_type, "stacked", sep = "_"),
                       model_result)

    if(i == 1){
      result_df_all <- result_df
    } else{
      result_df_all <- rbind(result_df_all, result_df)
    }
  }
  result_file_dir_path <- paste0(output_dir_path, stack_model, "/")
  result_file_name <- paste0(comparison, "_", o_type,".csv")
  if(!dir.exists(result_file_dir_path)){
    dir.create(result_file_dir_path, recursive = TRUE)
  }
  write.csv(format(result_df_all, digits = 3), paste0(result_file_dir_path, result_file_name), 
            row.names = FALSE)
}


