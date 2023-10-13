#to combine predictions from multiple models into 1 using l2log_reg
library(tidyverse)

# comparison = "POSTOPE_TPVsREC_TP"
# conditions = c("REC_TP", "POSTOPE_TP")
# data_file_path = "Data/prediction_result/integration/POSTOPE_TPVsREC_TP.csv"
# o_type = "prot"
# o_type = "tra"
# o_type = "both"
# k_prod_times = 30

level2_stack <- function(comparison,
                         conditions,
                         data_file_path,
                         o_type,
                         stack_model = "rf",
                         k_prod_times = 30){
  data <- read.csv(data_file_path) %>%
    dplyr::select(-c(PredictedLabel, cutoff))
  if(o_type != "both"){
    data <- data %>% 
      filter(omics_type == o_type)
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
      result_df <- cbind(iter = i,
                         omics_type = paste(o_type, "stacked", sep = "_"),
                         rf_model(data.train, label.train, data.test, label.test,
                                  classes = conditions))      
    } else if(stack_model == "l2_log_reg"){
      result_df <- cbind(iter = i,
                         omics_type = paste(o_type, "stacked", sep = "_"),
                         log_reg_model(data.train, label.train, data.test, label.test,
                                  classes = conditions, regularize = 'l2'))   
    } else if(stack_model == "l1_log_reg"){
      result_df <- cbind(iter = i,
                         omics_type = paste(o_type, "stacked", sep = "_"),
                         log_reg_model(data.train, label.train, data.test, label.test,
                                       classes = conditions, regularize = 'l1'))   
    }

    
    
    if(i == 1){
      result_df_all <- result_df
    } else{
      result_df_all <- rbind(result_df_all, result_df)
    }
  }
  result_file_dir_path <- paste0("integration_prediction_result/stacked/", stack_model, "/")
  result_file_name <- paste0(comparison, "_", o_type,".csv")
  if(!dir.exists(result_file_dir_path)){
    dir.create(result_file_dir_path, recursive = TRUE)
  }
  write.csv(format(result_df_all, digits = 3), paste0(result_file_dir_path, result_file_name), 
            row.names = FALSE)
}


