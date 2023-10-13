# comparison = "POSTOPE_TPVsREC_TP"
# conditions = c("REC_TP", "POSTOPE_TP")
# o_type = "prot"
# m = "Simple logistic regression"
# sample_type = "test"
# result_file_path = "Data/prediction_result/integration/POSTOPE_TPVsREC_TP.csv"
# metric_output_file_path = "Data/prediction_result/integration/metrics_base_models.csv"
# 
# comparison = "POSTOPE_TPVsREC_TP"
# conditions = c("REC_TP", "POSTOPE_TP")
# o_type = "tra"
# m = "L2 Regularized logistic regression"
# sample_type = "train"
# result_file_path = paste0("Data/prediction_result/integration/stacked/POSTOPE_TPVsREC_TP",
#                           "_", "tra",
#                           ".csv")
# metric_output_file_path = "Data/prediction_result/integration/stacked/metrics.csv"


compute_metrics <- function(comparison,
                            conditions,
                            o_type, 
                            m, 
                            sample_type,
                            result_file_path, 
                            metric_output_file_path,
                            k_prod_times = 30){
  print(comparison)
  result_df <- read.csv(result_file_path) %>%
    filter(omics_type == o_type, Type == sample_type, model == m)
  acc_vec <- c()
  auc_vec <- c()
  tpr_vec <- c()
  tnr_vec <- c()
  for(i in c(1:k_prod_times)){
    # i <- 1
    result_df_sub <- result_df %>%
      filter(iter == i)
    
    acc <- mean(result_df_sub$TrueLabel == result_df_sub$PredictedLabel)
    
    pr <- ROCR::prediction(result_df_sub$Pred_prob, result_df_sub$TrueLabel, label.ordering = conditions)
    auc <- ROCR::performance(pr, measure = "auc")@y.values[[1]]
    
    tpr <- caret::sensitivity(factor(result_df_sub$PredictedLabel), 
                              factor(result_df_sub$TrueLabel), positive = conditions[2])
    tnr <- caret::specificity(factor(result_df_sub$PredictedLabel), 
                              factor(result_df_sub$TrueLabel), negative = conditions[1])
    
    acc_vec[i] <- acc 
    auc_vec[i] <- auc
    tpr_vec[i] <- tpr
    tnr_vec[i] <- tnr
  }
  
  metrics <- data.frame(Comparison = comparison,
                        Omics_type = o_type,
                        model = m,
                        Type = sample_type,
                        MeanAccuracy = mean(acc_vec),
                        MeanAUC = mean(auc_vec),
                        MeanTPR = mean(tpr_vec),
                        MeanTNR = mean(tnr_vec))
  
  write.table(x = metrics, file = metric_output_file_path, append = TRUE, 
              col.names = !file.exists(metric_output_file_path), sep = ",",
              row.names = FALSE)
}


# comparison = "PREOPEVsREC_TP"
# conditions = c("REC_TP", "PREOPE")
# alt_type = "tra"
# source_type = "prot"
# sample_type = "test"
# result_file_path = paste0("Data/prediction_result/integration/delegate/PREOPEVsREC_TP",
#                           "_", alt_type,
#                           ".csv")
# metric_output_file_path = "Data/prediction_result/integration/delegate/source_type_metrics.csv"
# k_prod_times = 30

compute_metrics.integrated <- function(comparison,
                                       conditions,
                                       alt_type, 
                                       sample_type,
                                       result_file_path, 
                                       metric_output_file_path,
                                       source_type = NA,
                                       k_prod_times = 30){
  print(comparison)
  result_df <- read.csv(result_file_path) %>%
    filter(Type == sample_type)
  acc_vec <- c()
  auc_vec <- c()
  tpr_vec <- c()
  tnr_vec <- c()
  source_ratio_vec <- c()
  for(i in c(1:k_prod_times)){
    # i <- 21
    result_df_sub <- result_df %>%
      filter(iter == i)
    if(!is.na(source_type)){
      result_df_sub2 <- result_df_sub %>%
        filter(source == source_type)
      source_ratio <- nrow(result_df_sub2) / nrow(result_df_sub)
      result_df_sub <- result_df_sub2
    }
    
    acc <- mean(result_df_sub$TrueLabel == result_df_sub$PredictedLabel)
    
    auc <- NA
    tpr <- NA
    tnr <- NA
    try({
      pr <- ROCR::prediction(result_df_sub$Pred_prob, result_df_sub$TrueLabel, label.ordering = conditions)
      auc <- ROCR::performance(pr, measure = "auc")@y.values[[1]]
      
      tpr <- caret::sensitivity(factor(result_df_sub$PredictedLabel), 
                                factor(result_df_sub$TrueLabel), positive = conditions[2])
      tnr <- caret::specificity(factor(result_df_sub$PredictedLabel), 
                                factor(result_df_sub$TrueLabel), negative = conditions[1])
    })
    
    acc_vec[i] <- acc 
    auc_vec[i] <- auc
    tpr_vec[i] <- tpr
    tnr_vec[i] <- tnr
    
    if(!is.na(source_type)){
      source_ratio_vec[i] <- source_ratio
    }
  }
  
  metrics <- data.frame(Comparison = comparison,
                        Alt_type = alt_type,
                        Type = sample_type,
                        MeanAccuracy = mean(acc_vec, na.rm = TRUE),
                        MeanAUC = mean(auc_vec, na.rm = TRUE),
                        MeanTPR = mean(tpr_vec, na.rm = TRUE),
                        MeanTNR = mean(tnr_vec, na.rm = TRUE))
  if(!is.na(source_type)){
    metrics <- metrics %>%
      mutate(source = source_type, .after = Type) %>%
      mutate(MeanSourceRatio = mean(source_ratio_vec), .after = source)
  }
  
  write.table(x = metrics, file = metric_output_file_path, append = TRUE, 
              col.names = !file.exists(metric_output_file_path), sep = ",",
              row.names = FALSE)
}
