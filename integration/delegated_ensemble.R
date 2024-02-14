library(tidyverse)


#options - 
# 1. prot otherwise prot + tra
# 2. prot otherwise tra

#prot stacked result
#alternate stacked result

# prot_stacked_result_file_path = "integration_prediction_result/stacked/rf/CFRDVsIGT_prot.csv" 
# alternate_stacked_result_file_path = "integration_prediction_result/stacked/rf/CFRDVsIGT_both.csv"
# alt_type = "both"
# comparison = "CFRDVsIGT"
# conditions = c("IGT", "CFRD")
# k_prod_times = 30

delegated_ensemble <- function(comparison,
                               conditions,
                               prot_stacked_result_file_path,
                               alternate_stacked_result_file_path,
                               alt_type,
                               k_prod_times = 30){
  data.prot <- read.csv(prot_stacked_result_file_path) %>%
    dplyr::select(-c(PredictedLabel, cutoff, model))
  data.alt <- read.csv(alternate_stacked_result_file_path) %>%
    dplyr::select(-c(PredictedLabel, cutoff, model)) %>%
    mutate(omics_type = "alt_stacked")
  
  data.combined <- rbind(data.prot, data.alt)
  
  data.combined <- data.combined %>%
    pivot_wider(names_from = omics_type, values_from = Pred_prob)
  sum(is.na(data.combined$prot_stacked))
  sum(is.na(data.combined$alt_stacked))
  
  for(i in c(1:k_prod_times)){
    # i <- 2
    data_sub.combined <- data.combined %>%
      filter(iter == i)
    
    data_sub.combined.train <- data_sub.combined %>%
      filter(Type == "train")
    
    best_cutoff <- 0.9
    best_auc <- 0
    for(cutoff in c(0.9, 0.8, 0.7)){
      data_sub.combined.train.pred_prob <- data_sub.combined.train %>%
        mutate(Pred_prob = case_when(is.na(prot_stacked) ~ alt_stacked,
                                     is.na(alt_stacked) ~ prot_stacked,
                                     ((prot_stacked > cutoff) | (prot_stacked < (1-cutoff))) ~ prot_stacked,
                                     TRUE ~ alt_stacked)) %>%
        mutate(source = case_when(is.na(prot_stacked) ~ "alt",
                                  is.na(alt_stacked) ~ "prot",
                                  ((prot_stacked > cutoff) | (prot_stacked < (1-cutoff))) ~ "prot",
                                  TRUE ~ "alt"))
      
      pr <- ROCR::prediction(data_sub.combined.train.pred_prob$Pred_prob, 
                             data_sub.combined.train.pred_prob$TrueLabel, label.ordering = conditions)
      auc <- ROCR::performance(pr, measure = "auc")@y.values[[1]]
      
      if(auc > best_auc){
        best_auc <- auc
        best_cutoff <- cutoff
      }
    }
    
    result_df <- data_sub.combined %>%
      mutate(Pred_prob = case_when(is.na(prot_stacked) ~ alt_stacked,
                                   is.na(alt_stacked) ~ prot_stacked,
                                   ((prot_stacked > best_cutoff) | (prot_stacked < (1-best_cutoff))) ~ prot_stacked,
                                   TRUE ~ alt_stacked)) %>%
      mutate(PredictedLabel = case_when(Pred_prob > 0.5 ~ conditions[2],
                                        TRUE ~ conditions[1])) %>%
      mutate(source = case_when(is.na(prot_stacked) ~ "alt",
                                is.na(alt_stacked) ~ "prot",
                                ((prot_stacked > best_cutoff) | (prot_stacked < (1-best_cutoff))) ~ "prot",
                                TRUE ~ "alt")) %>%
      mutate(delegate_cutoff = best_cutoff)
    
    if(i == 1){
      result_df_all <- result_df
    } else{
      result_df_all <- rbind(result_df_all, result_df)
    }
  }
  result_file_dir_path <- "integration_prediction_result/delegate/"
  result_file_name <- paste0(comparison, "_", alt_type,".csv")
  if(!dir.exists(result_file_dir_path)){
    dir.create(result_file_dir_path, recursive = TRUE)
  }
  write.csv(result_df_all, paste0(result_file_dir_path, result_file_name), 
            row.names = FALSE)
}


