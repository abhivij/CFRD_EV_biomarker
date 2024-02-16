base_dir <- "/home/abhivij/UNSW/VafaeeLab/CysticFibrosisGroup/ExoCF/CFRD_EV_biomarker/"
setwd(base_dir)

source("integration/base_model_pipeline.R")

base_model_pipeline(comparison = "CFRDVsIGT",
                    conditions = c("IGT", "CFRD"),
                    data_file_path_prot <- "data/formatted/proteomics/CFRDVsIGT_imputed333_mf_quantile_combat.csv",
                    data_file_path_tra <- "data/formatted/CFRDVsIGT_umi_counts_combat_processed.csv")

base_model_pipeline(comparison = "CFRDVsNGT",
                    conditions = c("NGT", "CFRD"),
                    data_file_path_prot <- "data/formatted/proteomics/CFRDVsNGT_imputed333_mf_quantile_combat.csv",
                    data_file_path_tra <- "data/formatted/CFRDVsNGT_umi_counts_combat_processed.csv")

base_model_pipeline(comparison = "IGTVsNGT",
                    conditions = c("NGT", "IGT"),
                    data_file_path_prot <- "data/formatted/proteomics/IGTVsNGT_imputed333_mf_quantile_combat.csv",
                    data_file_path_tra <- "data/formatted/IGTVsNGT_umi_counts_combat_processed.csv")

# <start> with phenotypes and associated data from Jan 2024

base_model_pipeline(comparison = "CFRDVsIGT",
                    conditions = c("IGT", "CFRD"),
                    data_file_path_prot = "data/formatted/proteomics/CFRDVsIGT_imputed333_mf_quantile_combat.csv",
                    data_file_path_tra = "data/formatted/rna_all/umi_counts_filter90.csv",
                    prot_bec = "combat",
                    tra_bec = "none",
                    dataset_replace_string_prot = "CF_EV_prot_mf_quantile_combat_",
                    dataset_replace_string_tra = "CF_EV_tra_334_")
base_model_pipeline(comparison = "CFRDVsNGT",
                    conditions = c("NGT", "CFRD"),
                    data_file_path_prot = "data/formatted/proteomics/CFRDVsNGT_imputed333_mf_quantile_combat.csv",
                    data_file_path_tra = "data/formatted/rna_all/umi_counts_filter90.csv",
                    prot_bec = "combat",
                    tra_bec = "none",
                    dataset_replace_string_prot = "CF_EV_prot_mf_quantile_combat_",
                    dataset_replace_string_tra = "CF_EV_tra_334_")
base_model_pipeline(comparison = "IGTVsNGT",
                    conditions = c("NGT", "IGT"),
                    data_file_path_prot = "data/formatted/proteomics/IGTVsNGT_imputed333_mf_quantile_combat.csv",
                    data_file_path_tra = "data/formatted/rna_all/umi_counts_filter90.csv",
                    prot_bec = "combat",
                    tra_bec = "none",
                    dataset_replace_string_prot = "CF_EV_prot_mf_quantile_combat_",
                    dataset_replace_string_tra = "CF_EV_tra_334_")

# <end> with phenotypes and associated data from Jan 2024

# rest of the executions are same as previous
# integration_prediction_result has been renamed to integration_prediction_result_old

models <- read.csv("integration_prediction_result/CFRDVsIGT.csv") %>%
  dplyr::select(model) %>%
  distinct()
models <- models$model
omics_types <- c("prot", "tra")
sample_types <- c("train", "test") 

for(omics_type in omics_types){
  for(model in models){
    for(sample_type in sample_types){
      compute_metrics(comparison = "CFRDVsIGT",
                      conditions = c("IGT", "CFRD"),
                      o_type = omics_type,
                      m = model,
                      sample_type = sample_type,
                      result_file_path = "integration_prediction_result/CFRDVsIGT.csv",
                      metric_output_file_path = "integration_prediction_result/metrics_base_models.csv")
      compute_metrics(comparison = "CFRDVsNGT",
                      conditions = c("NGT", "CFRD"),
                      o_type = omics_type,
                      m = model,
                      sample_type = sample_type,
                      result_file_path = "integration_prediction_result/CFRDVsNGT.csv",
                      metric_output_file_path = "integration_prediction_result/metrics_base_models.csv")
      compute_metrics(comparison = "IGTVsNGT",
                      conditions = c("NGT", "IGT"),
                      o_type = omics_type,
                      m = model,
                      sample_type = sample_type,
                      result_file_path = "integration_prediction_result/IGTVsNGT.csv",
                      metric_output_file_path = "integration_prediction_result/metrics_base_models.csv")
    }
  }  
}

metrics <- read.csv("integration_prediction_result/metrics_base_models.csv") %>%
  arrange(Comparison, Omics_type, Type, model)
write.csv(metrics, "integration_prediction_result/metrics_base_models_arranged.csv", row.names = FALSE)



#level 2 stack

# for(o in c("prot", "tra", "both")){
#   level2_stack(comparison = "CFRDVsIGT",
#                conditions = c("IGT", "CFRD"),
#                data_file_path = "integration_prediction_result/CFRDVsIGT.csv",
#                o_type = o)
#   level2_stack(comparison = "CFRDVsNGT",
#                conditions = c("NGT", "CFRD"),
#                data_file_path = "integration_prediction_result/CFRDVsNGT.csv",
#                o_type = o)
#   level2_stack(comparison = "IGTVsNGT",
#                conditions = c("NGT", "IGT"),
#                data_file_path = "integration_prediction_result/IGTVsNGT.csv",
#                o_type = o)
# }
# 
# for(o in c("tra", "prot", "both")){
#   for(sample_type in sample_types){
#     print(paste(o, sample_type))
#     compute_metrics(comparison = "CFRDVsIGT",
#                     conditions = c("IGT", "CFRD"),
#                     o_type = paste(o, "stacked", sep = "_"),
#                     m = "Random Forest",
#                     sample_type = sample_type,
#                     result_file_path = paste0("integration_prediction_result/stacked/rf/CFRDVsIGT",
#                                               "_", o,
#                                               ".csv"),
#                     metric_output_file_path = "integration_prediction_result/stacked/rf/metrics.csv")
#     compute_metrics(comparison = "CFRDVsNGT",
#                     conditions = c("NGT", "CFRD"),
#                     o_type = paste(o, "stacked", sep = "_"),
#                     m = "Random Forest",
#                     sample_type = sample_type,
#                     result_file_path = paste0("integration_prediction_result/stacked/rf/CFRDVsNGT",
#                                               "_", o,
#                                               ".csv"),
#                     metric_output_file_path = "integration_prediction_result/stacked/rf/metrics.csv")
#     compute_metrics(comparison = "IGTVsNGT",
#                     conditions = c("NGT", "IGT"),
#                     o_type = paste(o, "stacked", sep = "_"),
#                     m = "Random Forest",
#                     sample_type = sample_type,
#                     result_file_path = paste0("integration_prediction_result/stacked/rf/IGTVsNGT",
#                                               "_", o,
#                                               ".csv"),
#                     metric_output_file_path = "integration_prediction_result/stacked/rf/metrics.csv")
#     print(paste(o, sample_type, "done------------------"))
#   }
# }


classification_models <- c(
  "rf" = "Random Forest",
  "l1_log_reg" = "L1 Regularized logistic regression",
  "l2_log_reg" = "L2 Regularized logistic regression",
  "el_log_reg" = "Elastic net logistic regression",
  "sig_svm" = "Sigmoid Kernel SVM",
  "rad_svm" = "Radial Kernel SVM"
)

for(s_model in names(classification_models)){
  print(s_model)
  print(classification_models[s_model])
}


for(s_model in names(classification_models)){
  
  print(s_model)
  
  for(o in c("prot", "tra", "both")){
    level2_stack(comparison = "CFRDVsIGT",
                 conditions = c("IGT", "CFRD"),
                 data_file_path = "integration_prediction_result/CFRDVsIGT.csv",
                 o_type = o,
                 stack_model = s_model)
    level2_stack(comparison = "CFRDVsNGT",
                 conditions = c("NGT", "CFRD"),
                 data_file_path = "integration_prediction_result/CFRDVsNGT.csv",
                 o_type = o,
                 stack_model = s_model)
    level2_stack(comparison = "IGTVsNGT",
                 conditions = c("NGT", "IGT"),
                 data_file_path = "integration_prediction_result/IGTVsNGT.csv",
                 o_type = o,
                 stack_model = s_model)
  }  

  for(o in c("tra", "prot", "both")){
    for(sample_type in sample_types){
      print(paste(o, sample_type))
      compute_metrics(comparison = "CFRDVsIGT",
                      conditions = c("IGT", "CFRD"),
                      o_type = paste(o, "stacked", sep = "_"),
                      m = classification_models[s_model],
                      sample_type = sample_type,
                      result_file_path = paste0("integration_prediction_result/stacked/", s_model, "/CFRDVsIGT",
                                                "_", o,
                                                ".csv"),
                      metric_output_file_path = paste0("integration_prediction_result/stacked/", s_model, "/metrics.csv"))
      compute_metrics(comparison = "CFRDVsNGT",
                      conditions = c("NGT", "CFRD"),
                      o_type = paste(o, "stacked", sep = "_"),
                      m = classification_models[s_model],
                      sample_type = sample_type,
                      result_file_path = paste0("integration_prediction_result/stacked/", s_model, "/CFRDVsNGT",
                                                "_", o,
                                                ".csv"),
                      metric_output_file_path = paste0("integration_prediction_result/stacked/", s_model, "/metrics.csv"))
      compute_metrics(comparison = "IGTVsNGT",
                      conditions = c("NGT", "IGT"),
                      o_type = paste(o, "stacked", sep = "_"),
                      m = classification_models[s_model],
                      sample_type = sample_type,
                      result_file_path = paste0("integration_prediction_result/stacked/", s_model, "/IGTVsNGT",
                                                "_", o,
                                                ".csv"),
                      metric_output_file_path = paste0("integration_prediction_result/stacked/", s_model, "/metrics.csv"))
      print(paste(o, sample_type, "done------------------"))
    }
  }

}





# 
# 
# for(o in c("tra", "prot", "both")){
#   for(sample_type in sample_types){
#     print(paste(o, sample_type))
#     compute_metrics(comparison = "CFRDVsIGT",
#                     conditions = c("IGT", "CFRD"),
#                     o_type = paste(o, "stacked", sep = "_"),
#                     m = "L2 Regularized logistic regression",
#                     sample_type = sample_type,
#                     result_file_path = paste0("integration_prediction_result/stacked/l2_log_reg/CFRDVsIGT",
#                                               "_", o,
#                                               ".csv"),
#                     metric_output_file_path = "integration_prediction_result/stacked/l2_log_reg/metrics.csv")
#     compute_metrics(comparison = "CFRDVsNGT",
#                     conditions = c("NGT", "CFRD"),
#                     o_type = paste(o, "stacked", sep = "_"),
#                     m = "L2 Regularized logistic regression",
#                     sample_type = sample_type,
#                     result_file_path = paste0("integration_prediction_result/stacked/l2_log_reg/CFRDVsNGT",
#                                               "_", o,
#                                               ".csv"),
#                     metric_output_file_path = "integration_prediction_result/stacked/l2_log_reg/metrics.csv")
#     compute_metrics(comparison = "IGTVsNGT",
#                     conditions = c("NGT", "IGT"),
#                     o_type = paste(o, "stacked", sep = "_"),
#                     m = "L2 Regularized logistic regression",
#                     sample_type = sample_type,
#                     result_file_path = paste0("integration_prediction_result/stacked/l2_log_reg/IGTVsNGT",
#                                               "_", o,
#                                               ".csv"),
#                     metric_output_file_path = "integration_prediction_result/stacked/l2_log_reg/metrics.csv")
#     print(paste(o, sample_type, "done------------------"))
#   }
# }
# 
# for(o in c("tra", "prot", "both")){
#   for(sample_type in sample_types){
#     print(paste(o, sample_type))
#     compute_metrics(comparison = "CFRDVsIGT",
#                     conditions = c("IGT", "CFRD"),
#                     o_type = paste(o, "stacked", sep = "_"),
#                     m = "L1 Regularized logistic regression",
#                     sample_type = sample_type,
#                     result_file_path = paste0("integration_prediction_result/stacked/l1_log_reg/CFRDVsIGT",
#                                               "_", o,
#                                               ".csv"),
#                     metric_output_file_path = "integration_prediction_result/stacked/l1_log_reg/metrics.csv")
#     compute_metrics(comparison = "CFRDVsNGT",
#                     conditions = c("NGT", "CFRD"),
#                     o_type = paste(o, "stacked", sep = "_"),
#                     m = "L1 Regularized logistic regression",
#                     sample_type = sample_type,
#                     result_file_path = paste0("integration_prediction_result/stacked/l1_log_reg/CFRDVsNGT",
#                                               "_", o,
#                                               ".csv"),
#                     metric_output_file_path = "integration_prediction_result/stacked/l1_log_reg/metrics.csv")
#     compute_metrics(comparison = "IGTVsNGT",
#                     conditions = c("NGT", "IGT"),
#                     o_type = paste(o, "stacked", sep = "_"),
#                     m = "L1 Regularized logistic regression",
#                     sample_type = sample_type,
#                     result_file_path = paste0("integration_prediction_result/stacked/l1_log_reg/IGTVsNGT",
#                                               "_", o,
#                                               ".csv"),
#                     metric_output_file_path = "integration_prediction_result/stacked/l1_log_reg/metrics.csv")
#     print(paste(o, sample_type, "done------------------"))
#   }
# }




# 6: In svm.default(data.train, factor(label.train$Label, levels = classes),  ... :
#                     Variable(s) ‘prot_L1.Regularized.logistic.regression’ constant. Cannot scale data.
# 7: In svm.default(data.train, factor(label.train$Label, levels = classes),  ... :
#                     Variable(s) ‘prot_L1.Regularized.logistic.regression’ and ‘prot_Elastic.net.logistic.regression’ constant. Cannot scale data.
# 8: In svm.default(data.train, factor(label.train$Label, levels = classes),  ... :
#                     Variable(s) ‘prot_L1.Regularized.logistic.regression’ and ‘prot_Elastic.net.logistic.regression’ constant. Cannot scale data.



delegated_ensemble(comparison = "CFRDVsIGT",
                   conditions = c("IGT", "CFRD"),
                   prot_stacked_result_file_path = "integration_prediction_result/stacked/rf/CFRDVsIGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result/stacked/rf/CFRDVsIGT_tra.csv",
                   alt_type = "tra")
delegated_ensemble(comparison = "CFRDVsIGT",
                   conditions = c("IGT", "CFRD"),
                   prot_stacked_result_file_path = "integration_prediction_result/stacked/rf/CFRDVsIGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result/stacked/rf/CFRDVsIGT_both.csv",
                   alt_type = "both")
delegated_ensemble(comparison = "CFRDVsNGT",
                   conditions = c("NGT", "CFRD"),
                   prot_stacked_result_file_path = "integration_prediction_result/stacked/rf/CFRDVsNGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result/stacked/rf/CFRDVsNGT_tra.csv",
                   alt_type = "tra")
delegated_ensemble(comparison = "CFRDVsNGT",
                   conditions = c("NGT", "CFRD"),
                   prot_stacked_result_file_path = "integration_prediction_result/stacked/rf/CFRDVsNGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result/stacked/rf/CFRDVsNGT_both.csv",
                   alt_type = "both")
delegated_ensemble(comparison = "IGTVsNGT",
                   conditions = c("NGT", "IGT"),
                   prot_stacked_result_file_path = "integration_prediction_result/stacked/rf/IGTVsNGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result/stacked/rf/IGTVsNGT_tra.csv",
                   alt_type = "tra")
delegated_ensemble(comparison = "IGTVsNGT",
                   conditions = c("NGT", "IGT"),
                   prot_stacked_result_file_path = "integration_prediction_result/stacked/rf/IGTVsNGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result/stacked/rf/IGTVsNGT_both.csv",
                   alt_type = "both")

for(sample_type in c("train", "test")){
  for(alt_type in c("tra", "both")){
    compute_metrics.integrated(comparison = "CFRDVsIGT",
                               conditions = c("IGT", "CFRD"),
                               alt_type = alt_type,
                               sample_type = sample_type,
                               result_file_path = paste0("integration_prediction_result/delegate/CFRDVsIGT",
                                                         "_", alt_type,
                                                         ".csv"),
                               metric_output_file_path = "integration_prediction_result/delegate/metrics.csv")
    compute_metrics.integrated(comparison = "CFRDVsNGT",
                               conditions = c("NGT", "CFRD"),
                               alt_type = alt_type,
                               sample_type = sample_type,
                               result_file_path = paste0("integration_prediction_result/delegate/CFRDVsNGT",
                                                         "_", alt_type,
                                                         ".csv"),
                               metric_output_file_path = "integration_prediction_result/delegate/metrics.csv")
    compute_metrics.integrated(comparison = "IGTVsNGT",
                               conditions = c("NGT", "IGT"),
                               alt_type = alt_type,
                               sample_type = sample_type,
                               result_file_path = paste0("integration_prediction_result/delegate/IGTVsNGT",
                                                         "_", alt_type,
                                                         ".csv"),
                               metric_output_file_path = "integration_prediction_result/delegate/metrics.csv")
  }
}

##############

#stack using top 2 models

classification_models <- c(
  "rf" = "Random Forest",
  "l1_log_reg" = "L1 Regularized logistic regression",
  "l2_log_reg" = "L2 Regularized logistic regression",
  "el_log_reg" = "Elastic net logistic regression",
  "sig_svm" = "Sigmoid Kernel SVM",
  "rad_svm" = "Radial Kernel SVM"
)

for(s_model in names(classification_models)){
  print(s_model)
  print(classification_models[s_model])
}

models <- read.csv("integration_prediction_result/CFRDVsIGT.csv") %>%
  dplyr::select(model) %>%
  distinct()
models <- models$model
omics_types <- c("prot", "tra")
sample_types <- c("train", "test") 



for(s_model in names(classification_models)){
  
  print(s_model)
  
  for(o in c("prot", "tra", "both")){
    print(o)
    print("CFRDVsIGT")
    level2_stack(comparison = "CFRDVsIGT",
                 conditions = c("IGT", "CFRD"),
                 data_file_path = "integration_prediction_result/CFRDVsIGT.csv",
                 o_type = o,
                 stack_model = s_model,
                 base_model_metrics_file_path = "integration_prediction_result/metrics_base_models.csv",
                 output_dir_path = "integration_prediction_result/stacked_top2/")
    print("CFRDVsNGT")
    if(o == "tra"){
      level2_stack(comparison = "CFRDVsNGT",
                   conditions = c("NGT", "CFRD"),
                   data_file_path = "integration_prediction_result/CFRDVsNGT.csv",
                   o_type = o,
                   stack_model = s_model,
                   base_model_metrics_file_path = "integration_prediction_result/metrics_base_models.csv",
                   output_dir_path = "integration_prediction_result/stacked_top2/", 
                   tra_n = 3)      #using top 3 models for tra - top 2 have similar predictions and so cause issues
    } else{
      level2_stack(comparison = "CFRDVsNGT",
                   conditions = c("NGT", "CFRD"),
                   data_file_path = "integration_prediction_result/CFRDVsNGT.csv",
                   o_type = o,
                   stack_model = s_model,
                   base_model_metrics_file_path = "integration_prediction_result/metrics_base_models.csv",
                   output_dir_path = "integration_prediction_result/stacked_top2/")
    }
    print("IGTVsNGT")
    level2_stack(comparison = "IGTVsNGT",
                 conditions = c("NGT", "IGT"),
                 data_file_path = "integration_prediction_result/IGTVsNGT.csv",
                 o_type = o,
                 stack_model = s_model,
                 base_model_metrics_file_path = "integration_prediction_result/metrics_base_models.csv",
                 output_dir_path = "integration_prediction_result/stacked_top2/")
    print(paste(o, "done"))
  }  
  
  print("computing metrics")
  for(o in c("tra", "prot", "both")){
    for(sample_type in sample_types){
      print(paste(o, sample_type))
      compute_metrics(comparison = "CFRDVsIGT",
                      conditions = c("IGT", "CFRD"),
                      o_type = paste(o, "stacked", sep = "_"),
                      m = classification_models[s_model],
                      sample_type = sample_type,
                      result_file_path = paste0("integration_prediction_result/stacked_top2/", s_model, "/CFRDVsIGT",
                                                "_", o,
                                                ".csv"),
                      metric_output_file_path = paste0("integration_prediction_result/stacked_top2/", s_model, "/metrics.csv"))
      compute_metrics(comparison = "CFRDVsNGT",
                      conditions = c("NGT", "CFRD"),
                      o_type = paste(o, "stacked", sep = "_"),
                      m = classification_models[s_model],
                      sample_type = sample_type,
                      result_file_path = paste0("integration_prediction_result/stacked_top2/", s_model, "/CFRDVsNGT",
                                                "_", o,
                                                ".csv"),
                      metric_output_file_path = paste0("integration_prediction_result/stacked_top2/", s_model, "/metrics.csv"))
      compute_metrics(comparison = "IGTVsNGT",
                      conditions = c("NGT", "IGT"),
                      o_type = paste(o, "stacked", sep = "_"),
                      m = classification_models[s_model],
                      sample_type = sample_type,
                      result_file_path = paste0("integration_prediction_result/stacked_top2/", s_model, "/IGTVsNGT",
                                                "_", o,
                                                ".csv"),
                      metric_output_file_path = paste0("integration_prediction_result/stacked_top2/", s_model, "/metrics.csv"))
      print(paste(o, sample_type, "done------------------"))
    }
  }
  
}


delegated_ensemble(comparison = "CFRDVsIGT",
                   conditions = c("IGT", "CFRD"),
                   prot_stacked_result_file_path = "integration_prediction_result/stacked_top2/l2_log_reg/CFRDVsIGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result/stacked_top2/l2_log_reg/CFRDVsIGT_tra.csv",
                   alt_type = "tra",
                   result_file_dir_path = "integration_prediction_result/delegate_top2_l2_log_reg/")
delegated_ensemble(comparison = "CFRDVsIGT",
                   conditions = c("IGT", "CFRD"),
                   prot_stacked_result_file_path = "integration_prediction_result/stacked_top2/l2_log_reg/CFRDVsIGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result/stacked_top2/l2_log_reg/CFRDVsIGT_both.csv",
                   alt_type = "both",
                   result_file_dir_path = "integration_prediction_result/delegate_top2_l2_log_reg/")
delegated_ensemble(comparison = "CFRDVsNGT",
                   conditions = c("NGT", "CFRD"),
                   prot_stacked_result_file_path = "integration_prediction_result/stacked_top2/l2_log_reg/CFRDVsNGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result/stacked_top2/l2_log_reg/CFRDVsNGT_tra.csv",
                   alt_type = "tra",
                   result_file_dir_path = "integration_prediction_result/delegate_top2_l2_log_reg/")
delegated_ensemble(comparison = "CFRDVsNGT",
                   conditions = c("NGT", "CFRD"),
                   prot_stacked_result_file_path = "integration_prediction_result/stacked_top2/l2_log_reg/CFRDVsNGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result/stacked_top2/l2_log_reg/CFRDVsNGT_both.csv",
                   alt_type = "both",
                   result_file_dir_path = "integration_prediction_result/delegate_top2_l2_log_reg/")
delegated_ensemble(comparison = "IGTVsNGT",
                   conditions = c("NGT", "IGT"),
                   prot_stacked_result_file_path = "integration_prediction_result/stacked_top2/l2_log_reg/IGTVsNGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result/stacked_top2/l2_log_reg/IGTVsNGT_tra.csv",
                   alt_type = "tra",
                   result_file_dir_path = "integration_prediction_result/delegate_top2_l2_log_reg/")
delegated_ensemble(comparison = "IGTVsNGT",
                   conditions = c("NGT", "IGT"),
                   prot_stacked_result_file_path = "integration_prediction_result/stacked_top2/l2_log_reg/IGTVsNGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result/stacked_top2/l2_log_reg/IGTVsNGT_both.csv",
                   alt_type = "both",
                   result_file_dir_path = "integration_prediction_result/delegate_top2_l2_log_reg/")

for(sample_type in c("train", "test")){
  for(alt_type in c("tra", "both")){
    compute_metrics.integrated(comparison = "CFRDVsIGT",
                               conditions = c("IGT", "CFRD"),
                               alt_type = alt_type,
                               sample_type = sample_type,
                               result_file_path = paste0("integration_prediction_result/delegate_top2_l2_log_reg/CFRDVsIGT",
                                                         "_", alt_type,
                                                         ".csv"),
                               metric_output_file_path = "integration_prediction_result/delegate_top2_l2_log_reg/metrics.csv")
    compute_metrics.integrated(comparison = "CFRDVsNGT",
                               conditions = c("NGT", "CFRD"),
                               alt_type = alt_type,
                               sample_type = sample_type,
                               result_file_path = paste0("integration_prediction_result/delegate_top2_l2_log_reg/CFRDVsNGT",
                                                         "_", alt_type,
                                                         ".csv"),
                               metric_output_file_path = "integration_prediction_result/delegate_top2_l2_log_reg/metrics.csv")
    compute_metrics.integrated(comparison = "IGTVsNGT",
                               conditions = c("NGT", "IGT"),
                               alt_type = alt_type,
                               sample_type = sample_type,
                               result_file_path = paste0("integration_prediction_result/delegate_top2_l2_log_reg/IGTVsNGT",
                                                         "_", alt_type,
                                                         ".csv"),
                               metric_output_file_path = "integration_prediction_result/delegate_top2_l2_log_reg/metrics.csv")
  }
}


delegated_ensemble(comparison = "CFRDVsIGT",
                   conditions = c("IGT", "CFRD"),
                   prot_stacked_result_file_path = "integration_prediction_result/stacked_top2/rad_svm/CFRDVsIGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result/stacked_top2/rad_svm/CFRDVsIGT_tra.csv",
                   alt_type = "tra",
                   result_file_dir_path = "integration_prediction_result/delegate_top2_rad_svm/")
delegated_ensemble(comparison = "CFRDVsIGT",
                   conditions = c("IGT", "CFRD"),
                   prot_stacked_result_file_path = "integration_prediction_result/stacked_top2/rad_svm/CFRDVsIGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result/stacked_top2/rad_svm/CFRDVsIGT_both.csv",
                   alt_type = "both",
                   result_file_dir_path = "integration_prediction_result/delegate_top2_rad_svm/")
delegated_ensemble(comparison = "CFRDVsNGT",
                   conditions = c("NGT", "CFRD"),
                   prot_stacked_result_file_path = "integration_prediction_result/stacked_top2/rad_svm/CFRDVsNGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result/stacked_top2/rad_svm/CFRDVsNGT_tra.csv",
                   alt_type = "tra",
                   result_file_dir_path = "integration_prediction_result/delegate_top2_rad_svm/")
delegated_ensemble(comparison = "CFRDVsNGT",
                   conditions = c("NGT", "CFRD"),
                   prot_stacked_result_file_path = "integration_prediction_result/stacked_top2/rad_svm/CFRDVsNGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result/stacked_top2/rad_svm/CFRDVsNGT_both.csv",
                   alt_type = "both",
                   result_file_dir_path = "integration_prediction_result/delegate_top2_rad_svm/")
delegated_ensemble(comparison = "IGTVsNGT",
                   conditions = c("NGT", "IGT"),
                   prot_stacked_result_file_path = "integration_prediction_result/stacked_top2/rad_svm/IGTVsNGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result/stacked_top2/rad_svm/IGTVsNGT_tra.csv",
                   alt_type = "tra",
                   result_file_dir_path = "integration_prediction_result/delegate_top2_rad_svm/")
delegated_ensemble(comparison = "IGTVsNGT",
                   conditions = c("NGT", "IGT"),
                   prot_stacked_result_file_path = "integration_prediction_result/stacked_top2/rad_svm/IGTVsNGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result/stacked_top2/rad_svm/IGTVsNGT_both.csv",
                   alt_type = "both",
                   result_file_dir_path = "integration_prediction_result/delegate_top2_rad_svm/")

for(sample_type in c("train", "test")){
  for(alt_type in c("tra", "both")){
    compute_metrics.integrated(comparison = "CFRDVsIGT",
                               conditions = c("IGT", "CFRD"),
                               alt_type = alt_type,
                               sample_type = sample_type,
                               result_file_path = paste0("integration_prediction_result/delegate_top2_rad_svm/CFRDVsIGT",
                                                         "_", alt_type,
                                                         ".csv"),
                               metric_output_file_path = "integration_prediction_result/delegate_top2_rad_svm/metrics.csv")
    compute_metrics.integrated(comparison = "CFRDVsNGT",
                               conditions = c("NGT", "CFRD"),
                               alt_type = alt_type,
                               sample_type = sample_type,
                               result_file_path = paste0("integration_prediction_result/delegate_top2_rad_svm/CFRDVsNGT",
                                                         "_", alt_type,
                                                         ".csv"),
                               metric_output_file_path = "integration_prediction_result/delegate_top2_rad_svm/metrics.csv")
    compute_metrics.integrated(comparison = "IGTVsNGT",
                               conditions = c("NGT", "IGT"),
                               alt_type = alt_type,
                               sample_type = sample_type,
                               result_file_path = paste0("integration_prediction_result/delegate_top2_rad_svm/IGTVsNGT",
                                                         "_", alt_type,
                                                         ".csv"),
                               metric_output_file_path = "integration_prediction_result/delegate_top2_rad_svm/metrics.csv")
  }
}



#updated delegate ensemble with 0.7 cutoff

delegated_ensemble(comparison = "CFRDVsIGT",
                   conditions = c("IGT", "CFRD"),
                   prot_stacked_result_file_path = "integration_prediction_result/stacked_top2/l2_log_reg/CFRDVsIGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result/stacked_top2/l2_log_reg/CFRDVsIGT_tra.csv",
                   alt_type = "tra",
                   result_file_dir_path = "integration_prediction_result/delegate_top2_l2_log_reg_updated/")
delegated_ensemble(comparison = "CFRDVsIGT",
                   conditions = c("IGT", "CFRD"),
                   prot_stacked_result_file_path = "integration_prediction_result/stacked_top2/l2_log_reg/CFRDVsIGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result/stacked_top2/l2_log_reg/CFRDVsIGT_both.csv",
                   alt_type = "both",
                   result_file_dir_path = "integration_prediction_result/delegate_top2_l2_log_reg_updated/")
delegated_ensemble(comparison = "CFRDVsNGT",
                   conditions = c("NGT", "CFRD"),
                   prot_stacked_result_file_path = "integration_prediction_result/stacked_top2/l2_log_reg/CFRDVsNGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result/stacked_top2/l2_log_reg/CFRDVsNGT_tra.csv",
                   alt_type = "tra",
                   result_file_dir_path = "integration_prediction_result/delegate_top2_l2_log_reg_updated/")
delegated_ensemble(comparison = "CFRDVsNGT",
                   conditions = c("NGT", "CFRD"),
                   prot_stacked_result_file_path = "integration_prediction_result/stacked_top2/l2_log_reg/CFRDVsNGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result/stacked_top2/l2_log_reg/CFRDVsNGT_both.csv",
                   alt_type = "both",
                   result_file_dir_path = "integration_prediction_result/delegate_top2_l2_log_reg_updated/")
delegated_ensemble(comparison = "IGTVsNGT",
                   conditions = c("NGT", "IGT"),
                   prot_stacked_result_file_path = "integration_prediction_result/stacked_top2/l2_log_reg/IGTVsNGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result/stacked_top2/l2_log_reg/IGTVsNGT_tra.csv",
                   alt_type = "tra",
                   result_file_dir_path = "integration_prediction_result/delegate_top2_l2_log_reg_updated/")
delegated_ensemble(comparison = "IGTVsNGT",
                   conditions = c("NGT", "IGT"),
                   prot_stacked_result_file_path = "integration_prediction_result/stacked_top2/l2_log_reg/IGTVsNGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result/stacked_top2/l2_log_reg/IGTVsNGT_both.csv",
                   alt_type = "both",
                   result_file_dir_path = "integration_prediction_result/delegate_top2_l2_log_reg_updated/")

for(sample_type in c("train", "test")){
  for(alt_type in c("tra", "both")){
    compute_metrics.integrated(comparison = "CFRDVsIGT",
                               conditions = c("IGT", "CFRD"),
                               alt_type = alt_type,
                               sample_type = sample_type,
                               result_file_path = paste0("integration_prediction_result/delegate_top2_l2_log_reg_updated/CFRDVsIGT",
                                                         "_", alt_type,
                                                         ".csv"),
                               metric_output_file_path = "integration_prediction_result/delegate_top2_l2_log_reg_updated/metrics.csv")
    compute_metrics.integrated(comparison = "CFRDVsNGT",
                               conditions = c("NGT", "CFRD"),
                               alt_type = alt_type,
                               sample_type = sample_type,
                               result_file_path = paste0("integration_prediction_result/delegate_top2_l2_log_reg_updated/CFRDVsNGT",
                                                         "_", alt_type,
                                                         ".csv"),
                               metric_output_file_path = "integration_prediction_result/delegate_top2_l2_log_reg_updated/metrics.csv")
    compute_metrics.integrated(comparison = "IGTVsNGT",
                               conditions = c("NGT", "IGT"),
                               alt_type = alt_type,
                               sample_type = sample_type,
                               result_file_path = paste0("integration_prediction_result/delegate_top2_l2_log_reg_updated/IGTVsNGT",
                                                         "_", alt_type,
                                                         ".csv"),
                               metric_output_file_path = "integration_prediction_result/delegate_top2_l2_log_reg_updated/metrics.csv")
  }
}


##########################
##########################


#both combat

base_model_pipeline(comparison = "CFRDVsIGT",
                    conditions = c("IGT", "CFRD"),
                    data_file_path_prot = "data/formatted/proteomics/CFRDVsIGT_imputed333_mf_quantile_combat.csv",
                    data_file_path_tra = "data/formatted/rna_all/CFRDVsIGT_log_cpm_combat.csv",
                    prot_bec = "combat",
                    tra_bec = "combat",
                    dataset_replace_string_prot = "CF_EV_prot_mf_quantile_combat_",
                    dataset_replace_string_tra = "CF_EV_tra_334_combat_",
                    result_file_dir_path = "integration_prediction_result_both_combat/")
base_model_pipeline(comparison = "CFRDVsNGT",
                    conditions = c("NGT", "CFRD"),
                    data_file_path_prot = "data/formatted/proteomics/CFRDVsNGT_imputed333_mf_quantile_combat.csv",
                    data_file_path_tra = "data/formatted/rna_all/CFRDVsNGT_log_cpm_combat.csv",
                    prot_bec = "combat",
                    tra_bec = "combat",
                    dataset_replace_string_prot = "CF_EV_prot_mf_quantile_combat_",
                    dataset_replace_string_tra = "CF_EV_tra_334_combat_",
                    result_file_dir_path = "integration_prediction_result_both_combat/")
base_model_pipeline(comparison = "IGTVsNGT",
                    conditions = c("NGT", "IGT"),
                    data_file_path_prot = "data/formatted/proteomics/IGTVsNGT_imputed333_mf_quantile_combat.csv",
                    data_file_path_tra = "data/formatted/rna_all/IGTVsNGT_log_cpm_combat.csv",
                    prot_bec = "combat",
                    tra_bec = "combat",
                    dataset_replace_string_prot = "CF_EV_prot_mf_quantile_combat_",
                    dataset_replace_string_tra = "CF_EV_tra_334_combat_",
                    result_file_dir_path = "integration_prediction_result_both_combat/")

models <- read.csv("integration_prediction_result_both_combat/CFRDVsIGT.csv") %>%
  dplyr::select(model) %>%
  distinct()
models <- models$model
omics_types <- c("prot", "tra")
sample_types <- c("train", "test") 

for(omics_type in omics_types){
  for(model in models){
    for(sample_type in sample_types){
      compute_metrics(comparison = "CFRDVsIGT",
                      conditions = c("IGT", "CFRD"),
                      o_type = omics_type,
                      m = model,
                      sample_type = sample_type,
                      result_file_path = "integration_prediction_result_both_combat/CFRDVsIGT.csv",
                      metric_output_file_path = "integration_prediction_result_both_combat/metrics_base_models.csv")
      compute_metrics(comparison = "CFRDVsNGT",
                      conditions = c("NGT", "CFRD"),
                      o_type = omics_type,
                      m = model,
                      sample_type = sample_type,
                      result_file_path = "integration_prediction_result_both_combat/CFRDVsNGT.csv",
                      metric_output_file_path = "integration_prediction_result_both_combat/metrics_base_models.csv")
      compute_metrics(comparison = "IGTVsNGT",
                      conditions = c("NGT", "IGT"),
                      o_type = omics_type,
                      m = model,
                      sample_type = sample_type,
                      result_file_path = "integration_prediction_result_both_combat/IGTVsNGT.csv",
                      metric_output_file_path = "integration_prediction_result_both_combat/metrics_base_models.csv")
    }
  }  
}

metrics <- read.csv("integration_prediction_result_both_combat/metrics_base_models.csv") %>%
  arrange(Comparison, Omics_type, Type, model)
write.csv(metrics, "integration_prediction_result_both_combat/metrics_base_models_arranged.csv", row.names = FALSE)


#stack using top 2 models

classification_models <- c(
  "rf" = "Random Forest",
  "l1_log_reg" = "L1 Regularized logistic regression",
  "l2_log_reg" = "L2 Regularized logistic regression",
  "el_log_reg" = "Elastic net logistic regression",
  "sig_svm" = "Sigmoid Kernel SVM",
  "rad_svm" = "Radial Kernel SVM"
)

for(s_model in names(classification_models)){
  print(s_model)
  print(classification_models[s_model])
}

models <- read.csv("integration_prediction_result_both_combat/CFRDVsIGT.csv") %>%
  dplyr::select(model) %>%
  distinct()
models <- models$model
omics_types <- c("prot", "tra")
sample_types <- c("train", "test") 



for(s_model in names(classification_models)){
  
  print(s_model)
  
  for(o in c("prot", "tra", "both")){
    print(o)
    print("CFRDVsIGT")
    level2_stack(comparison = "CFRDVsIGT",
                 conditions = c("IGT", "CFRD"),
                 data_file_path = "integration_prediction_result_both_combat/CFRDVsIGT.csv",
                 o_type = o,
                 stack_model = s_model,
                 base_model_metrics_file_path = "integration_prediction_result_both_combat/metrics_base_models.csv",
                 output_dir_path = "integration_prediction_result_both_combat/stacked_top2/")
    print("CFRDVsNGT")
    level2_stack(comparison = "CFRDVsNGT",
                   conditions = c("NGT", "CFRD"),
                   data_file_path = "integration_prediction_result_both_combat/CFRDVsNGT.csv",
                   o_type = o,
                   stack_model = s_model,
                   base_model_metrics_file_path = "integration_prediction_result_both_combat/metrics_base_models.csv",
                   output_dir_path = "integration_prediction_result_both_combat/stacked_top2/")
    print("IGTVsNGT")
    level2_stack(comparison = "IGTVsNGT",
                 conditions = c("NGT", "IGT"),
                 data_file_path = "integration_prediction_result_both_combat/IGTVsNGT.csv",
                 o_type = o,
                 stack_model = s_model,
                 base_model_metrics_file_path = "integration_prediction_result_both_combat/metrics_base_models.csv",
                 output_dir_path = "integration_prediction_result_both_combat/stacked_top2/")
    print(paste(o, "done"))
  }  
  
  print("computing metrics")
  for(o in c("tra", "prot", "both")){
    for(sample_type in sample_types){
      print(paste(o, sample_type))
      compute_metrics(comparison = "CFRDVsIGT",
                      conditions = c("IGT", "CFRD"),
                      o_type = paste(o, "stacked", sep = "_"),
                      m = classification_models[s_model],
                      sample_type = sample_type,
                      result_file_path = paste0("integration_prediction_result_both_combat/stacked_top2/", s_model, "/CFRDVsIGT",
                                                "_", o,
                                                ".csv"),
                      metric_output_file_path = paste0("integration_prediction_result_both_combat/stacked_top2/", s_model, "/metrics.csv"))
      compute_metrics(comparison = "CFRDVsNGT",
                      conditions = c("NGT", "CFRD"),
                      o_type = paste(o, "stacked", sep = "_"),
                      m = classification_models[s_model],
                      sample_type = sample_type,
                      result_file_path = paste0("integration_prediction_result_both_combat/stacked_top2/", s_model, "/CFRDVsNGT",
                                                "_", o,
                                                ".csv"),
                      metric_output_file_path = paste0("integration_prediction_result_both_combat/stacked_top2/", s_model, "/metrics.csv"))
      compute_metrics(comparison = "IGTVsNGT",
                      conditions = c("NGT", "IGT"),
                      o_type = paste(o, "stacked", sep = "_"),
                      m = classification_models[s_model],
                      sample_type = sample_type,
                      result_file_path = paste0("integration_prediction_result_both_combat/stacked_top2/", s_model, "/IGTVsNGT",
                                                "_", o,
                                                ".csv"),
                      metric_output_file_path = paste0("integration_prediction_result_both_combat/stacked_top2/", s_model, "/metrics.csv"))
      print(paste(o, sample_type, "done------------------"))
    }
  }
  
}


delegated_ensemble(comparison = "CFRDVsIGT",
                   conditions = c("IGT", "CFRD"),
                   prot_stacked_result_file_path = "integration_prediction_result_both_combat/stacked_top2/l2_log_reg/CFRDVsIGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result_both_combat/stacked_top2/l2_log_reg/CFRDVsIGT_tra.csv",
                   alt_type = "tra",
                   result_file_dir_path = "integration_prediction_result_both_combat/delegate_top2_l2_log_reg/")
delegated_ensemble(comparison = "CFRDVsIGT",
                   conditions = c("IGT", "CFRD"),
                   prot_stacked_result_file_path = "integration_prediction_result_both_combat/stacked_top2/l2_log_reg/CFRDVsIGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result_both_combat/stacked_top2/l2_log_reg/CFRDVsIGT_both.csv",
                   alt_type = "both",
                   result_file_dir_path = "integration_prediction_result_both_combat/delegate_top2_l2_log_reg/")
delegated_ensemble(comparison = "CFRDVsNGT",
                   conditions = c("NGT", "CFRD"),
                   prot_stacked_result_file_path = "integration_prediction_result_both_combat/stacked_top2/l2_log_reg/CFRDVsNGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result_both_combat/stacked_top2/l2_log_reg/CFRDVsNGT_tra.csv",
                   alt_type = "tra",
                   result_file_dir_path = "integration_prediction_result_both_combat/delegate_top2_l2_log_reg/")
delegated_ensemble(comparison = "CFRDVsNGT",
                   conditions = c("NGT", "CFRD"),
                   prot_stacked_result_file_path = "integration_prediction_result_both_combat/stacked_top2/l2_log_reg/CFRDVsNGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result_both_combat/stacked_top2/l2_log_reg/CFRDVsNGT_both.csv",
                   alt_type = "both",
                   result_file_dir_path = "integration_prediction_result_both_combat/delegate_top2_l2_log_reg/")
delegated_ensemble(comparison = "IGTVsNGT",
                   conditions = c("NGT", "IGT"),
                   prot_stacked_result_file_path = "integration_prediction_result_both_combat/stacked_top2/l2_log_reg/IGTVsNGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result_both_combat/stacked_top2/l2_log_reg/IGTVsNGT_tra.csv",
                   alt_type = "tra",
                   result_file_dir_path = "integration_prediction_result_both_combat/delegate_top2_l2_log_reg/")
delegated_ensemble(comparison = "IGTVsNGT",
                   conditions = c("NGT", "IGT"),
                   prot_stacked_result_file_path = "integration_prediction_result_both_combat/stacked_top2/l2_log_reg/IGTVsNGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result_both_combat/stacked_top2/l2_log_reg/IGTVsNGT_both.csv",
                   alt_type = "both",
                   result_file_dir_path = "integration_prediction_result_both_combat/delegate_top2_l2_log_reg/")

for(sample_type in c("train", "test")){
  for(alt_type in c("tra", "both")){
    compute_metrics.integrated(comparison = "CFRDVsIGT",
                               conditions = c("IGT", "CFRD"),
                               alt_type = alt_type,
                               sample_type = sample_type,
                               result_file_path = paste0("integration_prediction_result_both_combat/delegate_top2_l2_log_reg/CFRDVsIGT",
                                                         "_", alt_type,
                                                         ".csv"),
                               metric_output_file_path = "integration_prediction_result_both_combat/delegate_top2_l2_log_reg/metrics.csv")
    compute_metrics.integrated(comparison = "CFRDVsNGT",
                               conditions = c("NGT", "CFRD"),
                               alt_type = alt_type,
                               sample_type = sample_type,
                               result_file_path = paste0("integration_prediction_result_both_combat/delegate_top2_l2_log_reg/CFRDVsNGT",
                                                         "_", alt_type,
                                                         ".csv"),
                               metric_output_file_path = "integration_prediction_result_both_combat/delegate_top2_l2_log_reg/metrics.csv")
    compute_metrics.integrated(comparison = "IGTVsNGT",
                               conditions = c("NGT", "IGT"),
                               alt_type = alt_type,
                               sample_type = sample_type,
                               result_file_path = paste0("integration_prediction_result_both_combat/delegate_top2_l2_log_reg/IGTVsNGT",
                                                         "_", alt_type,
                                                         ".csv"),
                               metric_output_file_path = "integration_prediction_result_both_combat/delegate_top2_l2_log_reg/metrics.csv")
  }
}


#updated ones - are rerun with delegate with only 0.7 as cutoff
delegated_ensemble(comparison = "CFRDVsIGT",
                   conditions = c("IGT", "CFRD"),
                   prot_stacked_result_file_path = "integration_prediction_result_both_combat/stacked_top2/l2_log_reg/CFRDVsIGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result_both_combat/stacked_top2/l2_log_reg/CFRDVsIGT_tra.csv",
                   alt_type = "tra",
                   result_file_dir_path = "integration_prediction_result_both_combat/delegate_top2_l2_log_reg_updated/")
delegated_ensemble(comparison = "CFRDVsIGT",
                   conditions = c("IGT", "CFRD"),
                   prot_stacked_result_file_path = "integration_prediction_result_both_combat/stacked_top2/l2_log_reg/CFRDVsIGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result_both_combat/stacked_top2/l2_log_reg/CFRDVsIGT_both.csv",
                   alt_type = "both",
                   result_file_dir_path = "integration_prediction_result_both_combat/delegate_top2_l2_log_reg_updated/")
delegated_ensemble(comparison = "CFRDVsNGT",
                   conditions = c("NGT", "CFRD"),
                   prot_stacked_result_file_path = "integration_prediction_result_both_combat/stacked_top2/l2_log_reg/CFRDVsNGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result_both_combat/stacked_top2/l2_log_reg/CFRDVsNGT_tra.csv",
                   alt_type = "tra",
                   result_file_dir_path = "integration_prediction_result_both_combat/delegate_top2_l2_log_reg_updated/")
delegated_ensemble(comparison = "CFRDVsNGT",
                   conditions = c("NGT", "CFRD"),
                   prot_stacked_result_file_path = "integration_prediction_result_both_combat/stacked_top2/l2_log_reg/CFRDVsNGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result_both_combat/stacked_top2/l2_log_reg/CFRDVsNGT_both.csv",
                   alt_type = "both",
                   result_file_dir_path = "integration_prediction_result_both_combat/delegate_top2_l2_log_reg_updated/")
delegated_ensemble(comparison = "IGTVsNGT",
                   conditions = c("NGT", "IGT"),
                   prot_stacked_result_file_path = "integration_prediction_result_both_combat/stacked_top2/l2_log_reg/IGTVsNGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result_both_combat/stacked_top2/l2_log_reg/IGTVsNGT_tra.csv",
                   alt_type = "tra",
                   result_file_dir_path = "integration_prediction_result_both_combat/delegate_top2_l2_log_reg_updated/")
delegated_ensemble(comparison = "IGTVsNGT",
                   conditions = c("NGT", "IGT"),
                   prot_stacked_result_file_path = "integration_prediction_result_both_combat/stacked_top2/l2_log_reg/IGTVsNGT_prot.csv",
                   alternate_stacked_result_file_path = "integration_prediction_result_both_combat/stacked_top2/l2_log_reg/IGTVsNGT_both.csv",
                   alt_type = "both",
                   result_file_dir_path = "integration_prediction_result_both_combat/delegate_top2_l2_log_reg_updated/")

for(sample_type in c("train", "test")){
  for(alt_type in c("tra", "both")){
    compute_metrics.integrated(comparison = "CFRDVsIGT",
                               conditions = c("IGT", "CFRD"),
                               alt_type = alt_type,
                               sample_type = sample_type,
                               result_file_path = paste0("integration_prediction_result_both_combat/delegate_top2_l2_log_reg_updated/CFRDVsIGT",
                                                         "_", alt_type,
                                                         ".csv"),
                               metric_output_file_path = "integration_prediction_result_both_combat/delegate_top2_l2_log_reg_updated/metrics.csv")
    compute_metrics.integrated(comparison = "CFRDVsNGT",
                               conditions = c("NGT", "CFRD"),
                               alt_type = alt_type,
                               sample_type = sample_type,
                               result_file_path = paste0("integration_prediction_result_both_combat/delegate_top2_l2_log_reg_updated/CFRDVsNGT",
                                                         "_", alt_type,
                                                         ".csv"),
                               metric_output_file_path = "integration_prediction_result_both_combat/delegate_top2_l2_log_reg_updated/metrics.csv")
    compute_metrics.integrated(comparison = "IGTVsNGT",
                               conditions = c("NGT", "IGT"),
                               alt_type = alt_type,
                               sample_type = sample_type,
                               result_file_path = paste0("integration_prediction_result_both_combat/delegate_top2_l2_log_reg_updated/IGTVsNGT",
                                                         "_", alt_type,
                                                         ".csv"),
                               metric_output_file_path = "integration_prediction_result_both_combat/delegate_top2_l2_log_reg_updated/metrics.csv")
  }
}
