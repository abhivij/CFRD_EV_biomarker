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


