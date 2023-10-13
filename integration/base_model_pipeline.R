library(tidyverse)
library(sva)
library(readxl)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/CysticFibrosisGroup/ExoCF/CFRD_EV_biomarker/"
setwd(base_dir)

source("prediction_pipeline/cm_logistic_regression.R")
source("prediction_pipeline/cm_svm.R")
source("prediction_pipeline/cm_rf.R")
source("integration/run_all_models.R")

# data_file_path_prot <- "data/formatted/proteomics/CFRDVsIGT_imputed333_mf_quantile_combat.csv"
# data_file_path_tra <- "data/formatted/CFRDVsIGT_umi_counts_combat_processed.csv"
# comparison = "CFRDVsIGT"
# conditions = c("IGT", "CFRD")

convert_format_and_select_biomarkers <- function(data, biomarkers){
  data <- as.data.frame(t(as.matrix(data)))
  
  # / or - in biomarkers would have been changed to .
  # section with features_with_slash, is used to change 'biomarkers' back to the one with / or -
  features_with_slash <- colnames(data)[grepl("/", colnames(data), fixed = TRUE)] 
  for(f in features_with_slash){
    f_replaced <- gsub("/|-", ".", f) 
    if(f_replaced %in% biomarkers){
      biomarkers[biomarkers == f_replaced] = f
    }
  }
  # end of section with features_with_slash
  
  biomarkers <- gsub(".", "-", biomarkers, fixed = TRUE)
  data <- as.data.frame(t(as.matrix(data)))
  data <- data[biomarkers, ]
  
  return (data)
}

base_model_pipeline <- function(comparison, conditions, 
                                data_file_path_prot,
                                data_file_path_tra){
  
  phenotype_file_path_prot <- "data/formatted/prot_phenotype_333.txt"
  phenotype_file_path_tra <- "data/formatted/phenotype.txt"
  mapping_file_path <- "data/formatted/pt_mapping_and_common_sample_name.txt"
  
  best_features_file_path <- "data/selected_features/best_features_with_is_best.csv"
  dataset_replace_string_prot = "CF_EV_prot_mf_quantile_combat_"
  dataset_replace_string_tra = "CF_EV_tra_combat_"
  
  best_features <- read.csv(best_features_file_path)  
  
  best_features_prot <- best_features %>%
    mutate(dataset_id = gsub(dataset_replace_string_prot, "", dataset_id)) %>%
    filter(is_best == 1, dataset_id == comparison)
  best_features_prot <- strsplit(best_features_prot$biomarkers, split = "|", fixed = TRUE)[[1]] 
  best_features_tra <- best_features %>%
    mutate(dataset_id = gsub(dataset_replace_string_tra, "", dataset_id)) %>%
    filter(is_best == 1, dataset_id == comparison)
  best_features_tra <- strsplit(best_features_tra$biomarkers, split = "|", fixed = TRUE)[[1]] 
  
  data_prot <- read.table(data_file_path_prot, header=TRUE, sep=",", row.names=1, skip=0,
                          nrows=-1, comment.char="", fill=TRUE, na.strings = "NA") 
  data_prot <- convert_format_and_select_biomarkers(data_prot, best_features_prot)
  
  data_tra <- read.table(data_file_path_tra, header=TRUE, sep=",", row.names=1, skip=0,
                         nrows=-1, comment.char="", fill=TRUE, na.strings = "NA") 
  data_tra <- convert_format_and_select_biomarkers(data_tra, best_features_tra)
  
  label_prot <- read.table(phenotype_file_path_prot, header=TRUE, sep="\t") %>%
    rename("Label" = comparison) %>%
    filter(Label %in% conditions) %>%
    dplyr::select(Sample, Label, country) %>%
    arrange(Sample)
  label_tra <- read.table(phenotype_file_path_tra, header=TRUE, sep="\t") %>%
    rename("Label" = comparison) %>%
    filter(Label %in% conditions) %>%
    dplyr::select(Sample, Label, country) %>%
    arrange(Sample)
  
  prot_tra_mapping <- read.table(mapping_file_path, header=TRUE, sep="\t")
  
  label_prot <- label_prot %>%
    inner_join(prot_tra_mapping, by = c("Sample" = "prot_sample_name"))%>%
    dplyr::select(-c(tra_sample_name))
  
  #mapping file contains mapping from single tra sample to multiple prot samples - since in some cases,
  #   there are prot replicates. This causes repeats of tra sample name in the mapping file
  #   In these cases, choose sample that comes up first
  
  label_tra <- label_tra %>%
    inner_join(prot_tra_mapping, by = c("Sample" = "tra_sample_name")) %>%
    dplyr::select(-c(prot_sample_name))
  label_tra <- label_tra %>%
    group_by(Sample) %>%
    mutate(replicate = row_number()) %>%
    ungroup()
  label_tra <- label_tra %>%
    filter(replicate == 1) %>%
    dplyr::select(-c(replicate))

  data_prot <- data_prot[, label_prot$Sample]
  colnames(data_prot) <- label_prot$sample_name
  
  data_tra <- data_tra[, label_tra$Sample]
  colnames(data_tra) <- label_tra$sample_name
  
  label_prot <- label_prot %>%
    dplyr::select(-c(Sample)) %>%
    dplyr::rename(c("Sample" = "sample_name")) %>%
    arrange(Sample)
  label_tra <- label_tra %>%
    dplyr::select(-c(Sample)) %>%
    dplyr::rename(c("Sample" = "sample_name")) %>%
    arrange(Sample)
  
  #all_labels variable is necessary because there are some unique samples in prot and tra
  all_labels <- rbind(label_prot, label_tra) %>%
    distinct() %>%
    arrange(Sample) %>%
    mutate(Label_cohort = paste(Label, country, sep = "_"))
  
  all_labels$Sample[!all_labels$Sample %in% label_prot$Sample]
  all_labels$Sample[!all_labels$Sample %in% label_tra$Sample]
  
  length(unique(all_labels$Sample))
  length(unique(label_prot$Sample))
  length(unique(label_tra$Sample))
  
  all_labels[!all_labels$Sample %in% label_prot$Sample, ]
  all_labels[!all_labels$Sample %in% label_tra$Sample,]
  
  data_prot <- as.data.frame(t(data_prot))
  data_tra <- as.data.frame(t(data_tra))
  
  set.seed(1000)
  k_val <- 5
  times_val <- 6
  k_prod_times <- k_val * times_val
  train_index <- caret::createMultiFolds(y = all_labels$Label_cohort, k = k_val, times = times_val)
  
  for (i in c(1:k_prod_times)){
    # print(i)
    all_labels.train <- all_labels[train_index[[i]], , drop = FALSE]
    all_labels.test <- all_labels[-train_index[[i]], , drop = FALSE]
    
    label_prot.train <- label_prot[label_prot$Sample %in% all_labels.train$Sample, ]
    data_prot.train <- data_prot[label_prot.train$Sample, ]
    label_prot.test <- label_prot[label_prot$Sample %in% all_labels.test$Sample, ]
    data_prot.test <- data_prot[label_prot.test$Sample, ]
    
    label_tra.train <- label_tra[label_tra$Sample %in% all_labels.train$Sample, ]
    data_tra.train <- data_tra[label_tra.train$Sample, ]
    label_tra.test <- label_tra[label_tra$Sample %in% all_labels.test$Sample, ]
    data_tra.test <- data_tra[label_tra.test$Sample, ]
    result_df_prot <- cbind(omics_type = "prot", 
                            run_all_models(data_prot.train, label_prot.train,
                                           data_prot.test, label_prot.test,
                                           conditions))
    # print("prot done")
    result_df_tra <- cbind(omics_type = "tra", 
                           run_all_models(data_tra.train, label_tra.train,
                                          data_tra.test, label_tra.test,
                                          conditions))
    # print("tra done")
    result_df <- cbind(iter = i,
                       rbind(result_df_prot,
                             result_df_tra))
    
    if(i == 1){
      result_df_all <- result_df
    } else{
      result_df_all <- rbind(result_df_all, result_df)
    }
  }
  
  result_file_dir_path <- "integration_prediction_result/"
  result_file_name <- paste0(comparison, ".csv")
  if(!dir.exists(result_file_dir_path)){
    dir.create(result_file_dir_path, recursive = TRUE)
  }
  write.csv(format(result_df_all, digits = 3), paste0(result_file_dir_path, result_file_name), 
            row.names = FALSE)
}
