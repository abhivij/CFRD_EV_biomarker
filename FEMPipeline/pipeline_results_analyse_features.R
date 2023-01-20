setwd("~/UNSW/VafaeeLab/CysticFibrosisGroup/ExoCF/CFRD_EV_biomarker/FEMPipeline/")
source("dataset_pipeline_arguments.R")
source("utils.R")
library(tidyverse)
library(viridis)
library(ComplexHeatmap)
library(UpSetR)

explore_common_features <- function(dparg_id, best_fsm_vec, 
                                    min_iter_feature_presence,
                                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                                    results_dir = "fem_pipeline_results",
                                    dir_path = "plots/FEMPipeline/common_features_upset"){
  
  ds <- dataset_pipeline_arguments[[dparg_id]]
  dataset_id <- paste(ds$dataset_id, ds$classification_criteria, sep = "_")
  print(dataset_id)
  
  features_file <- paste(dataset_id, "features.csv", sep = "_")
  features_file <- paste(results_dir, features_file, sep = "/")
  print(features_file)
  
  selected_features <- list()
  max_size <- 0  #for upset plot max bar length
  for(best_fsm in best_fsm_vec){
    features_info <- read.table(features_file, sep = ',', header = TRUE)
    
    features_info <- features_info %>%
      filter(FSM == best_fsm) %>%
      select(-c(FSM, Iter))
    print(best_fsm)
    print(dim(features_info))
    
    print(dim(features_info[,colSums(features_info) >= min_iter_feature_presence, drop = FALSE]))
    selected_features[[best_fsm]] <- colnames(features_info[,colSums(features_info) >= min_iter_feature_presence, drop = FALSE])   
    
    size <- length(selected_features[[best_fsm]])
    if(size > max_size){
      max_size <- size
    }
  }
  ###########################################
  #upsetplot
  file_name <- paste(dataset_id, min_iter_feature_presence, 
                     "upset.png", sep = "_")
  
  if(!dir.exists(dir_path)){
    dir.create(dir_path, recursive = TRUE)
  }
  
  print("creating upset plot...")
  print(paste(dir_path, file_name, sep = "/"))
  # dev.copy(png, paste(dir_path, file_name, sep = "/"), 
  #          units = "cm", width = 30, height = 15, res = 1200)
  png(filename = paste(dir_path, file_name, sep = "/"),
      units = "cm", width = 30, height = 15, res = 1200)
  print({
    upset(fromList(selected_features), set_size.show = TRUE,  
          set_size.scale_max = max_size + 20)  
  })
  grid.text(paste(dataset_id, min_iter_feature_presence),
            x = 0.6, y = 0.95)
  dev.off()
  
  
  ####write best features to file
  output_dir <- "../data/selected_features/"
  file_name <- "best_FSM_common_features.csv"
  
  selected_features_df <- data.frame(matrix(nrow = 0, ncol = 4))
  for(fsm in names(selected_features)){
    print(fsm)
    row <- c(dataset_id,
             fsm,
             min_iter_feature_presence,
             paste(selected_features[[fsm]], collapse = "|")
    )
    selected_features_df <- rbind(selected_features_df, row)
  }
  colnames(selected_features_df) <- c("datasetid",
                                      "fsm",
                                      "miniter_feature_presence",
                                      "biomarkers")
  if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
  }
  
  if(length(best_fsm_vec) == dim(selected_features_df)[1]){
    file_path <- paste0(output_dir, file_name)
    write.table(selected_features_df, 
                file_path,
                sep = ",",
                row.names = FALSE, append = TRUE,
                col.names = !file.exists(file_path))    
  } else{
    print("Not writing selected features to file")
    print("all FSMs dont have common features across iter")
  }
}




get_features_from_df <- function(features_df, fsm){
  features <- features_df[features_df$fsm == fsm, "biomarkers"]
  if(length(features) != 0){
    features <- strsplit(features, split = "|", fixed = TRUE)[[1]]  
  } else{
    features <- c()
  }
  features
}

write_subset_file <- function(data, features, subset_file_path){
  
  #"-" and "/" both in features are replaced with "."
  #so replacing "." with "-" alone as done after the below 'for loop' is not sufficient
  #therefore special handling for features with "/"
  features_with_slash <- rownames(data)[grepl("/", rownames(data), fixed = TRUE)] 
  for(f in features_with_slash){
    f_replaced <- gsub("/|-", ".", f) 
    if(f_replaced %in% features){
      features[features == f_replaced] = f
    }
  }
  
  data_sub <- data[gsub(".", "-", features, fixed = TRUE),]
  print(dim(data_sub))
  print(sum(is.na(data_sub)))
  print(subset_file_path)
  write.csv(data_sub, subset_file_path)
}

# dparg_id = 42
# min_iter_feature_presence = 28
# subset_creation_criteria <- list("i"= c("ranger_impu_cor"))
# subset_file_name_substr = "ranger"
# 
# create_all_common = TRUE


create_data_subsets <- function(dparg_id, 
                                min_iter_feature_presence,
                                subset_creation_criteria,
                                subset_file_name_substr = "common3",
                                combat = FALSE,
                                create_all_common = TRUE,
                                dataset_pipeline_arguments = dataset_pipeline_arguments,
                                data_file_path = NA){
  
  ds <- dataset_pipeline_arguments[[dparg_id]]
  dataset_id <- paste(ds$dataset_id, ds$classification_criteria, sep = "_")
  print(dataset_id)
  
  output_dir <- "../data/selected_features/"
  file_name <- "best_FSM_common_features.csv"
  features_df <- read.csv(paste0(output_dir, file_name))
  
  features_df <- features_df %>%
    filter(datasetid == dataset_id) %>%
    filter(miniter_feature_presence == min_iter_feature_presence)
  
  if(is.na(data_file_path)){
    if(!combat){
      data_file_path = "../data/formatted/umi_counts.csv"  
    } else{
      data_file_path = "../data/formatted/umi_counts_combat_seq.csv"
    }    
  }
  
  ######create new data subsets
  
  data <- read.table(data_file_path, header=TRUE, sep=",", row.names=1, skip=0,
                       nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")   
  output_dir <- "../data/formatted/subset/"

  if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
  }
  
  best_features_df <- data.frame(matrix(nrow = 0, ncol = 5))
  
  #create datasubset with features common in all FSMs
  if(create_all_common == TRUE){
    intersect_list <- list()
    i <- 1
    for(i_fsm in features_df$fsm){
      intersect_list[[i]] <- get_features_from_df(features_df, i_fsm)
      i <- i + 1
    }
    features <- Reduce(intersect, intersect_list)  
    write_subset_file(data, features, 
                      subset_file_path = paste0(output_dir, dataset_id, "_all_common_",
                                                min_iter_feature_presence, ".csv"))
    
    best_features_df <- 
      rbind(
        best_features_df,
        data.frame(
          dataset_id = dataset_id,
          description = "all_common",
          min_iter_feature_presence = min_iter_feature_presence,
          biomarkers = paste(features, collapse = "|"),
          size = length(features)
        )        
      )
    
    if(sum(grepl("piR", features, fixed = TRUE)) > 0){
      print("all common features contain piRNA !")
      features_without_pir <- features[!grepl("piR", features, fixed = TRUE)]
      write_subset_file(data, features_without_pir, 
                        subset_file_path = paste0(output_dir, dataset_id, 
                                                  "_all_common_nopir_",
                                                  min_iter_feature_presence, ".csv"))
      best_features_df <- 
        rbind(
          best_features_df,
          data.frame(
            dataset_id = dataset_id,
            description = "all_common_nopir",
            min_iter_feature_presence = min_iter_feature_presence,
            biomarkers = paste(features_without_pir, collapse = "|"),
            size = length(features_without_pir)
          )        
        )
    }
  }
  
  if(length(subset_creation_criteria) != 0){
    # example subset_creation_criteria
    # subset_creation_criteria <- list("i"= c("t-test",
    #                                         "wilcoxontest",
    #                                         "wilcoxontest_pval_0.005"),
    #                                  "d"= c("mrmr75"))
    intersect_list <- list()
    i <- 1
    for(i_fsm in subset_creation_criteria[["i"]]){
      print(i_fsm)
      intersect_list[[i]] <- get_features_from_df(features_df, i_fsm)
      i <- i + 1
    }
    
    features <- Reduce(intersect, intersect_list)
    #currently handles case when "d" has single value only
    if(length(subset_creation_criteria[["d"]]) == 1){
      features <- setdiff(features,
                          get_features_from_df(features_df, subset_creation_criteria[["d"]])
      )
    }
    
    write_subset_file(data, features, 
                      subset_file_path = paste0(output_dir, dataset_id, 
                                                "_", subset_file_name_substr, "_",
                                                min_iter_feature_presence, ".csv"))
    best_features_df <- 
      rbind(
        best_features_df,
        data.frame(
          dataset_id = dataset_id,
          description = subset_file_name_substr,
          min_iter_feature_presence = min_iter_feature_presence,
          biomarkers = paste(features, collapse = "|"),
          size = length(features)
        )        
      )
    
    if(sum(grepl("piR", features, fixed = TRUE)) > 0){
      print("subset creation criteria features contain piRNA !")
      features_without_pir <- features[!grepl("piR", features, fixed = TRUE)]
      write_subset_file(data, features_without_pir, 
                        subset_file_path = paste0(output_dir, dataset_id, 
                                                  "_", subset_file_name_substr, 
                                                  "_nopir_",
                                                  min_iter_feature_presence, ".csv"))
      best_features_df <- 
        rbind(
          best_features_df,
          data.frame(
            dataset_id = dataset_id,
            description = paste0(subset_file_name_substr, "_nopir"),
            min_iter_feature_presence = min_iter_feature_presence,
            biomarkers = paste(features_without_pir, collapse = "|"),
            size = length(features_without_pir)
          )        
        )
    }
    
  }
  
  output_dir <- "../data/selected_features/"
  file_name <- "best_features.csv"
  if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
  }
  file_path <- paste0(output_dir, file_name)
  write.table(best_features_df, 
              file_path,
              sep = ",",
              row.names = FALSE, append = TRUE,
              col.names = !file.exists(file_path))
  
}



##############common features across all

# best_features <- read.csv("Data/selected_features/best_features_with_add_col.csv")
# best_features <- best_features %>%
#   filter(is_best == 1)
# 
# best_features_pr <- best_features %>%
#   filter(grepl("proteomic", dataset_id, fixed = TRUE))
# i <- 1
# intersect_list <- list()
# for(ds in best_features_pr$dataset_id){
#   print(ds)
#   features <- best_features_pr[best_features_pr$dataset_id == ds, "biomarkers"]
#   if(length(features) != 0){
#     features <- strsplit(features, split = "|", fixed = TRUE)[[1]]  
#   } else{
#     features <- c()
#   }
#   intersect_list[[i]] <- features
#   i <- i + 1
# }
# features <- Reduce(intersect, intersect_list)
# names(intersect_list) <- best_features_pr$dataset_id
# upset(fromList(intersect_list), set_size.show = TRUE)  
# 
# best_features_tr <- best_features %>%
#   filter(grepl("transcriptomic", dataset_id, fixed = TRUE))
# i <- 1
# intersect_list <- list()
# for(ds in best_features_tr$dataset_id){
#   print(ds)
#   features <- best_features_tr[best_features_tr$dataset_id == ds, "biomarkers"]
#   if(length(features) != 0){
#     features <- strsplit(features, split = "|", fixed = TRUE)[[1]]  
#   } else{
#     features <- c()
#   }
#   intersect_list[[i]] <- features
#   i <- i + 1
# }
# features <- Reduce(intersect, intersect_list)
# names(intersect_list) <- best_features_tr$dataset_id
# upset(fromList(intersect_list), set_size.show = TRUE)  
# dev.off()


#####################################

#initial cohort with new quantified results

explore_common_features(dparg_id = 25,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr_perc50"),
                        min_iter_feature_presence = 27,
                        results_dir = "../fem_pipeline_results_AU",
                        dir_path = "../plots/FEMPipeline/common_features_upset/AU")




#######################
#AU tmm : CFRDVsIGT

explore_common_features(dparg_id = 37,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr_perc50", "mrmr75",
                                         "ranger_pos_impu_cor"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_AU",
                        dir_path = "../plots/FEMPipeline_AU_tmm/common_features_upset")

create_data_subsets(dparg_id = 37,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75", "ranger_pos_impu_cor")),
                    subset_file_name_substr = "mrmr75_rangerpos",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 37,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 37,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("ranger_pos_impu_cor")),
                    subset_file_name_substr = "ranger_pos",
                    create_all_common = FALSE)




explore_common_features(dparg_id = 41,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr75", "mrmr_perc50",
                                         "wilcoxontest"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_AU",
                        dir_path = "../plots/FEMPipeline_AU_tmm/common_features_upset")
create_data_subsets(dparg_id = 41,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE)




explore_common_features(dparg_id = 45,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr75", "wilcoxontest", 
                                         "t-test"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_AU",
                        dir_path = "../plots/FEMPipeline_AU_tmm/common_features_upset")
create_data_subsets(dparg_id = 45,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE)




#######################
#All tmm : CFRDVsIGT
explore_common_features(dparg_id = 13,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("ga_rf", "mrmr75",
                                         "ranger_pos_impu_cor"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results",
                        dir_path = "../plots/FEMPipeline_tmm_all/common_features_upset")

#All tmm : CFRDVsNGT
explore_common_features(dparg_id = 17,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr10", "RF_RFE",
                                         "mrmr75", "ranger_pos_impu_cor"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results",
                        dir_path = "../plots/FEMPipeline_tmm_all/common_features_upset")

# included ranger_pos_impu_cor in above because out of other 3 only mrmr75 has non-zero features, 
#        and atleast 2 non-zero groups are required to create upset plot


#All tmm : IGTVsNGT
explore_common_features(dparg_id = 21,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("ga_rf", "mrmr_perc50",
                                         "ranger_pos_impu_cor", "RF_RFE", "mrmr10",
                                         "t-test", "mrmr75", "wilcoxontest"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results",
                        dir_path = "../plots/FEMPipeline_tmm_all/common_features_upset")
#using multiple methods because the top 3 didn't give good number of features.
#mrmr_perc50 had large number of features, others very few features

explore_common_features(dparg_id = 21,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("ga_rf", "mrmr_perc50",
                                         "ranger_pos_impu_cor", "RF_RFE", "mrmr10",
                                         "t-test", "mrmr75", "wilcoxontest"),
                        min_iter_feature_presence = 27,
                        results_dir = "../fem_pipeline_results",
                        dir_path = "../plots/FEMPipeline_tmm_all/common_features_upset")

explore_common_features(dparg_id = 21,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("ga_rf", "mrmr_perc50",
                                         "ranger_pos_impu_cor", "RF_RFE", "mrmr10",
                                         "t-test", "mrmr75", "wilcoxontest"),
                        min_iter_feature_presence = 29,
                        results_dir = "../fem_pipeline_results",
                        dir_path = "../plots/FEMPipeline_tmm_all/common_features_upset")



create_data_subsets(dparg_id = 13,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE)

create_data_subsets(dparg_id = 17,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE)

create_data_subsets(dparg_id = 21,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE)



#######################
#AU log tmm no scaling : CFRDVsIGT
explore_common_features(dparg_id = 65,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr_perc50", "mrmr75", "RF_RFE", 
                                         "ranger_pos_impu_cor"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_AU_logtmm/",
                        dir_path = "../plots/FEMPipeline_AU_noscaling_tmm/common_features_upset")

#AU log tmm no scaling : CFRDVsNGT
explore_common_features(dparg_id = 69,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr75", "mrmr_perc50", "RF_RFE", 
                                         "wilcoxontest"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_AU_logtmm/",
                        dir_path = "../plots/FEMPipeline_AU_noscaling_tmm/common_features_upset")

#AU log tmm no scaling : IGTVsNGT
explore_common_features(dparg_id = 73,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr75", "wilcoxontest", "t-test", 
                                         "ranger_pos_impu_cor"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_AU_logtmm/",
                        dir_path = "../plots/FEMPipeline_AU_noscaling_tmm/common_features_upset")


create_data_subsets(dparg_id = 65,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 69,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 73,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE)



#######################

#AU adult log tmm : CFRDVsIGT
explore_common_features(dparg_id = 83,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr_perc50", "mrmr75", "mrmr100"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_AU_adult_logtmm",
                        dir_path = "../plots/fem_pipeline_results_AU_adult_logtmm/common_features_upset")

#AU adult log tmm : CFRDVsNGT
explore_common_features(dparg_id = 87,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr75", "mrmr100", "mrmr_perc50"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_AU_adult_logtmm",
                        dir_path = "../plots/fem_pipeline_results_AU_adult_logtmm/common_features_upset")

#AU adult log tmm : IGTVsNGT
explore_common_features(dparg_id = 91,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr100", "mrmr75", "wilcoxontest"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_AU_adult_logtmm",
                        dir_path = "../plots/fem_pipeline_results_AU_adult_logtmm/common_features_upset")


create_data_subsets(dparg_id = 83,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 87,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 91,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 91,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE)



#######################

#DK adult log tmm : CFRDVsIGT
explore_common_features(dparg_id = 103,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("ga_rf", "ranger_pos_impu_cor", "mrmr_perc50", "t-test", "mrmr100"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_DK_adult_logtmm",
                        dir_path = "../plots/fem_pipeline_results_DK_adult_logtmm/common_features_upset")

#DK adult log tmm : CFRDVsNGT
explore_common_features(dparg_id = 107,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("ranger_pos_impu_cor", "mrmr10", "mrmr_perc50",
                                         "mrmr75"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_DK_adult_logtmm",
                        dir_path = "../plots/fem_pipeline_results_DK_adult_logtmm/common_features_upset")

#DK adult log tmm : IGTVsNGT
explore_common_features(dparg_id = 111,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("t-test", "wilcoxontest", "mrmr10", "mrmr_perc50",
                                         "ranger_pos_impu_cor", "mrmr100"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_DK_adult_logtmm",
                        dir_path = "../plots/fem_pipeline_results_DK_adult_logtmm/common_features_upset")


create_data_subsets(dparg_id = 103,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 107,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 111,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE)



#######################

#AU adult log cpm : CFRDVsIGT
explore_common_features(dparg_id = 121,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("ga_rf", "mrmr75", "t-test"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_AU_adult_logcpm",
                        dir_path = "../plots/fem_pipeline_results_AU_adult_logcpm/common_features_upset")

#AU adult log cpm : CFRDVsNGT
explore_common_features(dparg_id = 125,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr75", "mrmr100", "mrmr_perc50"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_AU_adult_logcpm",
                        dir_path = "../plots/fem_pipeline_results_AU_adult_logcpm/common_features_upset")

#AU adult log cpm : IGTVsNGT
explore_common_features(dparg_id = 129,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr75", "mrmr100", "mrmr10"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_AU_adult_logcpm",
                        dir_path = "../plots/fem_pipeline_results_AU_adult_logcpm/common_features_upset")


create_data_subsets(dparg_id = 121,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 125,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 129,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE)


#######################

#AU adult log  : CFRDVsIGT
explore_common_features(dparg_id = 133,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr100", "mrmr_perc50", "mrmr75"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_AU_adult_log",
                        dir_path = "../plots/fem_pipeline_results_AU_adult_log/common_features_upset")

#AU adult log  : CFRDVsNGT
explore_common_features(dparg_id = 137,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("t-test", "mrmr10", "ranger_pos_impu_cor", "ga_rf", "mrmr100"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_AU_adult_log",
                        dir_path = "../plots/fem_pipeline_results_AU_adult_log/common_features_upset")

#AU adult log  : IGTVsNGT
explore_common_features(dparg_id = 141,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("ranger_pos_impu_cor", "ga_rf", "mrmr_perc50", "RF_RFE", "mrmr75"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_AU_adult_log",
                        dir_path = "../plots/fem_pipeline_results_AU_adult_log/common_features_upset")


create_data_subsets(dparg_id = 133,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 137,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 141,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr_perc50")),
                    subset_file_name_substr = "mrmr_perc50",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 141,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE)

#######################

#AU adult combat seq log cpm : CFRDVsIGT
explore_common_features(dparg_id = 145,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("ranger_pos_impu_cor", "ga_rf", "mrmr_perc50", "mrmr100"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_AU_adult_combat_logcpm",
                        dir_path = "../plots/fem_pipeline_results_AU_adult_combat_logcpm/common_features_upset")

#AU adult combat seq log cpm : CFRDVsNGT
explore_common_features(dparg_id = 149,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr75", "mrmr100"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_AU_adult_combat_logcpm",
                        dir_path = "../plots/fem_pipeline_results_AU_adult_combat_logcpm/common_features_upset")

#AU adult combat seq log cpm : IGTVsNGT
explore_common_features(dparg_id = 153,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr75", "mrmr100"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_AU_adult_combat_logcpm",
                        dir_path = "../plots/fem_pipeline_results_AU_adult_combat_logcpm/common_features_upset")

create_data_subsets(dparg_id = 145,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr_perc50")),
                    subset_file_name_substr = "mrmr_perc50",
                    create_all_common = FALSE, combat = TRUE)
create_data_subsets(dparg_id = 145,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE, combat = TRUE)
create_data_subsets(dparg_id = 149,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE, combat = TRUE)
create_data_subsets(dparg_id = 153,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE, combat = TRUE)

#######################


#AU adult combat seq log tmm : CFRDVsIGT
explore_common_features(dparg_id = 157,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("ranger_pos_impu_cor", "mrmr_perc50", "ga_rf", "mrmr75"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_AU_adult_combat_logtmm",
                        dir_path = "../plots/fem_pipeline_results_AU_adult_combat_logtmm/common_features_upset")

#AU adult combat seq log tmm : CFRDVsNGT
explore_common_features(dparg_id = 161,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr100", "mrmr75"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_AU_adult_combat_logtmm",
                        dir_path = "../plots/fem_pipeline_results_AU_adult_combat_logtmm/common_features_upset")

#AU adult combat seq log tmm : IGTVsNGT
explore_common_features(dparg_id = 165,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr75", "mrmr100"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_AU_adult_combat_logtmm",
                        dir_path = "../plots/fem_pipeline_results_AU_adult_combat_logtmm/common_features_upset")


create_data_subsets(dparg_id = 157,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr_perc50")),
                    subset_file_name_substr = "mrmr_perc50",
                    create_all_common = FALSE, combat = TRUE)
create_data_subsets(dparg_id = 157,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE, combat = TRUE)
create_data_subsets(dparg_id = 161,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE, combat = TRUE)
create_data_subsets(dparg_id = 165,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE, combat = TRUE)

#######################

#AU adult prefiltered combat seq log tmm 
#CFRDVsIGT
explore_common_features(dparg_id = 199,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr100", "mrmr75", "wilcoxontest", "mrmr10"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_AU_prefiltered_adult_combat_logtmm",
                        dir_path = "../plots/fem_pipeline_results_AU_prefiltered_adult_combat_logtmm/common_features_upset")

#CFRDVsNGT
explore_common_features(dparg_id = 203,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr75", "mrmr100", "mrmr10", "RF_RFE"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_AU_prefiltered_adult_combat_logtmm",
                        dir_path = "../plots/fem_pipeline_results_AU_prefiltered_adult_combat_logtmm/common_features_upset")

#IGTVsNGT
explore_common_features(dparg_id = 207,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("RF_RFE", "mrmr75", "mrmr100", "mrmr10"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_AU_prefiltered_adult_combat_logtmm",
                        dir_path = "../plots/fem_pipeline_results_AU_prefiltered_adult_combat_logtmm/common_features_upset")

#######################

#AU adult prefiltered combat seq log cpm 
#CFRDVsIGT
explore_common_features(dparg_id = 211,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr75", "ranger_pos_impu_cor", "ga_rf", "mrmr100"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_AU_prefiltered_adult_combat_logcpm",
                        dir_path = "../plots/fem_pipeline_results_AU_prefiltered_adult_combat_logcpm/common_features_upset")

#CFRDVsNGT
explore_common_features(dparg_id = 215,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr75", "mrmr100", "mrmr_perc50", "ranger_pos_impu_cor"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_AU_prefiltered_adult_combat_logcpm",
                        dir_path = "../plots/fem_pipeline_results_AU_prefiltered_adult_combat_logcpm/common_features_upset")

#IGTVsNGT
explore_common_features(dparg_id = 219,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr100", "mrmr10", "wilcoxontest", "mrmr75", "RF_RFE"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_AU_prefiltered_adult_combat_logcpm",
                        dir_path = "../plots/fem_pipeline_results_AU_prefiltered_adult_combat_logcpm/common_features_upset")


create_data_subsets(dparg_id = 199,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE,
                    data_file_path = "../data/formatted/CFRDVsIGT_umi_counts_combat_seq.csv")
create_data_subsets(dparg_id = 203,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE,
                    data_file_path = "../data/formatted/CFRDVsNGT_umi_counts_combat_seq.csv")
create_data_subsets(dparg_id = 207,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE,
                    data_file_path = "../data/formatted/IGTVsNGT_umi_counts_combat_seq.csv")

create_data_subsets(dparg_id = 211,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE,
                    data_file_path = "../data/formatted/CFRDVsIGT_umi_counts_combat_seq.csv")
create_data_subsets(dparg_id = 215,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE,
                    data_file_path = "../data/formatted/CFRDVsNGT_umi_counts_combat_seq.csv")
create_data_subsets(dparg_id = 219,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE,
                    data_file_path = "../data/formatted/IGTVsNGT_umi_counts_combat_seq.csv")



###########
###########

#adult log

explore_common_features(dparg_id = 235,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr_perc50", "mrmr100"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_adult_log",
                        dir_path = "../plots/fem_pipeline_results_adult_log/common_features_upset")
explore_common_features(dparg_id = 239,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr10", "RF_RFE", "t-test", "mrmr100"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_adult_log",
                        dir_path = "../plots/fem_pipeline_results_adult_log/common_features_upset")
explore_common_features(dparg_id = 243,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("ranger_pos_impu_cor", "ga_rf", "mrmr_perc50",
                                         "t-test", "mrmr10", "RF_RFE", "mrmr100"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_adult_log",
                        dir_path = "../plots/fem_pipeline_results_adult_log/common_features_upset")

create_data_subsets(dparg_id = 235,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr_perc50")),
                    subset_file_name_substr = "mrmr_perc50",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 235,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE)

create_data_subsets(dparg_id = 239,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE)

create_data_subsets(dparg_id = 243,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr_perc50")),
                    subset_file_name_substr = "mrmr_perc50",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 243,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE)

###########
###########

#adult log_cpm

explore_common_features(dparg_id = 247,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("ga_rf", "mrmr100", "mrmr_perc50"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_adult_log_cpm",
                        dir_path = "../plots/fem_pipeline_results_adult_log_cpm/common_features_upset")
explore_common_features(dparg_id = 251,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr10", "wilcoxontest", "t-test", "RF_RFE", "mrmr75"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_adult_log_cpm",
                        dir_path = "../plots/fem_pipeline_results_adult_log_cpm/common_features_upset")
explore_common_features(dparg_id = 255,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("ga_rf", "ranger_pos_impu_cor", "mrmr_perc50", "mrmr100"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_adult_log_cpm",
                        dir_path = "../plots/fem_pipeline_results_adult_log_cpm/common_features_upset")

create_data_subsets(dparg_id = 247,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 251,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 255,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr_perc50")),
                    subset_file_name_substr = "mrmr_perc50",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 255,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE)
###########
###########

#adult log_tmm

explore_common_features(dparg_id = 259,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("ga_rf", "mrmr_perc50", "t-test", "wilcoxontest", 
                                         "mrmr10", "RF_RFE", "mrmr75"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_adult_log_tmm",
                        dir_path = "../plots/fem_pipeline_results_adult_log_tmm/common_features_upset")
explore_common_features(dparg_id = 263,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr10", "ga_rf", "t-test", "mrmr_perc50",
                                         "ranger_pos_impu_cor", "mrmr75"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_adult_log_tmm",
                        dir_path = "../plots/fem_pipeline_results_adult_log_tmm/common_features_upset")
explore_common_features(dparg_id = 267,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("ga_rf", "ranger_pos_impu_cor", "mrmr_perc50", "mrmr100"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_adult_log_tmm",
                        dir_path = "../plots/fem_pipeline_results_adult_log_tmm/common_features_upset")

create_data_subsets(dparg_id = 259,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr_perc50")),
                    subset_file_name_substr = "mrmr_perc50",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 259,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 263,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr_perc50")),
                    subset_file_name_substr = "mrmr_perc50",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 263,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 267,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr_perc50")),
                    subset_file_name_substr = "mrmr_perc50",
                    create_all_common = FALSE)
create_data_subsets(dparg_id = 267,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE)
###########
###########

#adult combat+log_cpm

explore_common_features(dparg_id = 271,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("ga_rf", "mrmr_perc50", "ranger_pos_impu_cor", "RF_RFE", "mrmr100"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_adult_combat_log_cpm",
                        dir_path = "../plots/fem_pipeline_results_adult_combat_log_cpm/common_features_upset")
explore_common_features(dparg_id = 275,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr10", "wilcoxontest", "t-test", "RF_RFE", "mrmr75", "mrmr100"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_adult_combat_log_cpm",
                        dir_path = "../plots/fem_pipeline_results_adult_combat_log_cpm/common_features_upset")
explore_common_features(dparg_id = 279,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("ga_rf", "mrmr_perc50", "ranger_pos_impu_cor", "mrmr100"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_adult_combat_log_cpm",
                        dir_path = "../plots/fem_pipeline_results_adult_combat_log_cpm/common_features_upset")

create_data_subsets(dparg_id = 271,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr_perc50")),
                    subset_file_name_substr = "mrmr_perc50",
                    create_all_common = FALSE, combat = TRUE)
create_data_subsets(dparg_id = 271,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE, combat = TRUE)

create_data_subsets(dparg_id = 275,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE, combat = TRUE)

create_data_subsets(dparg_id = 279,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr_perc50")),
                    subset_file_name_substr = "mrmr_perc50",
                    create_all_common = FALSE, combat = TRUE)
create_data_subsets(dparg_id = 279,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE, combat = TRUE)

###########
###########

#adult combat+log_tmm

explore_common_features(dparg_id = 283,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("ga_rf", "t-test", "mrmr_perc50", "ranger_pos_impu_cor", "mrmr100"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_adult_combat_log_tmm",
                        dir_path = "../plots/fem_pipeline_results_adult_combat_log_tmm/common_features_upset")
explore_common_features(dparg_id = 287,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("mrmr10", "wilcoxontest", "RF_RFE", "t-test", "ga_rf",
                                         "mrmr_perc50", "ranger_pos_impu_cor", "mrmr75"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_adult_combat_log_tmm",
                        dir_path = "../plots/fem_pipeline_results_adult_combat_log_tmm/common_features_upset")
explore_common_features(dparg_id = 291,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("ga_rf", "mrmr_perc50", "ranger_pos_impu_cor", "RF_RFE",
                                         "mrmr10", "mrmr100"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_adult_combat_log_tmm",
                        dir_path = "../plots/fem_pipeline_results_adult_combat_log_tmm/common_features_upset")

create_data_subsets(dparg_id = 283,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr_perc50")),
                    subset_file_name_substr = "mrmr_perc50",
                    create_all_common = FALSE, combat = TRUE)
create_data_subsets(dparg_id = 283,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE, combat = TRUE)

create_data_subsets(dparg_id = 287,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr_perc50")),
                    subset_file_name_substr = "mrmr_perc50",
                    create_all_common = FALSE, combat = TRUE)
create_data_subsets(dparg_id = 287,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE, combat = TRUE)

create_data_subsets(dparg_id = 291,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr_perc50")),
                    subset_file_name_substr = "mrmr_perc50",
                    create_all_common = FALSE, combat = TRUE)
create_data_subsets(dparg_id = 291,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE, combat = TRUE)


###########
###########

#adult f508del__f508del mutation

explore_common_features(dparg_id = 350,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("ga_rf", "t-test", "wilcoxontest",
                                         "mrmr_perc50", "ranger_pos_impu_cor", "mrmr100"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_adult_log_cpm_F508del___F508del",
                        dir_path = "../plots/fem_pipeline_results_adult_log_cpm_F508del___F508del/common_features_upset")
explore_common_features(dparg_id = 354,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("RF_RFE", "t-test", "wilcoxontest", "mrmr10",
                                         "mrmr100", "ga_rf", "mrmr75"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_adult_log_cpm_F508del___F508del",
                        dir_path = "../plots/fem_pipeline_results_adult_log_cpm_F508del___F508del/common_features_upset")
explore_common_features(dparg_id = 358,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("ga_rf", "ranger_pos_impu_cor", "mrmr_perc50",
                                         "mrmr75"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_adult_log_cpm_F508del___F508del",
                        dir_path = "../plots/fem_pipeline_results_adult_log_cpm_F508del___F508del/common_features_upset")



create_data_subsets(dparg_id = 350,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr_perc50")),
                    subset_file_name_substr = "mrmr_perc50",
                    create_all_common = FALSE, combat = FALSE)
create_data_subsets(dparg_id = 350,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE, combat = FALSE)

create_data_subsets(dparg_id = 354,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr100")),
                    subset_file_name_substr = "mrmr100",
                    create_all_common = FALSE, combat = FALSE)

create_data_subsets(dparg_id = 358,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr_perc50")),
                    subset_file_name_substr = "mrmr_perc50",
                    create_all_common = FALSE, combat = FALSE)
create_data_subsets(dparg_id = 358,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("mrmr75")),
                    subset_file_name_substr = "mrmr75",
                    create_all_common = FALSE, combat = FALSE)



###########
###########

#adult seurat normalize and find varible features no other norm

explore_common_features(dparg_id = 384,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("ga_rf", "wilcoxontest", "ranger_pos_impu_cor",
                                         "t-test", "RF_RFE",
                                         "mrmr_perc50", "mrmr75"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_adult_seurat3_norm_find_var_none",
                        dir_path = "../plots/fem_pipeline_results_adult_seurat3_norm_find_var_none/common_features_upset")
explore_common_features(dparg_id = 388,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("wilcoxontest", "RF_RFE", "t-test", "ga_rf",
                                         "ranger_pos_impu_cor", "mrmr100"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_adult_seurat3_norm_find_var_none",
                        dir_path = "../plots/fem_pipeline_results_adult_seurat3_norm_find_var_none/common_features_upset")
explore_common_features(dparg_id = 392,
                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                        best_fsm_vec = c("ga_rf", "ranger_pos_impu_cor", "RF_RFE", 
                                         "t-test", "wilcoxontest", "mrmr75"),
                        min_iter_feature_presence = 28,
                        results_dir = "../fem_pipeline_results_adult_seurat3_norm_find_var_none",
                        dir_path = "../plots/fem_pipeline_results_adult_seurat3_norm_find_var_none/common_features_upset")

create_data_subsets(dparg_id = 384,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("wilcoxontest")),
                    subset_file_name_substr = "wilcoxontest",
                    create_all_common = FALSE,
                    data_file_path = "../data/formatted/umi_counts_seurat3_with_norm_and_find_var_feat.csv")
create_data_subsets(dparg_id = 388,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("wilcoxontest")),
                    subset_file_name_substr = "wilcoxontest",
                    create_all_common = FALSE,
                    data_file_path = "../data/formatted/umi_counts_seurat3_with_norm_and_find_var_feat.csv")
create_data_subsets(dparg_id = 392,
                    dataset_pipeline_arguments = dataset_pipeline_arguments,
                    min_iter_feature_presence = 28,
                    subset_creation_criteria <- list("i"= c("ranger_pos_impu_cor")),
                    subset_file_name_substr = "ranger_pos_impu_cor",
                    create_all_common = FALSE,
                    data_file_path = "../data/formatted/umi_counts_seurat3_with_norm_and_find_var_feat.csv")
