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
  #so special handling for features with "/"
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
                                create_all_common = TRUE,
                                dataset_pipeline_arguments = dataset_pipeline_arguments){
  
  ds <- dataset_pipeline_arguments[[dparg_id]]
  dataset_id <- paste(ds$dataset_id, ds$classification_criteria, sep = "_")
  print(dataset_id)
  
  output_dir <- "../data/selected_features/"
  file_name <- "best_FSM_common_features.csv"
  features_df <- read.csv(paste0(output_dir, file_name))
  
  features_df <- features_df %>%
    filter(datasetid == dataset_id) %>%
    filter(miniter_feature_presence == min_iter_feature_presence)
  
  
  ######create new data subsets
  
  data <- read.table("../data/formatted/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
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
