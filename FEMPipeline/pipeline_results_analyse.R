setwd("~/UNSW/VafaeeLab/CysticFibrosisGroup/ExoCF/CFRD_EV_biomarker/FEMPipeline/")
library(tidyverse)
library(viridis)
library(ComplexHeatmap)
source("dataset_pipeline_arguments.R")
source("utils.R")

dparg_vec = c(199, 203, 207)
dataset_pipeline_arguments = dataset_pipeline_arguments
results_dir = "../fem_pipeline_results_AU_prefiltered_adult_combat_logtmm"
dir_path = "../plots/fem_pipeline_results_AU_prefiltered_adult_combat_logtmm/"
dataset_replace_string = "CF_EV_"
plot_heatmap <- function(dparg_vec,
                         dataset_pipeline_arguments = dataset_pipeline_arguments,
                         results_dir = "fem_pipeline_results",
                         dir_path = "",
                         dataset_replace_string = ""){
  data_info <- read.table(paste(results_dir, "data_info.csv", sep = "/"), 
                          sep = ',', header = TRUE)
  fsm_info <- read.table(paste(results_dir, "fsm_info.csv", sep = "/"),
                         sep = ',', header = TRUE)
  model_results <- read.table(paste(results_dir, "model_results_test.csv", sep = "/"),
                              sep = ',', header = TRUE)
  
  if(!dir.exists(dir_path)){
    dir.create(dir_path, recursive = TRUE)
  }
  
  dataset_id_vec <- c()
  for(id in dparg_vec){
    dataset_id <- paste(dataset_pipeline_arguments[[id]]$dataset_id,
                        dataset_pipeline_arguments[[id]]$classification_criteria,
                        sep = "_")
    print(dataset_id)
    dataset_id_vec <- append(dataset_id_vec, dataset_id)
  }
  dataset_id_vec 
  
  replaced_dataset_id_vec <- gsub(dataset_replace_string, "", 
                                  dataset_id_vec, fixed = TRUE)
  
  model_results <- model_results %>%
    filter(DataSetId %in% dataset_id_vec) %>%
    mutate(DataSetId = gsub(dataset_replace_string, "", DataSetId, fixed = TRUE))
  
  # model <- "L2 Regularized logistic regression"
  for(model in unique(model_results$Model)){
    print(model)
    
    model_results_sub <- model_results %>%
      filter(Model == model)
    
    data_to_plot <- model_results_sub %>%
      select(DataSetId, FSM, Mean_AUC) %>%
      # pivot_wider(names_from = DataSetId, values_from = Mean_AUC, values_fn = length) %>%          
      pivot_wider(names_from = DataSetId, values_from = Mean_AUC) %>%
      column_to_rownames("FSM")
    data_to_plot <- data.matrix(data_to_plot)
    
    file_name_curr <- paste0("Mean_AUC",
                             gsub(" ", "_", model),
                             ".png")
    file_name_curr <- paste0(dir_path, file_name_curr)
    
    png(file_name_curr, units = "cm", width = 20, height = 15, res = 1200)
    ht <- Heatmap(data_to_plot, name = "Mean AUC",
                  col = magma(5),
                  rect_gp = gpar(col = "white", lwd = 1),
                  cluster_columns = FALSE,
                  column_order = replaced_dataset_id_vec[replaced_dataset_id_vec %in% colnames(data_to_plot)],
                  column_names_rot = 60,
                  column_names_gp = gpar(fontsize = 10),
                  row_title = "Feature Selection Methods",
                  row_names_side = "left",
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.text(sprintf("%.3f", data_to_plot[i, j]), x, y, gp = gpar(fontsize = 7, col = "slateblue3"))
                  })
    draw(ht, column_title = paste("Mean AUC results from", model))
    dev.off()
    
  }
  
  ds_heatmap_dir_path <- paste0(dir_path, "/ds/") 
  if(!dir.exists(ds_heatmap_dir_path)){
    dir.create(ds_heatmap_dir_path)
  }
  for(ds in replaced_dataset_id_vec){
    print(ds)
    
    model_results_sub <- model_results %>%
      filter(Model != "Simple logistic regression") %>%
      mutate(Model = gsub("Regularized logistic regression", "regu log_reg", Model)) %>%
      filter(DataSetId == ds)
    
    data_to_plot <- model_results_sub %>%
      select(Model, FSM, Mean_AUC) %>%
      # pivot_wider(names_from = DataSetId, values_from = Mean_AUC, values_fn = length) %>%          
      pivot_wider(names_from = Model, values_from = Mean_AUC) %>%
      column_to_rownames("FSM")
    data_to_plot <- data.matrix(data_to_plot)
    
    file_name_curr <- paste0("Mean_AUC",
                             gsub(" ", "_", ds),
                             ".png")
    
    file_name_curr <- paste0(ds_heatmap_dir_path, file_name_curr)
    
    png(file_name_curr, units = "cm", width = 20, height = 15, res = 1200)
    ht <- Heatmap(data_to_plot, name = "Mean AUC",
                  col = magma(5),
                  rect_gp = gpar(col = "white", lwd = 1),
                  column_names_rot = 60,
                  column_names_gp = gpar(fontsize = 10),
                  row_title = "Feature Selection Methods",
                  row_names_side = "left",
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.text(sprintf("%.3f", data_to_plot[i, j]), x, y, gp = gpar(fontsize = 7, col = "slateblue3"))
                  })
    draw(ht, column_title = paste("Mean AUC results for", ds))
    dev.off()
    
  }
  
  
  dir_path <- paste0(dir_path, "common_barplot/")
  if(!dir.exists(dir_path)){
    dir.create(dir_path)
  }
  classification_models <- unique(model_results$Model)
  for (cm in classification_models) {
    individual_model <- model_results %>%
      filter(Model == cm)
    model_barplot <- ggplot(individual_model, aes(x=DataSetId, fill=FSM, y=Mean_AUC)) +
      geom_bar(stat="identity", position="dodge") +
      geom_errorbar( aes(x=DataSetId, ymin=X95.CI_AUC_lower, ymax=X95.CI_AUC_upper), position="dodge") +
      scale_fill_viridis_d() +
      theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
            axis.text.y = element_text(size=rel(1.2), face="italic", hjust=0.95),
            strip.text = element_text(size=rel(1.2), face="bold"),
            legend.title = element_text(size=rel(1.2)),
            legend.text = element_text(size=rel(1.1))
      ) +
      facet_wrap(facets = vars(Model))  
    
    file_name <- paste0(dir_path, str_replace(cm, " ", ""), "_barplot.png")
    ggsave(file_name, model_barplot, width=20, height=10, dpi=500)
  }  
}

plot_features_count_barplot <- function(dparg_vec,
                                        dataset_pipeline_arguments = dataset_pipeline_arguments,
                                        dir_path = "", results_dir = "fem_pipeline_results") {
  
  fsm_info <- read.table(paste(results_dir, "fsm_info.csv", sep = "/"),
                         sep = ',', header = TRUE)
  
  dataset_id_vec <- c()
  for(id in dparg_vec){
    dataset_id <- paste(dataset_pipeline_arguments[[id]]$dataset_id,
                        dataset_pipeline_arguments[[id]]$classification_criteria,
                        sep = "_")
    print(dataset_id)
    dataset_id_vec <- append(dataset_id_vec, dataset_id)
  }
  dataset_id_vec 
  
  fsm_info <- fsm_info %>%
    filter(DataSetId %in% dataset_id_vec)
  
  if (dir_path != "" & !dir.exists(dir_path)){
    dir.create(dir_path, recursive = TRUE)
  }
  
  features_barplot <- ggplot(fsm_info, aes(x=DataSetId, fill=FSM, y=Mean_Number.of.features)) +
    geom_bar(stat="identity", position="dodge") +
    geom_errorbar( aes(x=DataSetId, ymin=X95.CI_Number.of.features_lower, ymax=X95.CI_Number.of.features_upper), position="dodge") +
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
          axis.text.y = element_text(size=rel(1.2), face="italic", hjust=0.95),
          axis.title.x = element_text(size=rel(1.5)),
          axis.title.y = element_text(size=rel(1.5), angle=90),
          legend.title = element_text(size=rel(1.5)),
          legend.text = element_text(size=rel(1.2))) +
    scale_y_log10() +
    scale_fill_viridis_d() +
    ylab("Mean number of features used")
  ggsave(append_path(dir_path, "features_count_barplot.png"),
         features_barplot, width=12, height=12, dpi=500)
}

##########################


library(tidyverse)
library(viridis)
library(ComplexHeatmap)
source("scripts/R/utils.R")



# dparg_vec <- c(31:35)
# criteria <- 28

dparg_vec <- c(163:165)
criteria <- 29
results_dir = "fem_pipeline_results_subset"
dataset_replace_string = "transcriptomic_simple_norm_PREOPEVsMET_"
heatmap_file_name = "tr_PREOPEVsMET.png"

plot_common_feature_heatmap(c(163:165),
                            results_dir = "fem_pipeline_results_subset",
                            dataset_replace_string = "transcriptomic_simple_norm_PREOPEVsMET_",
                            heatmap_file_name = "tr_PREOPEVsMET.png"
)

plot_common_feature_heatmap <- function(dparg_vec, 
                                        results_dir = "fem_pipeline_results",
                                        dataset_replace_string = "",
                                        heatmap_file_name,
                                        plot_dir_path = "../plots/FEMPipeline_tmm_all/subset/"){
  data_info <- read.table(paste(results_dir, "data_info.csv", sep = "/"), 
                          sep = ',', header = TRUE)
  fsm_info <- read.table(paste(results_dir, "fsm_info.csv", sep = "/"),
                         sep = ',', header = TRUE)
  model_results <- read.table(paste(results_dir, "model_results_test.csv", sep = "/"),
                              sep = ',', header = TRUE)
  
  
  dataset_id_vec <- c()
  classification_criteria <- "" 
  for(id in dparg_vec){
    classification_criteria <- dataset_pipeline_arguments[[id]]$classification_criteria
    dataset_id <- paste(dataset_pipeline_arguments[[id]]$dataset_id,
                        classification_criteria,
                        sep = "_")
    print(dataset_id)
    dataset_id_vec <- append(dataset_id_vec, dataset_id)
  }
  dataset_id_vec 
  
  model_results <- model_results %>%
    filter(DataSetId %in% dataset_id_vec) %>%
    mutate(DataSetId = gsub(dataset_replace_string, "", DataSetId, fixed = TRUE)) %>%
    mutate(DataSetId = gsub(classification_criteria, "", DataSetId)) %>%
    mutate(DataSetId = gsub("^_|_$", "", DataSetId))
  
  data_to_plot <- model_results %>%
    select(DataSetId, Model, Mean_AUC) %>%
    pivot_wider(names_from = DataSetId, values_from = Mean_AUC) %>%
    column_to_rownames("Model")
  data_to_plot <- data.matrix(data_to_plot)
  
  if(!dir.exists(plot_dir_path)){
    dir.create(plot_dir_path, recursive = TRUE)
  }
  file_path <- paste0(plot_dir_path, heatmap_file_name)
  
  #classification_criteria obtained from dp_arg is used.
  #assuming that all classification criterias within the heatmap to be created are the same
  heatmap_title <- paste("Mean AUC results from selected common features for", classification_criteria)
  
  png(file_path, units = "cm", width = 30, height = 25, res = 1200)
  ht <- Heatmap(data_to_plot, name = "Mean AUC",
                col = magma(5),
                rect_gp = gpar(col = "white", lwd = 1),
                column_names_rot = 60,
                column_names_gp = gpar(fontsize = 10),
                row_title = "Classification Models",
                row_names_side = "left",
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.4f", data_to_plot[i, j]), x, y, gp = gpar(fontsize = 10, col = "slateblue3"))
                })
  draw(ht, column_title = heatmap_title)
  dev.off()
}



####################### initial cohort plots on newly quantified results
#### plot_heatmap(dparg_vec = c(41:43, 47:49, 53:55),
plot_heatmap(
  dparg_vec = c(1, 5, 9),
  dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
  results_dir = "fem_pipeline_results_tr",
  dir_path = "plots/FEMPipeline_new_quant/heatmap/transcriptomic/",
  dataset_replace_string = "_initial"
)

plot_JI_heatmap(filename = "all_ji.csv",
                dir = "fem_pipeline_results_tr/JI/",
                dataset_pipeline_arguments = dataset_pipeline_arguments_transcriptomic,
                heatmapfilename = "tr_JI.svg",
                heatmapfiledir = "plots/FEMPipeline_new_quant/",
                fsm_allowed = fsm_vector,
                dparg_vec = c(1, 5, 9),
                dataset_replace_string = "_initial")


##########################



plot_heatmap(
  dparg_vec = c(1, 5, 9),
  dataset_pipeline_arguments = dataset_pipeline_arguments,
  results_dir = "../fem_pipeline_results",
  dir_path = "../plots/FEMPipeline/",
  dataset_replace_string = "CF_EV_"
)

plot_heatmap(
  dparg_vec = c(13, 17, 21),
  dataset_pipeline_arguments = dataset_pipeline_arguments,
  results_dir = "../fem_pipeline_results",
  dir_path = "../plots/FEMPipeline_tmm_all/",
  dataset_replace_string = "CF_EV_"
)

plot_heatmap(
  dparg_vec = c(25, 29, 33),
  dataset_pipeline_arguments = dataset_pipeline_arguments,
  results_dir = "../fem_pipeline_results_AU",
  dir_path = "../plots/FEMPipeline_AU/",
  dataset_replace_string = "CF_EV_"
)

plot_heatmap(
  dparg_vec = c(37, 41, 45),
  dataset_pipeline_arguments = dataset_pipeline_arguments,
  results_dir = "../fem_pipeline_results_AU",
  dir_path = "../plots/FEMPipeline_AU_tmm/",
  dataset_replace_string = "CF_EV_"
)


plot_features_count_barplot(dparg_vec = c(37, 41, 45),
                            dataset_pipeline_arguments = dataset_pipeline_arguments,
                            dir_path = "../plots/FEMPipeline_AU_tmm/",
                            results_dir = "../fem_pipeline_results_AU")



plot_common_feature_heatmap(c(49:54),
                            results_dir = "../fem_pipeline_results_AU_subset",
                            dataset_replace_string = "CF_EV_AU_zlogtmm_",
                            heatmap_file_name = "AU_tmm_CFRDVsIGT.png"
)

plot_common_feature_heatmap(c(55:56),
                            results_dir = "../fem_pipeline_results_AU_subset",
                            dataset_replace_string = "CF_EV_AU_zlogtmm_",
                            heatmap_file_name = "AU_tmm_CFRDVsNGT.png"
)

plot_common_feature_heatmap(c(57:58),
                            results_dir = "../fem_pipeline_results_AU_subset",
                            dataset_replace_string = "CF_EV_AU_zlogtmm_",
                            heatmap_file_name = "AU_tmm_IGTVsNGT.png"
)


plot_common_feature_heatmap(c(59, 60),
                            results_dir = "../fem_pipeline_results_subset",
                            dataset_replace_string = "CF_EV_zlogtmm_",
                            heatmap_file_name = "tmm_CFRDVsIGT.png"
)

plot_common_feature_heatmap(c(61, 62),
                            results_dir = "../fem_pipeline_results_subset",
                            dataset_replace_string = "CF_EV_zlogtmm_",
                            heatmap_file_name = "tmm_CFRDVsNGT.png"
)

plot_common_feature_heatmap(c(63, 64),
                            results_dir = "../fem_pipeline_results_subset",
                            dataset_replace_string = "CF_EV_zlogtmm_",
                            heatmap_file_name = "tmm_IGTVsNGT.png"
)



plot_heatmap(
  dparg_vec = c(65, 69, 73),
  dataset_pipeline_arguments = dataset_pipeline_arguments,
  results_dir = "../fem_pipeline_results_AU_logtmm",
  dir_path = "../plots/FEMPipeline_AU_noscaling_tmm/",
  dataset_replace_string = "CF_EV_"
)

plot_common_feature_heatmap(c(77, 78),
                            results_dir = "../fem_pipeline_results_AU_logtmm_subset",
                            dataset_replace_string = "CF_EV_AU_zlogtmm_",
                            heatmap_file_name = "AU_noscale_tmm_CFRDVsIGT.png",
                            plot_dir_path = "../plots/FEMPipeline_AU_noscaling_tmm/subset/"
)

plot_common_feature_heatmap(c(79, 80),
                            results_dir = "../fem_pipeline_results_AU_logtmm_subset",
                            dataset_replace_string = "CF_EV_AU_zlogtmm_",
                            heatmap_file_name = "AU_noscale_tmm_CFRDVsNGT.png",
                            plot_dir_path = "../plots/FEMPipeline_AU_noscaling_tmm/subset/"
)

plot_common_feature_heatmap(c(81, 82),
                            results_dir = "../fem_pipeline_results_AU_logtmm_subset",
                            dataset_replace_string = "CF_EV_AU_zlogtmm_",
                            heatmap_file_name = "AU_noscale_tmm_IGTVsNGT.png",
                            plot_dir_path = "../plots/FEMPipeline_AU_noscaling_tmm/subset/"
)




plot_heatmap(
  dparg_vec = c(83, 87, 91),
  dataset_pipeline_arguments = dataset_pipeline_arguments,
  results_dir = "../fem_pipeline_results_AU_adult_logtmm",
  dir_path = "../plots/fem_pipeline_results_AU_adult_logtmm/",
  dataset_replace_string = "CF_EV_"
)


plot_common_feature_heatmap(c(95, 96),
                            results_dir = "../fem_pipeline_results_AU_adult_logtmm_subset",
                            dataset_replace_string = "CF_EV_AU_",
                            heatmap_file_name = "AU_noscale_tmm_CFRDVsIGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_adult_logtmm/subset/"
)

plot_common_feature_heatmap(c(97, 98),
                            results_dir = "../fem_pipeline_results_AU_adult_logtmm_subset",
                            dataset_replace_string = "CF_EV_AU_",
                            heatmap_file_name = "AU_noscale_tmm_CFRDVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_adult_logtmm/subset/"
)

plot_common_feature_heatmap(c(99, 100),
                            results_dir = "../fem_pipeline_results_AU_adult_logtmm_subset",
                            dataset_replace_string = "CF_EV_AU_",
                            heatmap_file_name = "AU_noscale_tmm_IGTVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_adult_logtmm/subset/"
)

plot_common_feature_heatmap(c(99, 100, 101, 102),
                            results_dir = "../fem_pipeline_results_AU_adult_logtmm_subset",
                            dataset_replace_string = "CF_EV_AU_",
                            heatmap_file_name = "AU_noscale_tmm_IGTVsNGT_all.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_adult_logtmm/subset/"
)



##########

#train on DK adult

plot_heatmap(
  dparg_vec = c(103, 107, 111),
  dataset_pipeline_arguments = dataset_pipeline_arguments,
  results_dir = "../fem_pipeline_results_DK_adult_logtmm",
  dir_path = "../plots/fem_pipeline_results_DK_adult_logtmm/",
  dataset_replace_string = "CF_EV_"
)

plot_common_feature_heatmap(c(115, 116),
                            results_dir = "../fem_pipeline_results_DK_adult_logtmm_subset",
                            dataset_replace_string = "CF_EV_DK_",
                            heatmap_file_name = "DK_noscale_tmm_CFRDVsIGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_DK_adult_logtmm/subset/"
)

plot_common_feature_heatmap(c(117, 118),
                            results_dir = "../fem_pipeline_results_DK_adult_logtmm_subset",
                            dataset_replace_string = "CF_EV_DK_",
                            heatmap_file_name = "DK_noscale_tmm_CFRDVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_DK_adult_logtmm/subset/"
)

plot_common_feature_heatmap(c(119, 120),
                            results_dir = "../fem_pipeline_results_DK_adult_logtmm_subset",
                            dataset_replace_string = "CF_EV_DK_",
                            heatmap_file_name = "DK_noscale_tmm_IGTVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_DK_adult_logtmm/subset/"
)


##########

#AU adult logcpm

plot_heatmap(
  dparg_vec = c(121, 125, 129),
  dataset_pipeline_arguments = dataset_pipeline_arguments,
  results_dir = "../fem_pipeline_results_AU_adult_logcpm",
  dir_path = "../plots/fem_pipeline_results_AU_adult_logcpm/",
  dataset_replace_string = "CF_EV_"
)

plot_common_feature_heatmap(c(169, 170),
                            results_dir = "../fem_pipeline_results_AU_adult_logcpm_subset",
                            dataset_replace_string = "CF_EV_AU_",
                            heatmap_file_name = "AU_logcpm_CFRDVsIGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_adult_logcpm/subset/"
)

plot_common_feature_heatmap(c(171, 172),
                            results_dir = "../fem_pipeline_results_AU_adult_logcpm_subset",
                            dataset_replace_string = "CF_EV_AU_",
                            heatmap_file_name = "AU_logcpm_CFRDVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_adult_logcpm/subset/"
)

plot_common_feature_heatmap(c(173, 174),
                            results_dir = "../fem_pipeline_results_AU_adult_logcpm_subset",
                            dataset_replace_string = "CF_EV_AU_",
                            heatmap_file_name = "AU_logcpm_IGTVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_adult_logcpm/subset/"
)


##########

#AU adult log

plot_heatmap(
  dparg_vec = c(133, 137, 141),
  dataset_pipeline_arguments = dataset_pipeline_arguments,
  results_dir = "../fem_pipeline_results_AU_adult_log",
  dir_path = "../plots/fem_pipeline_results_AU_adult_log/",
  dataset_replace_string = "CF_EV_"
)


plot_common_feature_heatmap(c(175, 176),
                            results_dir = "../fem_pipeline_results_AU_adult_log_subset",
                            dataset_replace_string = "CF_EV_AU_",
                            heatmap_file_name = "AU_log_CFRDVsIGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_adult_log/subset/"
)

plot_common_feature_heatmap(c(177, 178),
                            results_dir = "../fem_pipeline_results_AU_adult_log_subset",
                            dataset_replace_string = "CF_EV_AU_",
                            heatmap_file_name = "AU_log_CFRDVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_adult_log/subset/"
)

plot_common_feature_heatmap(c(179, 180, 181, 182),
                            results_dir = "../fem_pipeline_results_AU_adult_log_subset",
                            dataset_replace_string = "CF_EV_AU_",
                            heatmap_file_name = "AU_log_IGTVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_adult_log/subset/"
)

##########

#AU adult combatseq logcpm

plot_heatmap(
  dparg_vec = c(145, 149, 153),
  dataset_pipeline_arguments = dataset_pipeline_arguments,
  results_dir = "../fem_pipeline_results_AU_adult_combat_logcpm",
  dir_path = "../plots/fem_pipeline_results_AU_adult_combat_logcpm/",
  dataset_replace_string = "CF_EV_"
)

plot_common_feature_heatmap(c(183, 184, 185, 186),
                            results_dir = "../fem_pipeline_results_AU_adult_combat_logcpm_subset",
                            dataset_replace_string = "CF_EV_AU_",
                            heatmap_file_name = "AU_logcpm_CFRDVsIGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_adult_combat_logcpm/subset/"
)

plot_common_feature_heatmap(c(187, 188),
                            results_dir = "../fem_pipeline_results_AU_adult_combat_logcpm_subset",
                            dataset_replace_string = "CF_EV_AU_",
                            heatmap_file_name = "AU_logcpm_CFRDVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_adult_combat_logcpm/subset/"
)

plot_common_feature_heatmap(c(189, 190),
                            results_dir = "../fem_pipeline_results_AU_adult_combat_logcpm_subset",
                            dataset_replace_string = "CF_EV_AU_",
                            heatmap_file_name = "AU_logcpm_IGTVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_adult_combat_logcpm/subset/"
)


##########

#AU adult combatseq logtmm

plot_heatmap(
  dparg_vec = c(157, 161, 165),
  dataset_pipeline_arguments = dataset_pipeline_arguments,
  results_dir = "../fem_pipeline_results_AU_adult_combat_logtmm",
  dir_path = "../plots/fem_pipeline_results_AU_adult_combat_logtmm/",
  dataset_replace_string = "CF_EV_"
)

plot_common_feature_heatmap(c(191, 192, 193, 194),
                            results_dir = "../fem_pipeline_results_AU_adult_combat_logtmm_subset",
                            dataset_replace_string = "CF_EV_AU_",
                            heatmap_file_name = "AU_logtmm_CFRDVsIGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_adult_combat_logtmm/subset/"
)

plot_common_feature_heatmap(c(195, 196),
                            results_dir = "../fem_pipeline_results_AU_adult_combat_logtmm_subset",
                            dataset_replace_string = "CF_EV_AU_",
                            heatmap_file_name = "AU_logtmm_CFRDVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_adult_combat_logtmm/subset/"
)

plot_common_feature_heatmap(c(197, 198),
                            results_dir = "../fem_pipeline_results_AU_adult_combat_logtmm_subset",
                            dataset_replace_string = "CF_EV_AU_",
                            heatmap_file_name = "AU_logtmm_IGTVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_adult_combat_logtmm/subset/"
)

##################################################

#AU adult prefiltered combatseq logtmm

plot_heatmap(
  dparg_vec = c(199, 203, 207),
  dataset_pipeline_arguments = dataset_pipeline_arguments,
  results_dir = "../fem_pipeline_results_AU_prefiltered_adult_combat_logtmm",
  dir_path = "../plots/fem_pipeline_results_AU_prefiltered_adult_combat_logtmm/",
  dataset_replace_string = "CF_EV_"
)

plot_common_feature_heatmap(c(223, 224),
                            results_dir = "../fem_pipeline_results_AU_prefiltered_adult_combat_logtmm_subset",
                            dataset_replace_string = "CF_EV_AU_prefiltered_adult_combat",
                            heatmap_file_name = "AU_logtmm_CFRDVsIGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_prefiltered_adult_combat_logtmm/subset/"
)
plot_common_feature_heatmap(c(225, 226),
                            results_dir = "../fem_pipeline_results_AU_prefiltered_adult_combat_logtmm_subset",
                            dataset_replace_string = "CF_EV_AU_prefiltered_adult_combat",
                            heatmap_file_name = "AU_logtmm_CFRDVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_prefiltered_adult_combat_logtmm/subset/"
)
plot_common_feature_heatmap(c(227, 228),
                            results_dir = "../fem_pipeline_results_AU_prefiltered_adult_combat_logtmm_subset",
                            dataset_replace_string = "CF_EV_AU_prefiltered_adult_combat",
                            heatmap_file_name = "AU_logtmm_IGTVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_prefiltered_adult_combat_logtmm/subset/"
)


#AU adult prefiltered combatseq logcpm

plot_heatmap(
  dparg_vec = c(211, 215, 219),
  dataset_pipeline_arguments = dataset_pipeline_arguments,
  results_dir = "../fem_pipeline_results_AU_prefiltered_adult_combat_logcpm",
  dir_path = "../plots/fem_pipeline_results_AU_prefiltered_adult_combat_logcpm/",
  dataset_replace_string = "CF_EV_"
)

plot_common_feature_heatmap(c(229, 230),
                            results_dir = "../fem_pipeline_results_AU_prefiltered_adult_combat_logcpm_subset",
                            dataset_replace_string = "CF_EV_AU_prefiltered_adult_combat",
                            heatmap_file_name = "AU_logcpm_CFRDVsIGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_prefiltered_adult_combat_logcpm/subset/"
)
plot_common_feature_heatmap(c(231, 232),
                            results_dir = "../fem_pipeline_results_AU_prefiltered_adult_combat_logcpm_subset",
                            dataset_replace_string = "CF_EV_AU_prefiltered_adult_combat",
                            heatmap_file_name = "AU_logcpm_CFRDVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_prefiltered_adult_combat_logcpm/subset/"
)
plot_common_feature_heatmap(c(233, 234),
                            results_dir = "../fem_pipeline_results_AU_prefiltered_adult_combat_logcpm_subset",
                            dataset_replace_string = "CF_EV_AU_prefiltered_adult_combat",
                            heatmap_file_name = "AU_logcpm_IGTVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_prefiltered_adult_combat_logcpm/subset/"
)

##########

#adult log

plot_heatmap(
  dparg_vec = c(235, 239, 243),
  dataset_pipeline_arguments = dataset_pipeline_arguments,
  results_dir = "../fem_pipeline_results_adult_log",
  dir_path = "../plots/fem_pipeline_results_adult_log/",
  dataset_replace_string = "CF_EV_"
)

plot_common_feature_heatmap(c(295, 296, 297, 298),
                            results_dir = "../fem_pipeline_results_adult_log_subset",
                            dataset_replace_string = "CF_EV_",
                            heatmap_file_name = "log_CFRDVsIGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_adult_log/subset/"
)
plot_common_feature_heatmap(c(299, 300),
                            results_dir = "../fem_pipeline_results_adult_log_subset",
                            dataset_replace_string = "CF_EV_",
                            heatmap_file_name = "log_CFRDVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_adult_log/subset/"
)
plot_common_feature_heatmap(c(301, 302, 303, 304),
                            results_dir = "../fem_pipeline_results_adult_log_subset",
                            dataset_replace_string = "CF_EV_",
                            heatmap_file_name = "log_IGTVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_adult_log/subset/"
)

##########

#adult log_cpm

plot_heatmap(
  dparg_vec = c(247, 251, 255),
  dataset_pipeline_arguments = dataset_pipeline_arguments,
  results_dir = "../fem_pipeline_results_adult_log_cpm",
  dir_path = "../plots/fem_pipeline_results_adult_log_cpm/",
  dataset_replace_string = "CF_EV_"
)

plot_common_feature_heatmap(c(305, 306),
                            results_dir = "../fem_pipeline_results_adult_log_cpm_subset",
                            dataset_replace_string = "CF_EV_",
                            heatmap_file_name = "log_cpm_CFRDVsIGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_adult_log_cpm/subset/"
)
plot_common_feature_heatmap(c(307, 308),
                            results_dir = "../fem_pipeline_results_adult_log_cpm_subset",
                            dataset_replace_string = "CF_EV_",
                            heatmap_file_name = "log_cpm_CFRDVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_adult_log_cpm/subset/"
)
plot_common_feature_heatmap(c(309, 310, 311, 312),
                            results_dir = "../fem_pipeline_results_adult_log_cpm_subset",
                            dataset_replace_string = "CF_EV_",
                            heatmap_file_name = "log_cpm_IGTVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_adult_log_cpm/subset/"
)

##########

#adult log_tmm

plot_heatmap(
  dparg_vec = c(259, 263, 267),
  dataset_pipeline_arguments = dataset_pipeline_arguments,
  results_dir = "../fem_pipeline_results_adult_log_tmm",
  dir_path = "../plots/fem_pipeline_results_adult_log_tmm/",
  dataset_replace_string = "CF_EV_"
)

plot_common_feature_heatmap(c(313, 314, 315, 316),
                            results_dir = "../fem_pipeline_results_adult_log_tmm_subset",
                            dataset_replace_string = "CF_EV_",
                            heatmap_file_name = "log_tmm_CFRDVsIGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_adult_log_tmm/subset/"
)
plot_common_feature_heatmap(c(317, 318, 319, 320),
                            results_dir = "../fem_pipeline_results_adult_log_tmm_subset",
                            dataset_replace_string = "CF_EV_",
                            heatmap_file_name = "log_tmm_CFRDVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_adult_log_tmm/subset/"
)
plot_common_feature_heatmap(c(321, 322, 323, 324),
                            results_dir = "../fem_pipeline_results_adult_log_tmm_subset",
                            dataset_replace_string = "CF_EV_",
                            heatmap_file_name = "log_tmm_IGTVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_adult_log_tmm/subset/"
)

##########

#adult combat+log_cpm

plot_heatmap(
  dparg_vec = c(271, 275, 279),
  dataset_pipeline_arguments = dataset_pipeline_arguments,
  results_dir = "../fem_pipeline_results_adult_combat_log_cpm",
  dir_path = "../plots/fem_pipeline_results_adult_combat_log_cpm/",
  dataset_replace_string = "CF_EV_"
)

plot_common_feature_heatmap(c(325, 326, 327, 328),
                            results_dir = "../fem_pipeline_results_adult_combat_log_cpm_subset",
                            dataset_replace_string = "CF_EV_",
                            heatmap_file_name = "combat_log_cpm_CFRDVsIGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_adult_combat_log_cpm/subset/"
)
plot_common_feature_heatmap(c(329, 330),
                            results_dir = "../fem_pipeline_results_adult_combat_log_cpm_subset",
                            dataset_replace_string = "CF_EV_",
                            heatmap_file_name = "combat_log_cpm_CFRDVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_adult_combat_log_cpm/subset/"
)
plot_common_feature_heatmap(c(331, 332, 333, 334),
                            results_dir = "../fem_pipeline_results_adult_combat_log_cpm_subset",
                            dataset_replace_string = "CF_EV_",
                            heatmap_file_name = "combat_log_cpm_IGTVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_adult_combat_log_cpm/subset/"
)
##########

#adult combat+log_tmm

plot_heatmap(
  dparg_vec = c(283, 287, 291),
  dataset_pipeline_arguments = dataset_pipeline_arguments,
  results_dir = "../fem_pipeline_results_adult_combat_log_tmm",
  dir_path = "../plots/fem_pipeline_results_adult_combat_log_tmm/",
  dataset_replace_string = "CF_EV_"
)

plot_common_feature_heatmap(c(335, 336, 337, 338),
                            results_dir = "../fem_pipeline_results_adult_combat_log_tmm_subset",
                            dataset_replace_string = "CF_EV_",
                            heatmap_file_name = "combat_log_tmm_CFRDVsIGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_adult_combat_log_tmm/subset/"
)
plot_common_feature_heatmap(c(339, 340, 341, 342),
                            results_dir = "../fem_pipeline_results_adult_combat_log_tmm_subset",
                            dataset_replace_string = "CF_EV_",
                            heatmap_file_name = "combat_log_tmm_CFRDVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_adult_combat_log_tmm/subset/"
)
plot_common_feature_heatmap(c(343, 344, 345, 346),
                            results_dir = "../fem_pipeline_results_adult_combat_log_tmm_subset",
                            dataset_replace_string = "CF_EV_",
                            heatmap_file_name = "combat_log_tmm_IGTVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_adult_combat_log_tmm/subset/"
)

##########

#adult plsda

plot_heatmap(
  dparg_vec = c(347, 348, 349),
  dataset_pipeline_arguments = dataset_pipeline_arguments,
  results_dir = "../fem_pipeline_results_adult_log_cpm_plsda",
  dir_path = "../plots/fem_pipeline_results_adult_log_cpm_plsda/",
  dataset_replace_string = "CF_EV_"
)


##########

#adult F508del__F508del 

plot_heatmap(
  dparg_vec = c(350, 354, 358),
  dataset_pipeline_arguments = dataset_pipeline_arguments,
  results_dir = "../fem_pipeline_results_adult_log_cpm_F508del___F508del",
  dir_path = "../plots/fem_pipeline_results_adult_log_cpm_F508del___F508del/",
  dataset_replace_string = "CF_EV_"
)

plot_common_feature_heatmap(c(362, 363, 364, 365),
                            results_dir = "../fem_pipeline_results_adult_log_cpm_F508del___F508del_subset",
                            dataset_replace_string = "CF_EV_",
                            heatmap_file_name = "log_cpm_CFRDVsIGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_adult_log_cpm_F508del___F508del/subset/"
)
plot_common_feature_heatmap(c(366, 367),
                            results_dir = "../fem_pipeline_results_adult_log_cpm_F508del___F508del_subset",
                            dataset_replace_string = "CF_EV_",
                            heatmap_file_name = "log_cpm_CFRDVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_adult_log_cpm_F508del___F508del/subset/"
)
plot_common_feature_heatmap(c(368, 369, 370, 371),
                            results_dir = "../fem_pipeline_results_adult_log_cpm_F508del___F508del_subset",
                            dataset_replace_string = "CF_EV_",
                            heatmap_file_name = "log_cpm_IGTVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_adult_log_cpm_F508del___F508del/subset/"
)


##########

#adult seurat normalize and find varible features logtmm
#pipeline execution not successful

plot_heatmap(
  dparg_vec = c(372, 376, 380),
  dataset_pipeline_arguments = dataset_pipeline_arguments,
  results_dir = "../fem_pipeline_results_adult_seurat3_norm_find_var_log_tmm",
  dir_path = "../plots/fem_pipeline_results_adult_seurat3_norm_find_var_log_tmm/",
  dataset_replace_string = "CF_EV_"
)

##########

#adult seurat normalize and find varible features no other norm

plot_heatmap(
  dparg_vec = c(384, 388, 392),
  dataset_pipeline_arguments = dataset_pipeline_arguments,
  results_dir = "../fem_pipeline_results_adult_seurat3_norm_find_var_none",
  dir_path = "../plots/fem_pipeline_results_adult_seurat3_norm_find_var_none/",
  dataset_replace_string = "CF_EV_"
)

plot_common_feature_heatmap(c(408, 409),
                            results_dir = "../fem_pipeline_results_adult_seurat3_norm_find_var_none_subset",
                            dataset_replace_string = "CF_EV_adult_seurat3_norm_find_var_none_",
                            heatmap_file_name = "seurat3_norm_find_var_none_CFRDVsIGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_adult_seurat3_norm_find_var_none/subset/"
)
plot_common_feature_heatmap(c(410, 411),
                            results_dir = "../fem_pipeline_results_adult_seurat3_norm_find_var_none_subset",
                            dataset_replace_string = "CF_EV_adult_seurat3_norm_find_var_none_",
                            heatmap_file_name = "seurat3_norm_find_var_none_CFRDVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_adult_seurat3_norm_find_var_none/subset/"
)
plot_common_feature_heatmap(c(412),
                            results_dir = "../fem_pipeline_results_adult_seurat3_norm_find_var_none_subset",
                            dataset_replace_string = "CF_EV_adult_seurat3_norm_find_var_none_",
                            heatmap_file_name = "seurat3_norm_find_var_none_IGTVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_adult_seurat3_norm_find_var_none/subset/"
)

##########

#adult seurat without normalize and find varible features + logtmm
#pipeline execution not successful

plot_heatmap(
  dparg_vec = c(396, 400, 404),
  dataset_pipeline_arguments = dataset_pipeline_arguments,
  results_dir = "../fem_pipeline_results_adult_seurat3_log_tmm",
  dir_path = "../plots/fem_pipeline_results_adult_seurat3_log_tmm/",
  dataset_replace_string = "CF_EV_"
)


##########

#AU filter + seurat3 norm + batchcorrectedDK + none norm 

plot_heatmap(
  dparg_vec = c(413, 417, 421),
  dataset_pipeline_arguments = dataset_pipeline_arguments,
  results_dir = "../fem_pipeline_results_AU_adult_filt_seurat3norm_none",
  dir_path = "../plots/fem_pipeline_results_AU_adult_filt_seurat3norm_none/",
  dataset_replace_string = "CF_EV_"
)

plot_common_feature_heatmap(c(437, 438),
                            results_dir = "../fem_pipeline_results_AU_adult_filt_seurat3norm_none_subset",
                            dataset_replace_string = "CF_EV_AU_adult_filt_seurat3norm_none_",
                            heatmap_file_name = "filt_seurat3norm_none_CFRDVsIGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_adult_filt_seurat3norm_none/subset/"
)
plot_common_feature_heatmap(c(439, 440),
                            results_dir = "../fem_pipeline_results_AU_adult_filt_seurat3norm_none_subset",
                            dataset_replace_string = "CF_EV_AU_adult_filt_seurat3norm_none_",
                            heatmap_file_name = "filt_seurat3norm_none_CFRDVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_adult_filt_seurat3norm_none/subset/"
)
plot_common_feature_heatmap(c(441, 442),
                            results_dir = "../fem_pipeline_results_AU_adult_filt_seurat3norm_none_subset",
                            dataset_replace_string = "CF_EV_AU_adult_filt_seurat3norm_none_",
                            heatmap_file_name = "filt_seurat3norm_none_IGTVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_adult_filt_seurat3norm_none/subset/"
)

##########

#AU seurat3 norm + batchcorrectedDK + none norm 

plot_heatmap(
  dparg_vec = c(425, 429, 433),
  dataset_pipeline_arguments = dataset_pipeline_arguments,
  results_dir = "../fem_pipeline_results_AU_adult_nofilt_seurat3norm_none",
  dir_path = "../plots/fem_pipeline_results_AU_adult_nofilt_seurat3norm_none/",
  dataset_replace_string = "CF_EV_"
)

plot_common_feature_heatmap(c(443, 444),
                            results_dir = "../fem_pipeline_results_AU_adult_nofilt_seurat3norm_none_subset",
                            dataset_replace_string = "CF_EV_AU_adult_nofilt_seurat3norm_none_",
                            heatmap_file_name = "nofilt_seurat3norm_none_CFRDVsIGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_adult_nofilt_seurat3norm_none/subset/"
)
plot_common_feature_heatmap(c(445, 446),
                            results_dir = "../fem_pipeline_results_AU_adult_nofilt_seurat3norm_none_subset",
                            dataset_replace_string = "CF_EV_AU_adult_nofilt_seurat3norm_none_",
                            heatmap_file_name = "nofilt_seurat3norm_none_CFRDVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_adult_nofilt_seurat3norm_none/subset/"
)
plot_common_feature_heatmap(c(447, 448),
                            results_dir = "../fem_pipeline_results_AU_adult_nofilt_seurat3norm_none_subset",
                            dataset_replace_string = "CF_EV_AU_adult_nofilt_seurat3norm_none_",
                            heatmap_file_name = "nofilt_seurat3norm_none_IGTVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_AU_adult_nofilt_seurat3norm_none/subset/"
)

##########

#adult filtered + seurat normalize and find varible features no other norm

plot_heatmap(
  dparg_vec = c(449, 453, 457),
  dataset_pipeline_arguments = dataset_pipeline_arguments,
  results_dir = "../fem_pipeline_results_adult_filtered_then_seurat3_norm_find_var_none",
  dir_path = "../plots/fem_pipeline_results_adult_filtered_then_seurat3_norm_find_var_none/",
  dataset_replace_string = "CF_EV_"
)

plot_common_feature_heatmap(c(408, 409),
                            results_dir = "../fem_pipeline_results_adult_filtered_then_seurat3_norm_find_var_none_subset",
                            dataset_replace_string = "CF_EV_adult_filtered__seurat3_norm_find_var_none_",
                            heatmap_file_name = "seurat3_norm_find_var_none_CFRDVsIGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_adult_filtered_then_seurat3_norm_find_var_none/subset/"
)
plot_common_feature_heatmap(c(410, 411),
                            results_dir = "../fem_pipeline_results_adult_filtered_then_seurat3_norm_find_var_none_subset",
                            dataset_replace_string = "CF_EV_adult_filtered__seurat3_norm_find_var_none_",
                            heatmap_file_name = "seurat3_norm_find_var_none_CFRDVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_adult_filtered_then_seurat3_norm_find_var_none/subset/"
)
plot_common_feature_heatmap(c(412),
                            results_dir = "../fem_pipeline_results_adult_filtered_then_seurat3_norm_find_var_none_subset",
                            dataset_replace_string = "CF_EV_adult_filtered__seurat3_norm_find_var_none_",
                            heatmap_file_name = "seurat3_norm_find_var_none_IGTVsNGT.png",
                            plot_dir_path = "../plots/fem_pipeline_results_adult_filtered_then_seurat3_norm_find_var_none/subset/"
)



##########