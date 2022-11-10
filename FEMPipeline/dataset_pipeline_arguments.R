dataset_pipeline_arguments <- list(
  #1
  #CFRDVsIGT
  list(phenotype_file_name = "data/formatted/phenotype.txt",
       read_count_dir_path = "data/formatted",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "CF_EV_zlogcpm",
       classification_criteria = "CFRDVsIGT",
       classes = c("IGT", "CFRD"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       norm = "norm_log_cpm_simple",
       fems_to_run = c("all", 
                       "t-test", "wilcoxontest",
                       "ranger_impu_cor",
                       "ranger_pos_impu_cor",
                       "mrmr10", "mrmr75")),
  
  #2
  #CFRDVsIGT
  list(phenotype_file_name = "data/formatted/phenotype.txt",
       read_count_dir_path = "data/formatted",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "CF_EV_zlogcpm",
       classification_criteria = "CFRDVsIGT",
       classes = c("IGT", "CFRD"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       norm = "norm_log_cpm_simple",
       fems_to_run = c("mrmr_perc50")),
  
  #3
  #CFRDVsIGT
  list(phenotype_file_name = "data/formatted/phenotype.txt",
       read_count_dir_path = "data/formatted",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "CF_EV_zlogcpm",
       classification_criteria = "CFRDVsIGT",
       classes = c("IGT", "CFRD"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       norm = "norm_log_cpm_simple",
       fems_to_run = c("RF_RFE")),
  
  #4
  #CFRDVsIGT
  list(phenotype_file_name = "data/formatted/phenotype.txt",
       read_count_dir_path = "data/formatted",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "CF_EV_zlogcpm",
       classification_criteria = "CFRDVsIGT",
       classes = c("IGT", "CFRD"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       norm = "norm_log_cpm_simple",
       fems_to_run = c("ga_rf")),  
  

  #1
  #CFRDVsNGT
  list(phenotype_file_name = "data/formatted/phenotype.txt",
       read_count_dir_path = "data/formatted",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "CF_EV_zlogcpm",
       classification_criteria = "CFRDVsNGT",
       classes = c("NGT", "CFRD"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       norm = "norm_log_cpm_simple",
       fems_to_run = c("all", 
                       "t-test", "wilcoxontest",
                       "ranger_impu_cor",
                       "ranger_pos_impu_cor",
                       "mrmr10", "mrmr75")),
  
  #2
  #CFRDVsNGT
  list(phenotype_file_name = "data/formatted/phenotype.txt",
       read_count_dir_path = "data/formatted",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "CF_EV_zlogcpm",
       classification_criteria = "CFRDVsNGT",
       classes = c("NGT", "CFRD"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       norm = "norm_log_cpm_simple",
       fems_to_run = c("mrmr_perc50")),
  
  #3
  #CFRDVsNGT
  list(phenotype_file_name = "data/formatted/phenotype.txt",
       read_count_dir_path = "data/formatted",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "CF_EV_zlogcpm",
       classification_criteria = "CFRDVsNGT",
       classes = c("NGT", "CFRD"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       norm = "norm_log_cpm_simple",
       fems_to_run = c("RF_RFE")),
  
  #4
  #CFRDVsNGT
  list(phenotype_file_name = "data/formatted/phenotype.txt",
       read_count_dir_path = "data/formatted",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "CF_EV_zlogcpm",
       classification_criteria = "CFRDVsNGT",
       classes = c("NGT", "CFRD"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       norm = "norm_log_cpm_simple",
       fems_to_run = c("ga_rf")),    
  
  

  #1
  #IGTVsNGT
  list(phenotype_file_name = "data/formatted/phenotype.txt",
       read_count_dir_path = "data/formatted",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "CF_EV_zlogcpm",
       classification_criteria = "IGTVsNGT",
       classes = c("NGT", "IGT"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       norm = "norm_log_cpm_simple",
       fems_to_run = c("all", 
                       "t-test", "wilcoxontest",
                       "ranger_impu_cor",
                       "ranger_pos_impu_cor",
                       "mrmr10", "mrmr75")),
  
  #2
  #IGTVsNGT
  list(phenotype_file_name = "data/formatted/phenotype.txt",
       read_count_dir_path = "data/formatted",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "CF_EV_zlogcpm",
       classification_criteria = "IGTVsNGT",
       classes = c("NGT", "IGT"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       norm = "norm_log_cpm_simple",
       fems_to_run = c("mrmr_perc50")),
  
  #3
  #IGTVsNGT
  list(phenotype_file_name = "data/formatted/phenotype.txt",
       read_count_dir_path = "data/formatted",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "CF_EV_zlogcpm",
       classification_criteria = "IGTVsNGT",
       classes = c("NGT", "IGT"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       norm = "norm_log_cpm_simple",
       fems_to_run = c("RF_RFE")),
  
  #4
  #IGTVsNGT
  list(phenotype_file_name = "data/formatted/phenotype.txt",
       read_count_dir_path = "data/formatted",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "CF_EV_zlogcpm",
       classification_criteria = "IGTVsNGT",
       classes = c("NGT", "IGT"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       norm = "norm_log_cpm_simple",
       fems_to_run = c("ga_rf")),
  
  
  
  
  #1
  #CFRDVsIGT
  list(phenotype_file_name = "data/formatted/phenotype.txt",
       read_count_dir_path = "data/formatted",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "CF_EV_zlogtmm",
       classification_criteria = "CFRDVsIGT",
       classes = c("IGT", "CFRD"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       norm = "norm_log_tmm",
       fems_to_run = c("all", 
                       "t-test", "wilcoxontest",
                       "ranger_impu_cor",
                       "ranger_pos_impu_cor",
                       "mrmr10", "mrmr75")),
  
  #2
  #CFRDVsIGT
  list(phenotype_file_name = "data/formatted/phenotype.txt",
       read_count_dir_path = "data/formatted",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "CF_EV_zlogtmm",
       classification_criteria = "CFRDVsIGT",
       classes = c("IGT", "CFRD"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       norm = "norm_log_tmm",
       fems_to_run = c("mrmr_perc50")),
  
  #3
  #CFRDVsIGT
  list(phenotype_file_name = "data/formatted/phenotype.txt",
       read_count_dir_path = "data/formatted",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "CF_EV_zlogtmm",
       classification_criteria = "CFRDVsIGT",
       classes = c("IGT", "CFRD"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       norm = "norm_log_tmm",
       fems_to_run = c("RF_RFE")),
  
  #4
  #CFRDVsIGT
  list(phenotype_file_name = "data/formatted/phenotype.txt",
       read_count_dir_path = "data/formatted",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "CF_EV_zlogtmm",
       classification_criteria = "CFRDVsIGT",
       classes = c("IGT", "CFRD"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       norm = "norm_log_tmm",
       fems_to_run = c("ga_rf")),  
  
  
  #1
  #CFRDVsNGT
  list(phenotype_file_name = "data/formatted/phenotype.txt",
       read_count_dir_path = "data/formatted",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "CF_EV_zlogtmm",
       classification_criteria = "CFRDVsNGT",
       classes = c("NGT", "CFRD"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       norm = "norm_log_tmm",
       fems_to_run = c("all", 
                       "t-test", "wilcoxontest",
                       "ranger_impu_cor",
                       "ranger_pos_impu_cor",
                       "mrmr10", "mrmr75")),
  
  #2
  #CFRDVsNGT
  list(phenotype_file_name = "data/formatted/phenotype.txt",
       read_count_dir_path = "data/formatted",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "CF_EV_zlogtmm",
       classification_criteria = "CFRDVsNGT",
       classes = c("NGT", "CFRD"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       norm = "norm_log_tmm",
       fems_to_run = c("mrmr_perc50")),
  
  #3
  #CFRDVsNGT
  list(phenotype_file_name = "data/formatted/phenotype.txt",
       read_count_dir_path = "data/formatted",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "CF_EV_zlogtmm",
       classification_criteria = "CFRDVsNGT",
       classes = c("NGT", "CFRD"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       norm = "norm_log_tmm",
       fems_to_run = c("RF_RFE")),
  
  #4
  #CFRDVsNGT
  list(phenotype_file_name = "data/formatted/phenotype.txt",
       read_count_dir_path = "data/formatted",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "CF_EV_zlogtmm",
       classification_criteria = "CFRDVsNGT",
       classes = c("NGT", "CFRD"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       norm = "norm_log_tmm",
       fems_to_run = c("ga_rf")),    
  
  
  
  #1
  #IGTVsNGT
  list(phenotype_file_name = "data/formatted/phenotype.txt",
       read_count_dir_path = "data/formatted",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "CF_EV_zlogtmm",
       classification_criteria = "IGTVsNGT",
       classes = c("NGT", "IGT"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       norm = "norm_log_tmm",
       fems_to_run = c("all", 
                       "t-test", "wilcoxontest",
                       "ranger_impu_cor",
                       "ranger_pos_impu_cor",
                       "mrmr10", "mrmr75")),
  
  #2
  #IGTVsNGT
  list(phenotype_file_name = "data/formatted/phenotype.txt",
       read_count_dir_path = "data/formatted",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "CF_EV_zlogtmm",
       classification_criteria = "IGTVsNGT",
       classes = c("NGT", "IGT"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       norm = "norm_log_tmm",
       fems_to_run = c("mrmr_perc50")),
  
  #3
  #IGTVsNGT
  list(phenotype_file_name = "data/formatted/phenotype.txt",
       read_count_dir_path = "data/formatted",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "CF_EV_zlogtmm",
       classification_criteria = "IGTVsNGT",
       classes = c("NGT", "IGT"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       norm = "norm_log_tmm",
       fems_to_run = c("RF_RFE")),
  
  #4
  #IGTVsNGT
  list(phenotype_file_name = "data/formatted/phenotype.txt",
       read_count_dir_path = "data/formatted",
       read_count_file_name = "umi_counts.csv",
       sep = ",",
       dataset_id = "CF_EV_zlogtmm",
       classification_criteria = "IGTVsNGT",
       classes = c("NGT", "IGT"),
       cores = 16,
       results_dir_path = "fem_pipeline_results",
       norm = "norm_log_tmm",
       fems_to_run = c("ga_rf")),
  
  
  
######## AU samples

#1
#CFRDVsIGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogcpm",
     classification_criteria = "CFRDVsIGT",
     classes = c("IGT", "CFRD"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_cpm_simple",
     fems_to_run = c("all", 
                     "t-test", "wilcoxontest",
                     "ranger_impu_cor",
                     "ranger_pos_impu_cor",
                     "mrmr10", "mrmr75")),

#2
#CFRDVsIGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogcpm",
     classification_criteria = "CFRDVsIGT",
     classes = c("IGT", "CFRD"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_cpm_simple",
     fems_to_run = c("mrmr_perc50")),

#3
#CFRDVsIGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogcpm",
     classification_criteria = "CFRDVsIGT",
     classes = c("IGT", "CFRD"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_cpm_simple",
     fems_to_run = c("RF_RFE")),

#4
#CFRDVsIGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogcpm",
     classification_criteria = "CFRDVsIGT",
     classes = c("IGT", "CFRD"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_cpm_simple",
     fems_to_run = c("ga_rf")),  


#1
#CFRDVsNGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogcpm",
     classification_criteria = "CFRDVsNGT",
     classes = c("NGT", "CFRD"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_cpm_simple",
     fems_to_run = c("all", 
                     "t-test", "wilcoxontest",
                     "ranger_impu_cor",
                     "ranger_pos_impu_cor",
                     "mrmr10", "mrmr75")),

#2
#CFRDVsNGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogcpm",
     classification_criteria = "CFRDVsNGT",
     classes = c("NGT", "CFRD"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_cpm_simple",
     fems_to_run = c("mrmr_perc50")),

#3
#CFRDVsNGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogcpm",
     classification_criteria = "CFRDVsNGT",
     classes = c("NGT", "CFRD"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_cpm_simple",
     fems_to_run = c("RF_RFE")),

#4
#CFRDVsNGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogcpm",
     classification_criteria = "CFRDVsNGT",
     classes = c("NGT", "CFRD"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_cpm_simple",
     fems_to_run = c("ga_rf")),    



#1
#IGTVsNGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogcpm",
     classification_criteria = "IGTVsNGT",
     classes = c("NGT", "IGT"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_cpm_simple",
     fems_to_run = c("all", 
                     "t-test", "wilcoxontest",
                     "ranger_impu_cor",
                     "ranger_pos_impu_cor",
                     "mrmr10", "mrmr75")),

#2
#IGTVsNGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogcpm",
     classification_criteria = "IGTVsNGT",
     classes = c("NGT", "IGT"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_cpm_simple",
     fems_to_run = c("mrmr_perc50")),

#3
#IGTVsNGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogcpm",
     classification_criteria = "IGTVsNGT",
     classes = c("NGT", "IGT"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_cpm_simple",
     fems_to_run = c("RF_RFE")),

#4
#IGTVsNGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogcpm",
     classification_criteria = "IGTVsNGT",
     classes = c("NGT", "IGT"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_cpm_simple",
     fems_to_run = c("ga_rf")),




#1
#CFRDVsIGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogtmm",
     classification_criteria = "CFRDVsIGT",
     classes = c("IGT", "CFRD"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_tmm",
     fems_to_run = c("all", 
                     "t-test", "wilcoxontest",
                     "ranger_impu_cor",
                     "ranger_pos_impu_cor",
                     "mrmr10", "mrmr75")),

#2
#CFRDVsIGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogtmm",
     classification_criteria = "CFRDVsIGT",
     classes = c("IGT", "CFRD"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_tmm",
     fems_to_run = c("mrmr_perc50")),

#3
#CFRDVsIGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogtmm",
     classification_criteria = "CFRDVsIGT",
     classes = c("IGT", "CFRD"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_tmm",
     fems_to_run = c("RF_RFE")),

#4
#CFRDVsIGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogtmm",
     classification_criteria = "CFRDVsIGT",
     classes = c("IGT", "CFRD"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_tmm",
     fems_to_run = c("ga_rf")),  


#1
#CFRDVsNGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogtmm",
     classification_criteria = "CFRDVsNGT",
     classes = c("NGT", "CFRD"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_tmm",
     fems_to_run = c("all", 
                     "t-test", "wilcoxontest",
                     "ranger_impu_cor",
                     "ranger_pos_impu_cor",
                     "mrmr10", "mrmr75")),

#2
#CFRDVsNGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogtmm",
     classification_criteria = "CFRDVsNGT",
     classes = c("NGT", "CFRD"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_tmm",
     fems_to_run = c("mrmr_perc50")),

#3
#CFRDVsNGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogtmm",
     classification_criteria = "CFRDVsNGT",
     classes = c("NGT", "CFRD"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_tmm",
     fems_to_run = c("RF_RFE")),

#4
#CFRDVsNGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogtmm",
     classification_criteria = "CFRDVsNGT",
     classes = c("NGT", "CFRD"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_tmm",
     fems_to_run = c("ga_rf")),    



#1
#IGTVsNGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogtmm",
     classification_criteria = "IGTVsNGT",
     classes = c("NGT", "IGT"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_tmm",
     fems_to_run = c("all", 
                     "t-test", "wilcoxontest",
                     "ranger_impu_cor",
                     "ranger_pos_impu_cor",
                     "mrmr10", "mrmr75")),

#2
#IGTVsNGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogtmm",
     classification_criteria = "IGTVsNGT",
     classes = c("NGT", "IGT"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_tmm",
     fems_to_run = c("mrmr_perc50")),

#3
#IGTVsNGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogtmm",
     classification_criteria = "IGTVsNGT",
     classes = c("NGT", "IGT"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_tmm",
     fems_to_run = c("RF_RFE")),

#4
#IGTVsNGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogtmm",
     classification_criteria = "IGTVsNGT",
     classes = c("NGT", "IGT"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_tmm",
     fems_to_run = c("ga_rf")),   



#subsets for AU cohort with tmm norm
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted/subset",
     read_count_file_name = "CF_EV_AU_zlogtmm_CFRDVsIGT_mrmr75_rangerpos_28.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogtmm_CFRDVsIGT_mrmr75_rangerpos_28",
     classification_criteria = "CFRDVsIGT",
     classes = c("IGT", "CFRD"),
     cores = 4,
     results_dir_path = "fem_pipeline_results_AU_subset",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_tmm",
     fems_to_run = c("all"),
     random_seed = 2000),

list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted/subset",
     read_count_file_name = "CF_EV_AU_zlogtmm_CFRDVsIGT_mrmr75_rangerpos_nopir_28.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogtmm_CFRDVsIGT_mrmr75_rangerpos_nopir_28",
     classification_criteria = "CFRDVsIGT",
     classes = c("IGT", "CFRD"),
     cores = 4,
     results_dir_path = "fem_pipeline_results_AU_subset",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_tmm",
     fems_to_run = c("all"),
     random_seed = 2000),

list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted/subset",
     read_count_file_name = "CF_EV_AU_zlogtmm_CFRDVsIGT_mrmr75_28.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogtmm_CFRDVsIGT_mrmr75_28",
     classification_criteria = "CFRDVsIGT",
     classes = c("IGT", "CFRD"),
     cores = 4,
     results_dir_path = "fem_pipeline_results_AU_subset",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_tmm",
     fems_to_run = c("all"),
     random_seed = 2000),

list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted/subset",
     read_count_file_name = "CF_EV_AU_zlogtmm_CFRDVsIGT_mrmr75_nopir_28.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogtmm_CFRDVsIGT_mrmr75_nopir_28",
     classification_criteria = "CFRDVsIGT",
     classes = c("IGT", "CFRD"),
     cores = 4,
     results_dir_path = "fem_pipeline_results_AU_subset",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_tmm",
     fems_to_run = c("all"),
     random_seed = 2000),

list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted/subset",
     read_count_file_name = "CF_EV_AU_zlogtmm_CFRDVsIGT_ranger_pos_28.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogtmm_CFRDVsIGT_ranger_pos_28",
     classification_criteria = "CFRDVsIGT",
     classes = c("IGT", "CFRD"),
     cores = 4,
     results_dir_path = "fem_pipeline_results_AU_subset",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_tmm",
     fems_to_run = c("all"),
     random_seed = 2000),

list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted/subset",
     read_count_file_name = "CF_EV_AU_zlogtmm_CFRDVsIGT_ranger_pos_nopir_28.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogtmm_CFRDVsIGT_ranger_pos_nopir_28",
     classification_criteria = "CFRDVsIGT",
     classes = c("IGT", "CFRD"),
     cores = 4,
     results_dir_path = "fem_pipeline_results_AU_subset",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_tmm",
     fems_to_run = c("all"),
     random_seed = 2000),
  


list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted/subset",
     read_count_file_name = "CF_EV_AU_zlogtmm_CFRDVsNGT_mrmr75_28.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogtmm_CFRDVsNGT_mrmr75_28",
     classification_criteria = "CFRDVsNGT",
     classes = c("NGT", "CFRD"),
     cores = 4,
     results_dir_path = "fem_pipeline_results_AU_subset",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_tmm",
     fems_to_run = c("all"),
     random_seed = 2000),

list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted/subset",
     read_count_file_name = "CF_EV_AU_zlogtmm_CFRDVsNGT_mrmr75_nopir_28.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogtmm_CFRDVsNGT_mrmr75_nopir_28",
     classification_criteria = "CFRDVsNGT",
     classes = c("NGT", "CFRD"),
     cores = 4,
     results_dir_path = "fem_pipeline_results_AU_subset",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_tmm",
     fems_to_run = c("all"),
     random_seed = 2000),



list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted/subset",
     read_count_file_name = "CF_EV_AU_zlogtmm_IGTVsNGT_mrmr75_28.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogtmm_IGTVsNGT_mrmr75_28",
     classification_criteria = "IGTVsNGT",
     classes = c("NGT", "IGT"),
     cores = 4,
     results_dir_path = "fem_pipeline_results_AU_subset",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_tmm",
     fems_to_run = c("all"),
     random_seed = 2000),

list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted/subset",
     read_count_file_name = "CF_EV_AU_zlogtmm_IGTVsNGT_mrmr75_nopir_28.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_zlogtmm_IGTVsNGT_mrmr75_nopir_28",
     classification_criteria = "IGTVsNGT",
     classes = c("NGT", "IGT"),
     cores = 4,
     results_dir_path = "fem_pipeline_results_AU_subset",
     filter_expression = expression(country == 'AU'),
     norm = "norm_log_tmm",
     fems_to_run = c("all"),
     random_seed = 2000),



#subsets for AU+DK cohort with tmm norm
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted/subset",
     read_count_file_name = "CF_EV_zlogtmm_CFRDVsIGT_mrmr75_28.csv",
     sep = ",",
     dataset_id = "CF_EV_zlogtmm_CFRDVsIGT_mrmr75_28",
     classification_criteria = "CFRDVsIGT",
     classes = c("IGT", "CFRD"),
     cores = 4,
     results_dir_path = "fem_pipeline_results_subset",
     norm = "norm_log_tmm",
     fems_to_run = c("all"),
     random_seed = 2000),

list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted/subset",
     read_count_file_name = "CF_EV_zlogtmm_CFRDVsIGT_mrmr75_nopir_28.csv",
     sep = ",",
     dataset_id = "CF_EV_zlogtmm_CFRDVsIGT_mrmr75_nopir_28",
     classification_criteria = "CFRDVsIGT",
     classes = c("IGT", "CFRD"),
     cores = 4,
     results_dir_path = "fem_pipeline_results_subset",
     norm = "norm_log_tmm",
     fems_to_run = c("all"),
     random_seed = 2000),



list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted/subset",
     read_count_file_name = "CF_EV_zlogtmm_CFRDVsNGT_mrmr75_28.csv",
     sep = ",",
     dataset_id = "CF_EV_zlogtmm_CFRDVsNGT_mrmr75_28",
     classification_criteria = "CFRDVsNGT",
     classes = c("NGT", "CFRD"),
     cores = 4,
     results_dir_path = "fem_pipeline_results_subset",
     norm = "norm_log_tmm",
     fems_to_run = c("all"),
     random_seed = 2000),

list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted/subset",
     read_count_file_name = "CF_EV_zlogtmm_CFRDVsNGT_mrmr75_nopir_28.csv",
     sep = ",",
     dataset_id = "CF_EV_zlogtmm_CFRDVsNGT_mrmr75_nopir_28",
     classification_criteria = "CFRDVsNGT",
     classes = c("NGT", "CFRD"),
     cores = 4,
     results_dir_path = "fem_pipeline_results_subset",
     norm = "norm_log_tmm",
     fems_to_run = c("all"),
     random_seed = 2000),



list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted/subset",
     read_count_file_name = "CF_EV_zlogtmm_IGTVsNGT_mrmr75_28.csv",
     sep = ",",
     dataset_id = "CF_EV_zlogtmm_IGTVsNGT_mrmr75_28",
     classification_criteria = "IGTVsNGT",
     classes = c("NGT", "IGT"),
     cores = 4,
     results_dir_path = "fem_pipeline_results_subset",
     norm = "norm_log_tmm",
     fems_to_run = c("all"),
     random_seed = 2000),

list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted/subset",
     read_count_file_name = "CF_EV_zlogtmm_IGTVsNGT_mrmr75_nopir_28.csv",
     sep = ",",
     dataset_id = "CF_EV_zlogtmm_IGTVsNGT_mrmr75_nopir_28",
     classification_criteria = "IGTVsNGT",
     classes = c("NGT", "IGT"),
     cores = 4,
     results_dir_path = "fem_pipeline_results_subset",
     norm = "norm_log_tmm",
     fems_to_run = c("all"),
     random_seed = 2000),



#AU samples with just log tmm i.e. no scaling after log tmm 

#1
#CFRDVsIGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_logtmm",
     classification_criteria = "CFRDVsIGT",
     classes = c("IGT", "CFRD"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU_logtmm",
     filter_expression = expression(country == 'AU'),
     norm = "log_tmm",
     fems_to_run = c("all", 
                     "t-test", "wilcoxontest",
                     "ranger_impu_cor",
                     "ranger_pos_impu_cor",
                     "mrmr10", "mrmr75")),

#2
#CFRDVsIGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_logtmm",
     classification_criteria = "CFRDVsIGT",
     classes = c("IGT", "CFRD"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU_logtmm",
     filter_expression = expression(country == 'AU'),
     norm = "log_tmm",
     fems_to_run = c("mrmr_perc50")),

#3
#CFRDVsIGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_logtmm",
     classification_criteria = "CFRDVsIGT",
     classes = c("IGT", "CFRD"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU_logtmm",
     filter_expression = expression(country == 'AU'),
     norm = "log_tmm",
     fems_to_run = c("RF_RFE")),

#4
#CFRDVsIGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_logtmm",
     classification_criteria = "CFRDVsIGT",
     classes = c("IGT", "CFRD"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU_logtmm",
     filter_expression = expression(country == 'AU'),
     norm = "log_tmm",
     fems_to_run = c("ga_rf")),  


#1
#CFRDVsNGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_logtmm",
     classification_criteria = "CFRDVsNGT",
     classes = c("NGT", "CFRD"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU_logtmm",
     filter_expression = expression(country == 'AU'),
     norm = "log_tmm",
     fems_to_run = c("all", 
                     "t-test", "wilcoxontest",
                     "ranger_impu_cor",
                     "ranger_pos_impu_cor",
                     "mrmr10", "mrmr75")),

#2
#CFRDVsNGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_logtmm",
     classification_criteria = "CFRDVsNGT",
     classes = c("NGT", "CFRD"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU_logtmm",
     filter_expression = expression(country == 'AU'),
     norm = "log_tmm",
     fems_to_run = c("mrmr_perc50")),

#3
#CFRDVsNGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_logtmm",
     classification_criteria = "CFRDVsNGT",
     classes = c("NGT", "CFRD"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU_logtmm",
     filter_expression = expression(country == 'AU'),
     norm = "log_tmm",
     fems_to_run = c("RF_RFE")),

#4
#CFRDVsNGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_logtmm",
     classification_criteria = "CFRDVsNGT",
     classes = c("NGT", "CFRD"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU_logtmm",
     filter_expression = expression(country == 'AU'),
     norm = "log_tmm",
     fems_to_run = c("ga_rf")),    



#1
#IGTVsNGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_logtmm",
     classification_criteria = "IGTVsNGT",
     classes = c("NGT", "IGT"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU_logtmm",
     filter_expression = expression(country == 'AU'),
     norm = "log_tmm",
     fems_to_run = c("all", 
                     "t-test", "wilcoxontest",
                     "ranger_impu_cor",
                     "ranger_pos_impu_cor",
                     "mrmr10", "mrmr75")),

#2
#IGTVsNGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_logtmm",
     classification_criteria = "IGTVsNGT",
     classes = c("NGT", "IGT"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU_logtmm",
     filter_expression = expression(country == 'AU'),
     norm = "log_tmm",
     fems_to_run = c("mrmr_perc50")),

#3
#IGTVsNGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_logtmm",
     classification_criteria = "IGTVsNGT",
     classes = c("NGT", "IGT"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU_logtmm",
     filter_expression = expression(country == 'AU'),
     norm = "log_tmm",
     fems_to_run = c("RF_RFE")),

#4
#IGTVsNGT
list(phenotype_file_name = "data/formatted/phenotype.txt",
     read_count_dir_path = "data/formatted",
     read_count_file_name = "umi_counts.csv",
     sep = ",",
     dataset_id = "CF_EV_AU_logtmm",
     classification_criteria = "IGTVsNGT",
     classes = c("NGT", "IGT"),
     cores = 16,
     results_dir_path = "fem_pipeline_results_AU_logtmm",
     filter_expression = expression(country == 'AU'),
     norm = "log_tmm",
     fems_to_run = c("ga_rf"))


)  