library(sva)

base_dir <- "~/UNSW/VafaeeLab/CysticFibrosisGroup/ExoCF/CFRD_EV_biomarker/"
setwd(base_dir)

source("utils_diff.R")
source("utils.R")


data <- read.table(file = "data/proteomics/all/proteinGroups.txt", sep = "\t", header = TRUE)

data <- data[data$Reverse != "+" & data$Potential.contaminant != "+" & data$Unique.peptides >= 2,]
data <- data[data$Gene.names != "",] #removes rows with no Gene names
data <- data[complete.cases(data),] #removes NAs at the end 
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";") 

meta_data <- read.csv("data/proteomics/prot_metadata_all_2023Oct.csv") %>%
  mutate("PreModulatorVsPostModulator" = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ "PostModulator",
                                                   TRUE ~ "PreModulator")) %>%
  mutate(modstatus_condition_country = paste(PreModulatorVsPostModulator, condition, country, sep = "_")) %>%
  mutate(modstatus_condition = paste(PreModulatorVsPostModulator, condition, sep = "_")) %>%
  dplyr::rename(c("disease_status" = "condition"))

summary(factor(meta_data$modstatus_condition_country))

# PostModulator_CFRD_AU               PostModulator_CFRD_DK                PostModulator_IGT_AU 
# 11                                  35                                   3 
# PostModulator_IGT_DK                PostModulator_NGT_AU                PostModulator_NGT_DK 
# 30                                  13                                  62 
# PostModulator_UNKNOWN TO PREDICT_AU                PreModulator_CFRD_AU                PreModulator_CFRD_DK 
# 12                                  24                                  14 
# PreModulator_HC_AU                 PreModulator_IGT_AU                 PreModulator_IGT_DK 
# 13                                  12                                  26 
# PreModulator_NGT_AU                 PreModulator_NGT_DK  PreModulator_UNKNOWN TO PREDICT_AU 
# 30                                  26                                  22 

meta_data <- meta_data %>%
  filter(disease_status != "UNKNOWN TO PREDICT" & disease_status != "HC")
summary(factor(meta_data$modstatus_condition_country))

# PostModulator_CFRD_AU PostModulator_CFRD_DK  PostModulator_IGT_AU  PostModulator_IGT_DK  PostModulator_NGT_AU  PostModulator_NGT_DK 
# 11                    35                     3                    30                    13                    62 
# PreModulator_CFRD_AU  PreModulator_CFRD_DK   PreModulator_IGT_AU   PreModulator_IGT_DK   PreModulator_NGT_AU   PreModulator_NGT_DK 
# 24                    14                    12                    26                    30                    26 


#filter proteins using modstatus_condition_country and then use that for modstatus_condition


meta_data <- add_replicate_column(meta_data, "modstatus_condition_country")
dim(meta_data)
# [1] 286  37

exp_design <- meta_data %>%
  dplyr::select(label, condition, replicate) %>%
  mutate(label = as.character(label))
exp_design$label <- paste("LFQ.intensity.",exp_design$label,sep="")

columns <- grep("LFQ.intensity.", colnames(data_unique))

se <- make_se(data_unique, columns, exp_design)

# plot_coverage(se)

dim(assays(se)[[1]])
# [1] 1658  286

summary(summary(factor(meta_data$condition)) / dim(meta_data)[1])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01049 0.04458 0.08741 0.08333 0.10490 0.21678


filt <- filter_proteins(se, "fraction", min = 0.01049/2)
dim(assays(filt)[[1]])
# [1] 1536  286

filt_data <- get_df_long(filt) %>%
  dplyr::select(c(label, name, intensity)) %>%
  pivot_wider(names_from = name, values_from = intensity) %>%
  mutate(label = sub("LFQ.intensity.", "lfq_", label)) %>%
  column_to_rownames("label") 


sum(!is.na(filt_data)) / (nrow(filt_data) * ncol(filt_data))
#0.2354267

#to be noted : percentage of non-NAs are really low

filt_data_from_se <- as.data.frame(assays(filt)[[1]])
dim(filt_data_from_se)
# [1] 1536  286

########################################
#now use modstatus_condition and get filtered proteins
meta_data <- meta_data %>%
  dplyr::rename(c("modstatus_condition_country" = "condition"))
meta_data <- add_replicate_column(meta_data, "modstatus_condition")

dim(meta_data)
# [1] 286  37

exp_design <- meta_data %>%
  dplyr::select(label, condition, replicate) %>%
  mutate(label = as.character(label))
exp_design$label <- paste("LFQ.intensity.",exp_design$label,sep="")

columns <- grep("LFQ.intensity.", colnames(data_unique))

se <- make_se(data_unique %>% dplyr::filter(name %in% rownames(filt_data_from_se)), columns, exp_design)

# plot_coverage(se)

dim(assays(se)[[1]])
# [1] 1536  286
########################################

norm <- normalize_vsn(se)
meanSdPlot(norm)

#Note - don't run thsi again. Use the saved rds object
# imputed <- impute(norm, fun = "bpca")

plot_imputation(norm, imputed)
ggsave("imputation_plot_regular_DEP.png")

data.se <- get_df_long(se) %>%
  select(label, condition, replicate, name, intensity) %>%
  mutate(label = sub("LFQ.intensity.", "", label, fixed = TRUE)) %>%
  tidyr::separate(label, into = c(NA, NA, NA, "label"), sep = "-")
data.imp <- get_df_long(imputed) %>%
  select(label, condition, replicate, name, intensity) %>%
  mutate(label = sub("LFQ.intensity.", "", label, fixed = TRUE)) %>%
  tidyr::separate(label, into = c(NA, NA, NA, "label"), sep = "-")

sum(is.na(data.se))
# 335874
sum(is.na(data.imp))
# 0

######################################################################

#DE analysis
diffex(use_adj_pval = TRUE, imputed, cond1 = "PreModulator_CFRD", cond2 = "PreModulator_IGT", 
       file_path = "data/de_results_2024/proteomics/DEP_regular/padj/", protein_names_file_path = "data/proteomics/protein_names.csv")
diffex(imputed, cond1 = "PreModulator_CFRD", cond2 = "PreModulator_IGT",
       file_path = "data/de_results_2024/proteomics/DEP_regular/p/", 
       protein_names_file_path = "data/proteomics/protein_names.csv", 
       fc = 1.5)

diffex(use_adj_pval = TRUE, imputed, cond1 = "PreModulator_CFRD", cond2 = "PreModulator_NGT", 
       file_path = "data/de_results_2024/proteomics/DEP_regular/padj/", protein_names_file_path = "data/proteomics/protein_names.csv")
diffex(imputed, cond1 = "PreModulator_CFRD", cond2 = "PreModulator_NGT",
       file_path = "data/de_results_2024/proteomics/DEP_regular/p/", 
       protein_names_file_path = "data/proteomics/protein_names.csv", 
       fc = 1.5)

diffex(use_adj_pval = TRUE, imputed, cond1 = "PreModulator_IGT", cond2 = "PreModulator_NGT", 
       file_path = "data/de_results_2024/proteomics/DEP_regular/padj/", protein_names_file_path = "data/proteomics/protein_names.csv")
diffex(imputed, cond1 = "PreModulator_IGT", cond2 = "PreModulator_NGT",
       file_path = "data/de_results_2024/proteomics/DEP_regular/p/", 
       protein_names_file_path = "data/proteomics/protein_names.csv", 
       fc = 1.5)


diffex(use_adj_pval = TRUE, imputed, cond1 = "PostModulator_CFRD", cond2 = "PreModulator_CFRD", 
       file_path = "data/de_results_2024/proteomics/DEP_regular/padj/", protein_names_file_path = "data/proteomics/protein_names.csv")
diffex(imputed, cond1 = "PostModulator_CFRD", cond2 = "PreModulator_CFRD",
       file_path = "data/de_results_2024/proteomics/DEP_regular/p/", 
       protein_names_file_path = "data/proteomics/protein_names.csv", 
       fc = 1.5)

diffex(use_adj_pval = TRUE, imputed, cond1 = "PostModulator_IGT", cond2 = "PreModulator_IGT", 
       file_path = "data/de_results_2024/proteomics/DEP_regular/padj/", protein_names_file_path = "data/proteomics/protein_names.csv")
diffex(imputed, cond1 = "PostModulator_IGT", cond2 = "PreModulator_IGT",
       file_path = "data/de_results_2024/proteomics/DEP_regular/p/", 
       protein_names_file_path = "data/proteomics/protein_names.csv", 
       fc = 1.5)

diffex(use_adj_pval = TRUE, imputed, cond1 = "PostModulator_NGT", cond2 = "PreModulator_NGT", 
       file_path = "data/de_results_2024/proteomics/DEP_regular/padj/", protein_names_file_path = "data/proteomics/protein_names.csv")
diffex(imputed, cond1 = "PostModulator_NGT", cond2 = "PreModulator_NGT",
       file_path = "data/de_results_2024/proteomics/DEP_regular/p/", 
       protein_names_file_path = "data/proteomics/protein_names.csv", 
       fc = 1.5)

saveRDS(imputed, "data/de_results_2024/proteomics/DEP_regular/imputed.rds")

######################################

#DE analysis using preprocessing done in the prediction pipeline

data <- read.csv("data/proteomics/data_333samples_imputed_mf.csv", row.names = 1)

# filter out HC and unknown
# quantile norm
# combat

phenotype <- read.table("data/formatted/prot_phenotype_333_2024Jan.txt", header = TRUE, sep = "\t") %>%
  mutate(modstatus_condition_country = paste(PreModulatorVsPostModulator, condition, country, sep = "_")) %>%
  mutate(modstatus_condition = paste(PreModulatorVsPostModulator, condition, sep = "_"))


# %>%
#   dplyr::rename(c("disease_status" = "condition"))

phenotype <- phenotype %>%
  filter(condition != "UNKNOWN TO PREDICT" & condition != "HC")
summary(factor(phenotype$modstatus_condition_country))
summary(factor(phenotype$modstatus_condition))

# phenotype <- add_replicate_column(phenotype, "modstatus_condition")
# summary(factor(phenotype$condition))  

data <- data[, phenotype$Sample]

###########
#normalize

#adapted from https://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/
data.rank <- apply(data, 2, rank, ties.method="average")
data.sorted <- data.frame(apply(data, 2, sort))
data.mean <- apply(data.sorted, 1, mean)
index_to_mean <- function(index, data_mean){
  #index can be int or int+0.5
  #if int+0.5, take average of the numbers in those positions
  int.result <- data_mean[index]
  index.int <- floor(index)
  #some of the values in point5.result might be NA
  #but they won't be chosen
  point5.result <- (data_mean[index.int] + data_mean[index.int+1])/2
  point5.indices <- index%%1 != 0
  result <- int.result
  result[point5.indices] <- point5.result[point5.indices]
  return (result)
}
data.norm <- apply(data.rank, 2, index_to_mean, data_mean = data.mean)
rownames(data.norm) <- rownames(data)
data <- data.norm    

all.equal(colnames(data), phenotype$Sample)

##############################################################################

create_dim_red_plots(comparison = NA,
                     classes = c("CFRD", "IGT"), 
                     class_colours = c("red", "orange"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics",
                     perform_filter = FALSE,
                     colour_column = "PreModulatorVsPostModulator", point_border_colours = c("black", "green"),
                     data = data,
                     phenotype = phenotype,
                     filter_post_modulator = FALSE,
                     custom_title = "1 Proteomics CFRD, IGT pre and post modulator samples with Quantile norm", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("CFRD", "NGT"), 
                     class_colours = c("red", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics",
                     perform_filter = FALSE,
                     colour_column = "PreModulatorVsPostModulator", point_border_colours = c("black", "green"),
                     data = data,
                     phenotype = phenotype,
                     filter_post_modulator = FALSE,
                     custom_title = "2 Proteomics CFRD, NGT pre and post modulator samples with Quantile norm", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("IGT", "NGT"), 
                     class_colours = c("orange", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics",
                     perform_filter = FALSE,
                     colour_column = "PreModulatorVsPostModulator", point_border_colours = c("black", "green"),
                     data = data,
                     phenotype = phenotype,
                     filter_post_modulator = FALSE,
                     custom_title = "3 Proteomics IGT, NGT pre and post modulator samples with Quantile norm", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PostModulator_CFRD", "PreModulator_NGT"), 
                     class_colours = c("red", "indianred", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics",
                     dimred_plot_width_cm = 36,
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "4 Proteomics shift from CFRD to NGT with Quantile norm", 
                     combat = FALSE)

create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PreModulator_IGT"), 
                     class_colours = c("red", "orange"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "5 Proteomics premod CFRD IGT with Quantile norm", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PreModulator_NGT"), 
                     class_colours = c("red", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "6 Proteomics premod CFRD NGT with Quantile norm", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_IGT", "PreModulator_NGT"), 
                     class_colours = c("orange", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "7 Proteomics premod IGT NGT with Quantile norm", 
                     combat = FALSE)

create_dim_red_plots(comparison = NA,
                     classes = c("PostModulator_CFRD", "PostModulator_IGT"), 
                     class_colours = c("red", "orange"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "8 Proteomics postmod CFRD IGT with Quantile norm", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PostModulator_CFRD", "PostModulator_NGT"), 
                     class_colours = c("red", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "9 Proteomics postmod CFRD NGT with Quantile norm", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PostModulator_IGT", "PostModulator_NGT"), 
                     class_colours = c("orange", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "10 Proteomics postmod IGT NGT with Quantile norm", 
                     combat = FALSE)

##############################################################################

#perform combat
data_of_interest <- as.data.frame(as.matrix(data))

all.equal(colnames(data_of_interest), phenotype$Sample)

data_of_interest.combat = ComBat(dat = data_of_interest, 
                                 batch = phenotype$country)
data_of_interest.combat <- as.data.frame(as.matrix(data_of_interest.combat))
data <- data_of_interest.combat




create_dim_red_plots(comparison = NA,
                     classes = c("CFRD", "IGT"), 
                     class_colours = c("red", "orange"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected",
                     perform_filter = FALSE,
                     colour_column = "PreModulatorVsPostModulator", point_border_colours = c("black", "green"),
                     data = data,
                     phenotype = phenotype,
                     filter_post_modulator = FALSE,
                     custom_title = "1 Proteomics CFRD, IGT pre and post modulator samples with Quantile norm and ComBat", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("CFRD", "NGT"), 
                     class_colours = c("red", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected",
                     perform_filter = FALSE,
                     colour_column = "PreModulatorVsPostModulator", point_border_colours = c("black", "green"),
                     data = data,
                     phenotype = phenotype,
                     filter_post_modulator = FALSE,
                     custom_title = "2 Proteomics CFRD, NGT pre and post modulator samples with Quantile norm and ComBat", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("IGT", "NGT"), 
                     class_colours = c("orange", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected",
                     perform_filter = FALSE,
                     colour_column = "PreModulatorVsPostModulator", point_border_colours = c("black", "green"),
                     data = data,
                     phenotype = phenotype,
                     filter_post_modulator = FALSE,
                     custom_title = "3 Proteomics IGT, NGT pre and post modulator samples with Quantile norm and ComBat", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PostModulator_CFRD", "PreModulator_NGT"), 
                     class_colours = c("red", "indianred", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected",
                     dimred_plot_width_cm = 36,
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "4 Proteomics shift from CFRD to NGT with Quantile norm and ComBat", 
                     combat = FALSE)

create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PreModulator_IGT"), 
                     class_colours = c("red", "orange"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "5 Proteomics premod CFRD IGT with Quantile norm and ComBat", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PreModulator_NGT"), 
                     class_colours = c("red", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "6 Proteomics premod CFRD NGT with Quantile norm and ComBat", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_IGT", "PreModulator_NGT"), 
                     class_colours = c("orange", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "7 Proteomics premod IGT NGT with Quantile norm and ComBat", 
                     combat = FALSE)

create_dim_red_plots(comparison = NA,
                     classes = c("PostModulator_CFRD", "PostModulator_IGT"), 
                     class_colours = c("red", "orange"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "8 Proteomics postmod CFRD IGT with Quantile norm and ComBat", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PostModulator_CFRD", "PostModulator_NGT"), 
                     class_colours = c("red", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "9 Proteomics postmod CFRD NGT with Quantile norm and ComBat", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PostModulator_IGT", "PostModulator_NGT"), 
                     class_colours = c("orange", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "10 Proteomics postmod IGT NGT with Quantile norm and ComBat", 
                     combat = FALSE)


#perform combat with batch_name
data_of_interest <- as.data.frame(as.matrix(data))

all.equal(colnames(data_of_interest), phenotype$Sample)

data_of_interest.combat = ComBat(dat = data_of_interest, 
                                 batch = phenotype$batch_name)
data_of_interest.combat <- as.data.frame(as.matrix(data_of_interest.combat))
data <- data_of_interest.combat

create_dim_red_plots(comparison = NA,
                     classes = c("CFRD", "IGT"), 
                     class_colours = c("red", "orange"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected_twice",
                     perform_filter = FALSE,
                     colour_column = "PreModulatorVsPostModulator", point_border_colours = c("black", "green"),
                     data = data,
                     phenotype = phenotype,
                     filter_post_modulator = FALSE,
                     custom_title = "1 Proteomics CFRD, IGT pre and post modulator samples with Quantile norm and ComBat twice", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("CFRD", "NGT"), 
                     class_colours = c("red", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected_twice",
                     perform_filter = FALSE,
                     colour_column = "PreModulatorVsPostModulator", point_border_colours = c("black", "green"),
                     data = data,
                     phenotype = phenotype,
                     filter_post_modulator = FALSE,
                     custom_title = "2 Proteomics CFRD, NGT pre and post modulator samples with Quantile norm and ComBat twice", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("IGT", "NGT"), 
                     class_colours = c("orange", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected_twice",
                     perform_filter = FALSE,
                     colour_column = "PreModulatorVsPostModulator", point_border_colours = c("black", "green"),
                     data = data,
                     phenotype = phenotype,
                     filter_post_modulator = FALSE,
                     custom_title = "3 Proteomics IGT, NGT pre and post modulator samples with Quantile norm and ComBat twice", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PostModulator_CFRD", "PreModulator_NGT"), 
                     class_colours = c("red", "indianred", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected_twice",
                     dimred_plot_width_cm = 36,
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "4 Proteomics shift from CFRD to NGT with Quantile norm and ComBat twice", 
                     combat = FALSE)

create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PreModulator_IGT"), 
                     class_colours = c("red", "orange"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected_twice",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "5 Proteomics premod CFRD IGT with Quantile norm and ComBat twice", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PreModulator_NGT"), 
                     class_colours = c("red", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected_twice",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "6 Proteomics premod CFRD NGT with Quantile norm and ComBat twice", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_IGT", "PreModulator_NGT"), 
                     class_colours = c("orange", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected_twice",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "7 Proteomics premod IGT NGT with Quantile norm and ComBat twice", 
                     combat = FALSE)

create_dim_red_plots(comparison = NA,
                     classes = c("PostModulator_CFRD", "PostModulator_IGT"), 
                     class_colours = c("red", "orange"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected_twice",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "8 Proteomics postmod CFRD IGT with Quantile norm and ComBat twice", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PostModulator_CFRD", "PostModulator_NGT"), 
                     class_colours = c("red", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected_twice",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "9 Proteomics postmod CFRD NGT with Quantile norm and ComBat twice", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PostModulator_IGT", "PostModulator_NGT"), 
                     class_colours = c("orange", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected_twice",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "10 Proteomics postmod IGT NGT with Quantile norm and ComBat twice", 
                     combat = FALSE)


################################
# combat with country+batch_name

data <- read.csv("data/proteomics/data_333samples_imputed_mf.csv", row.names = 1)

# filter out HC and unknown
# quantile norm
# combat

phenotype <- read.table("data/formatted/prot_phenotype_333_2024Jan.txt", header = TRUE, sep = "\t") %>%
  mutate(modstatus_condition_country = paste(PreModulatorVsPostModulator, condition, country, sep = "_")) %>%
  mutate(modstatus_condition = paste(PreModulatorVsPostModulator, condition, sep = "_"))


# %>%
#   dplyr::rename(c("disease_status" = "condition"))

phenotype <- phenotype %>%
  filter(condition != "UNKNOWN TO PREDICT" & condition != "HC")
summary(factor(phenotype$modstatus_condition_country))
summary(factor(phenotype$modstatus_condition))

# phenotype <- add_replicate_column(phenotype, "modstatus_condition")
# summary(factor(phenotype$condition))  

data <- data[, phenotype$Sample]

###########
#normalize

#adapted from https://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/
data.rank <- apply(data, 2, rank, ties.method="average")
data.sorted <- data.frame(apply(data, 2, sort))
data.mean <- apply(data.sorted, 1, mean)
index_to_mean <- function(index, data_mean){
  #index can be int or int+0.5
  #if int+0.5, take average of the numbers in those positions
  int.result <- data_mean[index]
  index.int <- floor(index)
  #some of the values in point5.result might be NA
  #but they won't be chosen
  point5.result <- (data_mean[index.int] + data_mean[index.int+1])/2
  point5.indices <- index%%1 != 0
  result <- int.result
  result[point5.indices] <- point5.result[point5.indices]
  return (result)
}
data.norm <- apply(data.rank, 2, index_to_mean, data_mean = data.mean)
rownames(data.norm) <- rownames(data)
data <- data.norm    

all.equal(colnames(data), phenotype$Sample)

phenotype <- phenotype %>%
  mutate(country_batch_name = paste(country, batch_name, sep = "_"))
summary(factor(phenotype$country_batch_name))
# AU_main  DK_main   DK_new DK_other 
# 93      112       63       18

#perform combat with country+batch_name
data_of_interest <- as.data.frame(as.matrix(data))

all.equal(colnames(data_of_interest), phenotype$Sample)

data_of_interest.combat = ComBat(dat = data_of_interest, 
                                 batch = phenotype$country_batch_name)
data_of_interest.combat <- as.data.frame(as.matrix(data_of_interest.combat))
data <- data_of_interest.combat


create_dim_red_plots(comparison = NA,
                     classes = c("CFRD", "IGT"), 
                     class_colours = c("red", "orange"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected_two_in_one_shot",
                     perform_filter = FALSE,
                     colour_column = "PreModulatorVsPostModulator", point_border_colours = c("black", "green"),
                     data = data,
                     phenotype = phenotype,
                     filter_post_modulator = FALSE,
                     custom_title = "1 Proteomics CFRD, IGT pre and post modulator samples with Quantile norm and ComBat two_in_one_shot", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("CFRD", "NGT"), 
                     class_colours = c("red", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected_two_in_one_shot",
                     perform_filter = FALSE,
                     colour_column = "PreModulatorVsPostModulator", point_border_colours = c("black", "green"),
                     data = data,
                     phenotype = phenotype,
                     filter_post_modulator = FALSE,
                     custom_title = "2 Proteomics CFRD, NGT pre and post modulator samples with Quantile norm and ComBat two_in_one_shot", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("IGT", "NGT"), 
                     class_colours = c("orange", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected_two_in_one_shot",
                     perform_filter = FALSE,
                     colour_column = "PreModulatorVsPostModulator", point_border_colours = c("black", "green"),
                     data = data,
                     phenotype = phenotype,
                     filter_post_modulator = FALSE,
                     custom_title = "3 Proteomics IGT, NGT pre and post modulator samples with Quantile norm and ComBat two_in_one_shot", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PostModulator_CFRD", "PreModulator_NGT"), 
                     class_colours = c("red", "indianred", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected_two_in_one_shot",
                     dimred_plot_width_cm = 36,
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "4 Proteomics shift from CFRD to NGT with Quantile norm and ComBat two_in_one_shot", 
                     combat = FALSE)


create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PostModulator_CFRD", "PreModulator_NGT"), 
                     class_colours = c("red", "blue", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/shift",
                     dimred_plot_width_cm = 36,
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "Proteomics shift from CFRD to NGT with Quantile norm and ComBat two_in_one_shot with mean", 
                     combat = FALSE, simplified_with_mean = TRUE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PostModulator_CFRD", "PreModulator_NGT"), 
                     class_colours = c("red", "blue", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/shift",
                     dimred_plot_width_cm = 36,
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "Proteomics shift from CFRD to NGT with Quantile norm and ComBat two_in_one_shot", 
                     combat = FALSE, simplified = TRUE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PostModulator_CFRD", "PreModulator_NGT"), 
                     class_colours = c("red", "blue", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/shift",
                     dimred_plot_width_cm = 36,
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "Proteomics shift from CFRD to NGT with Quantile norm and ComBat two_in_one_shot with problem sample names", 
                     combat = FALSE, simplified = TRUE, shownames = TRUE)

create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PostModulator_CFRD", 
                                 "PreModulator_NGT", "PostModulator_NGT",
                                 "PreModulator_IGT", "PostModulator_IGT"), 
                     class_colours = c("red", "blue", 
                                       "yellow", "green",
                                       "orange", "purple"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/shift",
                     dimred_plot_width_cm = 36,
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "Proteomics pre and post modulator with Quantile norm and ComBat two_in_one_shot", 
                     combat = FALSE, simplified = TRUE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PostModulator_CFRD", 
                                 "PreModulator_NGT", "PostModulator_NGT",
                                 "PreModulator_IGT", "PostModulator_IGT"), 
                     class_colours = c("red", "blue", 
                                       "yellow", "green",
                                       "orange", "purple"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/shift",
                     dimred_plot_width_cm = 36,
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "Proteomics pre and post modulator with Quantile norm and ComBat two_in_one_shot with median", 
                     combat = FALSE, simplified_with_mean = TRUE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PostModulator_CFRD", 
                                 "PreModulator_NGT", "PostModulator_NGT",
                                 "PreModulator_IGT", "PostModulator_IGT"), 
                     class_colours = c("red", "blue", 
                                       "yellow", "green",
                                       "orange", "purple"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/shift",
                     dimred_plot_width_cm = 36,
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "Proteomics pre and post modulator with Quantile norm and ComBat two_in_one_shot with mean", 
                     combat = FALSE, simplified_with_mean = TRUE)


create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PreModulator_IGT"), 
                     class_colours = c("red", "orange"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected_two_in_one_shot",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "5 Proteomics premod CFRD IGT with Quantile norm and ComBat two_in_one_shot", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PreModulator_NGT"), 
                     class_colours = c("red", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected_two_in_one_shot",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "6 Proteomics premod CFRD NGT with Quantile norm and ComBat two_in_one_shot", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_IGT", "PreModulator_NGT"), 
                     class_colours = c("orange", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected_two_in_one_shot",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "7 Proteomics premod IGT NGT with Quantile norm and ComBat two_in_one_shot", 
                     combat = FALSE)

create_dim_red_plots(comparison = NA,
                     classes = c("PostModulator_CFRD", "PostModulator_IGT"), 
                     class_colours = c("red", "orange"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected_two_in_one_shot",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "8 Proteomics postmod CFRD IGT with Quantile norm and ComBat two_in_one_shot", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PostModulator_CFRD", "PostModulator_NGT"), 
                     class_colours = c("red", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected_two_in_one_shot",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "9 Proteomics postmod CFRD NGT with Quantile norm and ComBat two_in_one_shot", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PostModulator_IGT", "PostModulator_NGT"), 
                     class_colours = c("orange", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/proteomics/batch_corrected_two_in_one_shot",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "10 Proteomics postmod IGT NGT with Quantile norm and ComBat two_in_one_shot", 
                     combat = FALSE)


#use this combat 2 in 1 shot for DE analysis

# using DEP with this filtered data and manually creating SummarizedExperiment object 
#           without using make_se (so as to avoid log) was turning out to be difficult.
#           So directly using limma without voom as used within DEP

all.equal(phenotype$Sample, colnames(data))
# TRUE

summary(factor(phenotype$condition))
# CFRD  IGT  NGT 
# 84   71  131 

summary(factor(phenotype$modstatus_condition))
# PostModulator_CFRD  PostModulator_IGT  PostModulator_NGT  PreModulator_CFRD   PreModulator_IGT   PreModulator_NGT 
# 46                 33                 75                 38                 38                 56 

model_matrix <- model.matrix(~ 0 + phenotype$modstatus_condition)
colnames(model_matrix) <- sub("phenotype$modstatus_condition", "", colnames(model_matrix), fixed = TRUE)

contr_matrix <- makeContrasts(contrasts = "PreModulator_CFRD - PreModulator_IGT",
                              levels = colnames(model_matrix))
fit <- lmFit(data, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("protein")
result <- top.table %>%
  dplyr::select(protein, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = protein, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)

plot_volcano_and_save_DE(result, plot_title = "PreModulator_CFRD Vs PreModulator_IGT",
                         output_dir_path = "de_results_2024/proteomics/premod/p/",
                         plot_file_name = "PreModulator_CFRDVsPreModulator_IGT.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         molecule_names_file_path = "data/proteomics/protein_names.csv",
                         molecule_names_file_columns = c(4, 2),
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PreModulator_CFRD Vs PreModulator_IGT",
                         output_dir_path = "de_results_2024/proteomics/premod/padj/",
                         plot_file_name = "PreModulator_CFRDVsPreModulator_IGT.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         molecule_names_file_path = "data/proteomics/protein_names.csv",
                         molecule_names_file_columns = c(4, 2),
                         plot_width_cm = 25)


contr_matrix <- makeContrasts(contrasts = "PreModulator_CFRD - PreModulator_NGT",
                              levels = colnames(model_matrix))
fit <- lmFit(data, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("protein")
result <- top.table %>%
  dplyr::select(protein, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = protein, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)

plot_volcano_and_save_DE(result, plot_title = "PreModulator_CFRD Vs PreModulator_NGT",
                         output_dir_path = "de_results_2024/proteomics/premod/p/",
                         plot_file_name = "PreModulator_CFRDVsPreModulator_NGT.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         molecule_names_file_path = "data/proteomics/protein_names.csv",
                         molecule_names_file_columns = c(4, 2),
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PreModulator_CFRD Vs PreModulator_NGT",
                         output_dir_path = "de_results_2024/proteomics/premod/padj/",
                         plot_file_name = "PreModulator_CFRDVsPreModulator_NGT.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         molecule_names_file_path = "data/proteomics/protein_names.csv",
                         molecule_names_file_columns = c(4, 2),
                         plot_width_cm = 25)


contr_matrix <- makeContrasts(contrasts = "PreModulator_IGT - PreModulator_NGT",
                              levels = colnames(model_matrix))
fit <- lmFit(data, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("protein")
result <- top.table %>%
  dplyr::select(protein, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = protein, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)

plot_volcano_and_save_DE(result, plot_title = "PreModulator_IGT Vs PreModulator_NGT",
                         output_dir_path = "de_results_2024/proteomics/premod/p/",
                         plot_file_name = "PreModulator_IGTVsPreModulator_NGT.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         molecule_names_file_path = "data/proteomics/protein_names.csv",
                         molecule_names_file_columns = c(4, 2),
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PreModulator_IGT Vs PreModulator_NGT",
                         output_dir_path = "de_results_2024/proteomics/premod/padj/",
                         plot_file_name = "PreModulator_IGTVsPreModulator_NGT.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         molecule_names_file_path = "data/proteomics/protein_names.csv",
                         molecule_names_file_columns = c(4, 2),
                         plot_width_cm = 25)

################################################################################################

contr_matrix <- makeContrasts(contrasts = "PostModulator_CFRD - PostModulator_IGT",
                              levels = colnames(model_matrix))
fit <- lmFit(data, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("protein")
result <- top.table %>%
  dplyr::select(protein, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = protein, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)

plot_volcano_and_save_DE(result, plot_title = "PostModulator_CFRD Vs PostModulator_IGT",
                         output_dir_path = "de_results_2024/proteomics/postmod/p/",
                         plot_file_name = "PostModulator_CFRDVsPostModulator_IGT.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         molecule_names_file_path = "data/proteomics/protein_names.csv",
                         molecule_names_file_columns = c(4, 2),
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PostModulator_CFRD Vs PostModulator_IGT",
                         output_dir_path = "de_results_2024/proteomics/postmod/padj/",
                         plot_file_name = "PostModulator_CFRDVsPostModulator_IGT.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         molecule_names_file_path = "data/proteomics/protein_names.csv",
                         molecule_names_file_columns = c(4, 2),
                         plot_width_cm = 25)


contr_matrix <- makeContrasts(contrasts = "PostModulator_CFRD - PostModulator_NGT",
                              levels = colnames(model_matrix))
fit <- lmFit(data, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("protein")
result <- top.table %>%
  dplyr::select(protein, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = protein, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)

plot_volcano_and_save_DE(result, plot_title = "PostModulator_CFRD Vs PostModulator_NGT",
                         output_dir_path = "de_results_2024/proteomics/postmod/p/",
                         plot_file_name = "PostModulator_CFRDVsPostModulator_NGT.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         molecule_names_file_path = "data/proteomics/protein_names.csv",
                         molecule_names_file_columns = c(4, 2),
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PostModulator_CFRD Vs PostModulator_NGT",
                         output_dir_path = "de_results_2024/proteomics/postmod/padj/",
                         plot_file_name = "PostModulator_CFRDVsPostModulator_NGT.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         molecule_names_file_path = "data/proteomics/protein_names.csv",
                         molecule_names_file_columns = c(4, 2),
                         plot_width_cm = 25)


contr_matrix <- makeContrasts(contrasts = "PostModulator_IGT - PostModulator_NGT",
                              levels = colnames(model_matrix))
fit <- lmFit(data, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("protein")
result <- top.table %>%
  dplyr::select(protein, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = protein, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)

plot_volcano_and_save_DE(result, plot_title = "PostModulator_IGT Vs PostModulator_NGT",
                         output_dir_path = "de_results_2024/proteomics/postmod/p/",
                         plot_file_name = "PostModulator_IGTVsPostModulator_NGT.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         molecule_names_file_path = "data/proteomics/protein_names.csv",
                         molecule_names_file_columns = c(4, 2),
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PostModulator_IGT Vs PostModulator_NGT",
                         output_dir_path = "de_results_2024/proteomics/postmod/padj/",
                         plot_file_name = "PostModulator_IGTVsPostModulator_NGT.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         molecule_names_file_path = "data/proteomics/protein_names.csv",
                         molecule_names_file_columns = c(4, 2),
                         plot_width_cm = 25)

################################################################################################


contr_matrix <- makeContrasts(contrasts = "PostModulator_CFRD - PreModulator_CFRD",
                              levels = colnames(model_matrix))
fit <- lmFit(data, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("protein")
result <- top.table %>%
  dplyr::select(protein, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = protein, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)

plot_volcano_and_save_DE(result, plot_title = "PostModulator_CFRD Vs PreModulator_CFRD",
                         output_dir_path = "de_results_2024/proteomics/postmod_premod/p/",
                         plot_file_name = "PostModulator_CFRDVsPreModulator_CFRD.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         molecule_names_file_path = "data/proteomics/protein_names.csv",
                         molecule_names_file_columns = c(4, 2),
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PostModulator_CFRD Vs PreModulator_CFRD",
                         output_dir_path = "de_results_2024/proteomics/postmod_premod/padj/",
                         plot_file_name = "PostModulator_CFRDVsPreModulator_CFRD.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         molecule_names_file_path = "data/proteomics/protein_names.csv",
                         molecule_names_file_columns = c(4, 2),
                         plot_width_cm = 25)


contr_matrix <- makeContrasts(contrasts = "PostModulator_NGT - PreModulator_NGT",
                              levels = colnames(model_matrix))
fit <- lmFit(data, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("protein")
result <- top.table %>%
  dplyr::select(protein, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = protein, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)

plot_volcano_and_save_DE(result, plot_title = "PostModulator_NGT Vs PreModulator_NGT",
                         output_dir_path = "de_results_2024/proteomics/postmod_premod/p/",
                         plot_file_name = "PostModulator_NGTVsPreModulator_NGT.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         molecule_names_file_path = "data/proteomics/protein_names.csv",
                         molecule_names_file_columns = c(4, 2),
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PostModulator_NGT Vs PreModulator_NGT",
                         output_dir_path = "de_results_2024/proteomics/postmod_premod/padj/",
                         plot_file_name = "PostModulator_NGTVsPreModulator_NGT.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         molecule_names_file_path = "data/proteomics/protein_names.csv",
                         molecule_names_file_columns = c(4, 2),
                         plot_width_cm = 25)


contr_matrix <- makeContrasts(contrasts = "PostModulator_IGT - PreModulator_IGT",
                              levels = colnames(model_matrix))
fit <- lmFit(data, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("protein")
result <- top.table %>%
  dplyr::select(protein, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = protein, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)

plot_volcano_and_save_DE(result, plot_title = "PostModulator_IGT Vs PreModulator_IGT",
                         output_dir_path = "de_results_2024/proteomics/postmod_premod/p/",
                         plot_file_name = "PostModulator_IGTVsPreModulator_IGT.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         molecule_names_file_path = "data/proteomics/protein_names.csv",
                         molecule_names_file_columns = c(4, 2),
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PostModulator_IGT Vs PreModulator_IGT",
                         output_dir_path = "de_results_2024/proteomics/postmod_premod/padj/",
                         plot_file_name = "PostModulator_IGTVsPreModulator_IGT.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         molecule_names_file_path = "data/proteomics/protein_names.csv",
                         molecule_names_file_columns = c(4, 2),
                         plot_width_cm = 25)



model_matrix <- model.matrix(~ 0 + phenotype$PreModulatorVsPostModulator)
colnames(model_matrix) <- sub("phenotype$PreModulatorVsPostModulator", "", colnames(model_matrix), fixed = TRUE)

contr_matrix <- makeContrasts(contrasts = "PostModulator - PreModulator",
                              levels = colnames(model_matrix))
fit <- lmFit(data, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("protein")
result <- top.table %>%
  dplyr::select(protein, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = protein, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)

plot_volcano_and_save_DE(result, plot_title = "PostModulator Vs PreModulator",
                         output_dir_path = "de_results_2024/proteomics/postmod_premod/p/",
                         plot_file_name = "PostModulatorVsPreModulator.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         molecule_names_file_path = "data/proteomics/protein_names.csv",
                         molecule_names_file_columns = c(4, 2),
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PostModulator Vs PreModulator",
                         output_dir_path = "de_results_2024/proteomics/postmod_premod/padj/",
                         plot_file_name = "PostModulatorVsPreModulator.png",
                         fc_cutoff = 1.5,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         molecule_names_file_path = "data/proteomics/protein_names.csv",
                         molecule_names_file_columns = c(4, 2),
                         plot_width_cm = 25)


####################################
#create boxplots of upreg and downreg proteins in premod cfrd vs ngt, across premod cfrd, premod ngt, postmod cfrd

create_DE_boxplot(data, phenotype, 
                  conditions_of_interest = c("PreModulator_CFRD", "PreModulator_NGT", "PostModulator_CFRD"),
                  x_lab = "Proteins", output_dir_path = "plots_updated/post_mod/shift_cfrd_to_ngt/",
                  de_file_path = "de_results_2024/proteomics/premod/p/sig_no_name_PreModulator_CFRDVsPreModulator_NGT.csv", 
                  k = 10)

####################################
#venn diagram of post vs pre overlap DE

de.CFRD <- read.table("de_results_2024/proteomics/postmod_premod/p/sig_no_name_PostModulator_CFRDVsPreModulator_CFRD.csv", 
                      sep = "\t", header = TRUE)  
de.CFRD.up <- de.CFRD %>%
  filter(logFC > 0)
de.CFRD.down <- de.CFRD %>%
  filter(logFC < 0)

de.IGT <- read.table("de_results_2024/proteomics/postmod_premod/p/sig_no_name_PostModulator_IGTVsPreModulator_IGT.csv", 
                     sep = "\t", header = TRUE)  
de.IGT.up <- de.IGT %>%
  filter(logFC > 0)
de.IGT.down <- de.IGT %>%
  filter(logFC < 0)

de.NGT <- read.table("de_results_2024/proteomics/postmod_premod/p/sig_no_name_PostModulator_NGTVsPreModulator_NGT.csv", 
                     sep = "\t", header = TRUE)  
de.NGT.up <- de.NGT %>%
  filter(logFC > 0)
de.NGT.down <- de.NGT %>%
  filter(logFC < 0)

ggvenn(list("CFRD" = de.CFRD.up$Molecule,
            "IGT" = de.IGT.up$Molecule,
            "NGT" = de.NGT.up$Molecule),
       stroke_size = 0.1,
       set_name_size = 4,
       text_size = 3,
       fill_alpha = 0.5,
       fill_color = c("red", "orange", "yellow")) +
  ggtitle("PostModulator Vs PreModulator upregulated protein overlap") +
  theme(plot.title = element_text(vjust = 0, hjust = 0.5, size = rel(1.2), face = "bold"))
ggsave("de_results_2024/proteomics/postmod_premod/overlap_up.png")

intersect(de.CFRD.up$Molecule, de.IGT.up$Molecule)
intersect(de.IGT.up$Molecule, de.NGT.up$Molecule)

ggvenn(list("CFRD" = de.CFRD.down$Molecule,
            "IGT" = de.IGT.down$Molecule,
            "NGT" = de.NGT.down$Molecule),
       stroke_size = 0.1,
       set_name_size = 4,
       text_size = 3,
       fill_alpha = 0.5,
       fill_color = c("red", "orange", "yellow")) +
  ggtitle("PostModulator Vs PreModulator downregulated protein overlap") +
  theme(plot.title = element_text(vjust = 0, hjust = 0.5, size = rel(1.2), face = "bold"))
ggsave("de_results_2024/proteomics/postmod_premod/overlap_down.png")

intersect(de.CFRD.down$Molecule, de.IGT.down$Molecule)
intersect(de.CFRD.down$Molecule, de.NGT.down$Molecule)
intersect(de.IGT.down$Molecule, de.NGT.down$Molecule)


####################################
#venn diagram of pre CFRD, IGT, NGT DE results overlap with identified biomarkers

biomarkers.CFRD_IGT <- read.xlsx("data/selected_features/prot_combat_biomarkers.xlsx", sheetName = "prot_CFRDVsIGT")
biomarkers.CFRD_NGT <- read.xlsx("data/selected_features/prot_combat_biomarkers.xlsx", sheetName = "prot_CFRDVsNGT")
biomarkers.IGT_NGT <- read.xlsx("data/selected_features/prot_combat_biomarkers.xlsx", sheetName = "prot_IGTVsNGT")

de.CFRD_IGT <- read.table("de_results_2024/proteomics/premod/p/sig_no_name_PreModulator_CFRDVsPreModulator_IGT.csv", 
                          sep = "\t", header = TRUE)  
de.CFRD_IGT.up <- de.CFRD_IGT %>%
  filter(logFC > 0)
de.CFRD_IGT.down <- de.CFRD_IGT %>%
  filter(logFC < 0)

common <- intersect(de.CFRD_IGT.up$Molecule, biomarkers.CFRD_IGT$biomarkers)
if(length(common) > 0){
  caption_text <- paste0("Common: ", paste(common, collapse = ", ", sep = ""))
} else{
  caption_text <- ""
}
ggvenn(list("Upregulated proteins" = de.CFRD_IGT.up$Molecule,
            "Identified biomarkers" = biomarkers.CFRD_IGT$biomarkers),
       stroke_size = 0.1,
       set_name_size = 4,
       text_size = 3,
       fill_alpha = 0.5,
       fill_color = c("red", "gold")) +
  ggtitle("CFRD Vs IGT proteomics") +
  labs(caption = caption_text) +
  theme(plot.title = element_text(vjust = 0, hjust = 0.5, size = rel(1.2), face = "bold"),
        plot.caption = element_text(hjust = 0.5))
ggsave("biomarker_de_overlap/protein_CFRD_IGT_up.png")

common <- intersect(de.CFRD_IGT.down$Molecule, biomarkers.CFRD_IGT$biomarkers)
if(length(common) > 0){
  caption_text <- paste0("Common: ", paste(common, collapse = ", ", sep = ""))
} else{
  caption_text <- ""
}
ggvenn(list("Downregulated proteins" = de.CFRD_IGT.down$Molecule,
            "Identified biomarkers" = biomarkers.CFRD_IGT$biomarkers),
       stroke_size = 0.1,
       set_name_size = 4,
       text_size = 3,
       fill_alpha = 0.5,
       fill_color = c("blue", "gold")) +
  ggtitle("CFRD Vs IGT proteomics") +
  labs(caption = caption_text) +
  theme(plot.title = element_text(vjust = 0, hjust = 0.5, size = rel(1.2), face = "bold"),
        plot.caption = element_text(hjust = 0.5))
ggsave("biomarker_de_overlap/protein_CFRD_IGT_down.png")


de.CFRD_NGT <- read.table("de_results_2024/proteomics/premod/p/sig_no_name_PreModulator_CFRDVsPreModulator_NGT.csv", 
                          sep = "\t", header = TRUE)  
de.CFRD_NGT.up <- de.CFRD_NGT %>%
  filter(logFC > 0)
de.CFRD_NGT.down <- de.CFRD_NGT %>%
  filter(logFC < 0)

common <- intersect(de.CFRD_NGT.up$Molecule, biomarkers.CFRD_NGT$biomarkers)
if(length(common) > 0){
  caption_text <- paste0("Common: ", paste(common, collapse = ", ", sep = ""))
} else{
  caption_text <- ""
}
ggvenn(list("Upregulated proteins" = de.CFRD_NGT.up$Molecule,
            "Identified biomarkers" = biomarkers.CFRD_NGT$biomarkers),
       stroke_size = 0.1,
       set_name_size = 4,
       text_size = 3,
       fill_alpha = 0.5,
       fill_color = c("red", "gold")) +
  ggtitle("CFRD Vs NGT proteomics") +
  labs(caption = caption_text) +
  theme(plot.title = element_text(vjust = 0, hjust = 0.5, size = rel(1.2), face = "bold"),
        plot.caption = element_text(hjust = 0.5))
ggsave("biomarker_de_overlap/protein_CFRD_NGT_up.png")

common <- intersect(de.CFRD_NGT.down$Molecule, biomarkers.CFRD_NGT$biomarkers)
if(length(common) > 0){
  caption_text <- paste0("Common: ", paste(common, collapse = ", ", sep = ""))
} else{
  caption_text <- ""
}
ggvenn(list("Downregulated proteins" = de.CFRD_NGT.down$Molecule,
            "Identified biomarkers" = biomarkers.CFRD_NGT$biomarkers),
       stroke_size = 0.1,
       set_name_size = 4,
       text_size = 3,
       fill_alpha = 0.5,
       fill_color = c("blue", "gold")) +
  ggtitle("CFRD Vs NGT proteomics") +
  labs(caption = caption_text) +
  theme(plot.title = element_text(vjust = 0, hjust = 0.5, size = rel(1.2), face = "bold"),
        plot.caption = element_text(hjust = 0.5))
ggsave("biomarker_de_overlap/protein_CFRD_NGT_down.png")


de.IGT_NGT <- read.table("de_results_2024/proteomics/premod/p/sig_no_name_PreModulator_IGTVsPreModulator_NGT.csv", 
                         sep = "\t", header = TRUE)  
de.IGT_NGT.up <- de.IGT_NGT %>%
  filter(logFC > 0)
de.IGT_NGT.down <- de.IGT_NGT %>%
  filter(logFC < 0)

common <- intersect(de.IGT_NGT.up$Molecule, biomarkers.IGT_NGT$biomarkers)
if(length(common) > 0){
  caption_text <- paste0("Common: ", paste(common, collapse = ", ", sep = ""))
} else{
  caption_text <- ""
}
ggvenn(list("Upregulated proteins" = de.IGT_NGT.up$Molecule,
            "Identified biomarkers" = biomarkers.IGT_NGT$biomarkers),
       stroke_size = 0.1,
       set_name_size = 4,
       text_size = 3,
       fill_alpha = 0.5,
       fill_color = c("red", "gold")) +
  ggtitle("IGT Vs NGT proteomics") +
  labs(caption = caption_text) +
  theme(plot.title = element_text(vjust = 0, hjust = 0.5, size = rel(1.2), face = "bold"),
        plot.caption = element_text(hjust = 0.5))
ggsave("biomarker_de_overlap/protein_IGT_NGT_up.png")

common <- intersect(de.IGT_NGT.down$Molecule, biomarkers.IGT_NGT$biomarkers)
if(length(common) > 0){
  caption_text <- paste0("Common: ", paste(common, collapse = ", ", sep = ""))
} else{
  caption_text <- ""
}
ggvenn(list("Downregulated proteins" = de.IGT_NGT.down$Molecule,
            "Identified biomarkers" = biomarkers.IGT_NGT$biomarkers),
       stroke_size = 0.1,
       set_name_size = 4,
       text_size = 3,
       fill_alpha = 0.5,
       fill_color = c("blue", "gold")) +
  ggtitle("IGT Vs NGT proteomics") +
  labs(caption = caption_text) +
  theme(plot.title = element_text(vjust = 0, hjust = 0.5, size = rel(1.2), face = "bold"),
        plot.caption = element_text(hjust = 0.5))
ggsave("biomarker_de_overlap/protein_IGT_NGT_down.png")


#######################################################################

#venn diagram within premod overlap with within postmod

create_DE_up_down_venn(de_file_path1 = "de_results_2024/proteomics/premod/p/sig_no_name_PreModulator_CFRDVsPreModulator_IGT.csv",
                       de_file_path2 = "de_results_2024/proteomics/postmod/p/sig_no_name_PostModulator_CFRDVsPostModulator_IGT.csv",
                       output_dir_path = "de_results_2024/proteomics/pre_post_overlap/",
                       file_name_upreg = "proteomics_pre_post_CFRD_IGT_up.png",
                       file_name_downreg = "proteomics_pre_post_CFRD_IGT_down.png",
                       comparison1_name = "Pre Modulator",
                       comparison2_name = "Post Modulator",
                       title <- "CFRD Vs IGT")
create_DE_up_down_venn(de_file_path1 = "de_results_2024/proteomics/premod/p/sig_no_name_PreModulator_CFRDVsPreModulator_NGT.csv",
                       de_file_path2 = "de_results_2024/proteomics/postmod/p/sig_no_name_PostModulator_CFRDVsPostModulator_NGT.csv",
                       output_dir_path = "de_results_2024/proteomics/pre_post_overlap/",
                       file_name_upreg = "proteomics_pre_post_CFRD_NGT_up.png",
                       file_name_downreg = "proteomics_pre_post_CFRD_NGT_down.png",
                       comparison1_name = "Pre Modulator",
                       comparison2_name = "Post Modulator",
                       title <- "CFRD Vs NGT")
create_DE_up_down_venn(de_file_path1 = "de_results_2024/proteomics/premod/p/sig_no_name_PreModulator_IGTVsPreModulator_NGT.csv",
                       de_file_path2 = "de_results_2024/proteomics/postmod/p/sig_no_name_PostModulator_IGTVsPostModulator_NGT.csv",
                       output_dir_path = "de_results_2024/proteomics/pre_post_overlap/",
                       file_name_upreg = "proteomics_pre_post_IGT_NGT_up.png",
                       file_name_downreg = "proteomics_pre_post_IGT_NGT_down.png",
                       comparison1_name = "Pre Modulator",
                       comparison2_name = "Post Modulator",
                       title <- "IGT Vs NGT")
