library(sva)

base_dir <- "~/UNSW/VafaeeLab/CysticFibrosisGroup/ExoCF/CFRD_EV_biomarker/"
setwd(base_dir)

source("utils_diff.R")
source("utils.R")


data <- read.csv("data/formatted/rna_all/umi_counts_filter90.csv", row.names = 1)

phenotype <- read.table("data/formatted/tra_phenotype_2024Jan.txt", header = TRUE, sep = "\t") %>%
  mutate(modstatus_condition_country = paste(PreModulatorVsPostModulator, condition, country, sep = "_")) %>%
  mutate(modstatus_condition = paste(PreModulatorVsPostModulator, condition, sep = "_"))

summary(factor(phenotype$condition))
# CFRD                 HC                IGT                NGT UNKNOWN TO PREDICT 
# 82                  6                 80                138                 28


# %>%
#   dplyr::rename(c("disease_status" = "condition"))

phenotype <- phenotype %>%
  filter(condition != "UNKNOWN TO PREDICT" & condition != "HC")
summary(factor(phenotype$modstatus_condition_country))
summary(factor(phenotype$modstatus_condition))

data <- data[, phenotype$Sample]

all.equal(colnames(data), phenotype$Sample)

#filter
keep <- edgeR::filterByExpr(data, group = phenotype$modstatus_condition)
data <- data[keep, ]

#normalize
data <- edgeR::cpm(data, log = TRUE) 

all.equal(colnames(data), phenotype$Sample)
#TRUE

phenotype <- phenotype %>%
  mutate(country_batch_name = paste(country, batch_name, sep = "_"))
summary(factor(phenotype$country_batch_name))

# AU_initial DK_initial     DK_new 
# 105        133         62 



#########################################################

create_dim_red_plots(comparison = NA,
                     classes = c("CFRD", "IGT"), 
                     class_colours = c("red", "orange"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics",
                     perform_filter = FALSE,
                     colour_column = "PreModulatorVsPostModulator", point_border_colours = c("black", "green"),
                     data = data,
                     phenotype = phenotype,
                     filter_post_modulator = FALSE,
                     custom_title = "1 Transcriptomics CFRD, IGT pre and post modulator samples with LogCPM norm", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("CFRD", "NGT"), 
                     class_colours = c("red", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics",
                     perform_filter = FALSE,
                     colour_column = "PreModulatorVsPostModulator", point_border_colours = c("black", "green"),
                     data = data,
                     phenotype = phenotype,
                     filter_post_modulator = FALSE,
                     custom_title = "2 Transcriptomics CFRD, NGT pre and post modulator samples with LogCPM norm", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("IGT", "NGT"), 
                     class_colours = c("orange", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics",
                     perform_filter = FALSE,
                     colour_column = "PreModulatorVsPostModulator", point_border_colours = c("black", "green"),
                     data = data,
                     phenotype = phenotype,
                     filter_post_modulator = FALSE,
                     custom_title = "3 Transcriptomics IGT, NGT pre and post modulator samples with LogCPM norm", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PostModulator_CFRD", "PreModulator_NGT"), 
                     class_colours = c("red", "indianred", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics",
                     dimred_plot_width_cm = 36,
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "4 Transcriptomics shift from CFRD to NGT with LogCPM norm", 
                     combat = FALSE)

create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PreModulator_IGT"), 
                     class_colours = c("red", "orange"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "5 Transcriptomics premod CFRD IGT with LogCPM norm", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PreModulator_NGT"), 
                     class_colours = c("red", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "6 Transcriptomics premod CFRD NGT with LogCPM norm", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_IGT", "PreModulator_NGT"), 
                     class_colours = c("orange", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "7 Transcriptomics premod IGT NGT with LogCPM norm", 
                     combat = FALSE)

create_dim_red_plots(comparison = NA,
                     classes = c("PostModulator_CFRD", "PostModulator_IGT"), 
                     class_colours = c("red", "orange"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "8 Transcriptomics postmod CFRD IGT with LogCPM norm", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PostModulator_CFRD", "PostModulator_NGT"), 
                     class_colours = c("red", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "9 Transcriptomics postmod CFRD NGT with LogCPM norm", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PostModulator_IGT", "PostModulator_NGT"), 
                     class_colours = c("orange", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "10 Transcriptomics postmod IGT NGT with LogCPM norm", 
                     combat = FALSE)

#########################################################

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
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected",
                     perform_filter = FALSE,
                     colour_column = "PreModulatorVsPostModulator", point_border_colours = c("black", "green"),
                     data = data,
                     phenotype = phenotype,
                     filter_post_modulator = FALSE,
                     custom_title = "1 Transcriptomics CFRD, IGT pre and post modulator samples with Log CPM norm and ComBat", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("CFRD", "NGT"), 
                     class_colours = c("red", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected",
                     perform_filter = FALSE,
                     colour_column = "PreModulatorVsPostModulator", point_border_colours = c("black", "green"),
                     data = data,
                     phenotype = phenotype,
                     filter_post_modulator = FALSE,
                     custom_title = "2 Transcriptomics CFRD, NGT pre and post modulator samples with Log CPM norm and ComBat", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("IGT", "NGT"), 
                     class_colours = c("orange", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected",
                     perform_filter = FALSE,
                     colour_column = "PreModulatorVsPostModulator", point_border_colours = c("black", "green"),
                     data = data,
                     phenotype = phenotype,
                     filter_post_modulator = FALSE,
                     custom_title = "3 Transcriptomics IGT, NGT pre and post modulator samples with Log CPM norm and ComBat", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PostModulator_CFRD", "PreModulator_NGT"), 
                     class_colours = c("red", "indianred", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected",
                     dimred_plot_width_cm = 36,
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "4 Transcriptomics shift from CFRD to NGT with Log CPM norm and ComBat", 
                     combat = FALSE)

create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PreModulator_IGT"), 
                     class_colours = c("red", "orange"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "5 Transcriptomics premod CFRD IGT with Log CPM norm and ComBat", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PreModulator_NGT"), 
                     class_colours = c("red", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "6 Transcriptomics premod CFRD NGT with Log CPM norm and ComBat", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_IGT", "PreModulator_NGT"), 
                     class_colours = c("orange", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "7 Transcriptomics premod IGT NGT with Log CPM norm and ComBat", 
                     combat = FALSE)

create_dim_red_plots(comparison = NA,
                     classes = c("PostModulator_CFRD", "PostModulator_IGT"), 
                     class_colours = c("red", "orange"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "8 Transcriptomics postmod CFRD IGT with Log CPM norm and ComBat", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PostModulator_CFRD", "PostModulator_NGT"), 
                     class_colours = c("red", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "9 Transcriptomics postmod CFRD NGT with Log CPM norm and ComBat", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PostModulator_IGT", "PostModulator_NGT"), 
                     class_colours = c("orange", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "10 Transcriptomics postmod IGT NGT with Log CPM norm and ComBat", 
                     combat = FALSE)

#########################################################

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
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected_twice",
                     perform_filter = FALSE,
                     colour_column = "PreModulatorVsPostModulator", point_border_colours = c("black", "green"),
                     data = data,
                     phenotype = phenotype,
                     filter_post_modulator = FALSE,
                     custom_title = "1 Transcriptomics CFRD, IGT pre and post modulator samples with Log CPM norm and ComBat twice", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("CFRD", "NGT"), 
                     class_colours = c("red", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected_twice",
                     perform_filter = FALSE,
                     colour_column = "PreModulatorVsPostModulator", point_border_colours = c("black", "green"),
                     data = data,
                     phenotype = phenotype,
                     filter_post_modulator = FALSE,
                     custom_title = "2 Transcriptomics CFRD, NGT pre and post modulator samples with Log CPM norm and ComBat twice", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("IGT", "NGT"), 
                     class_colours = c("orange", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected_twice",
                     perform_filter = FALSE,
                     colour_column = "PreModulatorVsPostModulator", point_border_colours = c("black", "green"),
                     data = data,
                     phenotype = phenotype,
                     filter_post_modulator = FALSE,
                     custom_title = "3 Transcriptomics IGT, NGT pre and post modulator samples with Log CPM norm and ComBat twice", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PostModulator_CFRD", "PreModulator_NGT"), 
                     class_colours = c("red", "indianred", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected_twice",
                     dimred_plot_width_cm = 36,
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "4 Transcriptomics shift from CFRD to NGT with Log CPM norm and ComBat twice", 
                     combat = FALSE)

create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PreModulator_IGT"), 
                     class_colours = c("red", "orange"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected_twice",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "5 Transcriptomics premod CFRD IGT with Log CPM norm and ComBat twice", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PreModulator_NGT"), 
                     class_colours = c("red", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected_twice",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "6 Transcriptomics premod CFRD NGT with Log CPM norm and ComBat twice", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_IGT", "PreModulator_NGT"), 
                     class_colours = c("orange", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected_twice",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "7 Transcriptomics premod IGT NGT with Log CPM norm and ComBat twice", 
                     combat = FALSE)

create_dim_red_plots(comparison = NA,
                     classes = c("PostModulator_CFRD", "PostModulator_IGT"), 
                     class_colours = c("red", "orange"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected_twice",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "8 Transcriptomics postmod CFRD IGT with Log CPM norm and ComBat twice", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PostModulator_CFRD", "PostModulator_NGT"), 
                     class_colours = c("red", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected_twice",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "9 Transcriptomics postmod CFRD NGT with Log CPM norm and ComBat twice", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PostModulator_IGT", "PostModulator_NGT"), 
                     class_colours = c("orange", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected_twice",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "10 Transcriptomics postmod IGT NGT with Log CPM norm and ComBat twice", 
                     combat = FALSE)


#########################################################

# combat with country+batch_name

data <- read.csv("data/formatted/rna_all/umi_counts_filter90.csv", row.names = 1)

phenotype <- read.table("data/formatted/tra_phenotype_2024Jan.txt", header = TRUE, sep = "\t") %>%
  mutate(modstatus_condition_country = paste(PreModulatorVsPostModulator, condition, country, sep = "_")) %>%
  mutate(modstatus_condition = paste(PreModulatorVsPostModulator, condition, sep = "_"))

summary(factor(phenotype$condition))
# CFRD                 HC                IGT                NGT UNKNOWN TO PREDICT 
# 82                  6                 80                138                 28


# %>%
#   dplyr::rename(c("disease_status" = "condition"))

phenotype <- phenotype %>%
  filter(condition != "UNKNOWN TO PREDICT" & condition != "HC")
summary(factor(phenotype$modstatus_condition_country))
summary(factor(phenotype$modstatus_condition))

data <- data[, phenotype$Sample]

all.equal(colnames(data), phenotype$Sample)

#filter
keep <- edgeR::filterByExpr(data, group = phenotype$modstatus_condition)
data <- data[keep, ]

#normalize
data <- edgeR::cpm(data, log = TRUE) 

all.equal(colnames(data), phenotype$Sample)
#TRUE

phenotype <- phenotype %>%
  mutate(country_batch_name = paste(country, batch_name, sep = "_"))
summary(factor(phenotype$country_batch_name))

# AU_initial DK_initial     DK_new 
# 105        133         62

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
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected_two_in_one_shot",
                     perform_filter = FALSE,
                     colour_column = "PreModulatorVsPostModulator", point_border_colours = c("black", "green"),
                     data = data,
                     phenotype = phenotype,
                     filter_post_modulator = FALSE,
                     custom_title = "1 Transcriptomics CFRD, IGT pre and post modulator samples with Log CPM norm and ComBat two_in_one_shot", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("CFRD", "NGT"), 
                     class_colours = c("red", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected_two_in_one_shot",
                     perform_filter = FALSE,
                     colour_column = "PreModulatorVsPostModulator", point_border_colours = c("black", "green"),
                     data = data,
                     phenotype = phenotype,
                     filter_post_modulator = FALSE,
                     custom_title = "2 Transcriptomics CFRD, NGT pre and post modulator samples with Log CPM norm and ComBat two_in_one_shot", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("IGT", "NGT"), 
                     class_colours = c("orange", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected_two_in_one_shot",
                     perform_filter = FALSE,
                     colour_column = "PreModulatorVsPostModulator", point_border_colours = c("black", "green"),
                     data = data,
                     phenotype = phenotype,
                     filter_post_modulator = FALSE,
                     custom_title = "3 Transcriptomics IGT, NGT pre and post modulator samples with Log CPM norm and ComBat two_in_one_shot", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PostModulator_CFRD", "PreModulator_NGT"), 
                     class_colours = c("red", "indianred", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected_two_in_one_shot",
                     dimred_plot_width_cm = 36,
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "4 Transcriptomics shift from CFRD to NGT with Log CPM norm and ComBat two_in_one_shot", 
                     combat = FALSE)

create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PostModulator_CFRD", "PreModulator_NGT"), 
                     class_colours = c("red", "indianred", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/shift",
                     dimred_plot_width_cm = 36,
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "Transcriptomics shift from CFRD to NGT with Log CPM norm and ComBat two_in_one_shot", 
                     combat = FALSE, simplified = TRUE)

create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PreModulator_IGT"), 
                     class_colours = c("red", "orange"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected_two_in_one_shot",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "5 Transcriptomics premod CFRD IGT with Log CPM norm and ComBat two_in_one_shot", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_CFRD", "PreModulator_NGT"), 
                     class_colours = c("red", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected_two_in_one_shot",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "6 Transcriptomics premod CFRD NGT with Log CPM norm and ComBat two_in_one_shot", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PreModulator_IGT", "PreModulator_NGT"), 
                     class_colours = c("orange", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected_two_in_one_shot",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "7 Transcriptomics premod IGT NGT with Log CPM norm and ComBat two_in_one_shot", 
                     combat = FALSE)

create_dim_red_plots(comparison = NA,
                     classes = c("PostModulator_CFRD", "PostModulator_IGT"), 
                     class_colours = c("red", "orange"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected_two_in_one_shot",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "8 Transcriptomics postmod CFRD IGT with Log CPM norm and ComBat two_in_one_shot", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PostModulator_CFRD", "PostModulator_NGT"), 
                     class_colours = c("red", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected_two_in_one_shot",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "9 Transcriptomics postmod CFRD NGT with Log CPM norm and ComBat two_in_one_shot", 
                     combat = FALSE)
create_dim_red_plots(comparison = NA,
                     classes = c("PostModulator_IGT", "PostModulator_NGT"), 
                     class_colours = c("orange", "yellow"),
                     dim_red = "UMAP",
                     norm = "non-normalized",
                     dir_path = "plots_updated/post_mod/transcriptomics/batch_corrected_two_in_one_shot",
                     perform_filter = FALSE,
                     colour_column = "batch_name", point_border_colours = c("black", "green", "purple"),
                     data = data,
                     phenotype = phenotype %>% dplyr::rename("disease_status" = "condition") %>% dplyr::rename("condition" = "modstatus_condition"),
                     filter_post_modulator = FALSE,
                     custom_title = "10 Transcriptomics postmod IGT NGT with Log CPM norm and ComBat two_in_one_shot", 
                     combat = FALSE)


#########################################################


#use this combat 2 in 1 shot for DE analysis

all.equal(phenotype$Sample, colnames(data))
# TRUE

summary(factor(phenotype$condition))
# CFRD  IGT  NGT 
# 82   80  138 

summary(factor(phenotype$modstatus_condition))
# PostModulator_CFRD  PostModulator_IGT  PostModulator_NGT  PreModulator_CFRD   PreModulator_IGT   PreModulator_NGT 
# 48                 36                 76                 34                 44                 62 

model_matrix <- model.matrix(~ 0 + phenotype$modstatus_condition)
colnames(model_matrix) <- sub("phenotype$modstatus_condition", "", colnames(model_matrix), fixed = TRUE)


contr_matrix <- makeContrasts(contrasts = "PreModulator_CFRD - PreModulator_IGT",
                              levels = colnames(model_matrix))

# v <- voom(data, model_matrix, plot = TRUE)
# Error in voom(data, model_matrix, plot = TRUE) : 
#   Negative counts not allowed

fit <- lmFit(data, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("rna")
result <- top.table %>%
  dplyr::select(rna, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = rna, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)

plot_volcano_and_save_DE(result, plot_title = "PreModulator_CFRD Vs PreModulator_IGT",
                         output_dir_path = "de_results_2024/transcriptomics/premod/p/",
                         plot_file_name = "PreModulator_CFRDVsPreModulator_IGT.png",
                         fc_cutoff = 1.2,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PreModulator_CFRD Vs PreModulator_IGT",
                         output_dir_path = "de_results_2024/transcriptomics/premod/padj/",
                         plot_file_name = "PreModulator_CFRDVsPreModulator_IGT.png",
                         fc_cutoff = 1.2,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         plot_width_cm = 25)


contr_matrix <- makeContrasts(contrasts = "PreModulator_CFRD - PreModulator_NGT",
                              levels = colnames(model_matrix))
fit <- lmFit(data, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("rna")
result <- top.table %>%
  dplyr::select(rna, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = rna, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)

plot_volcano_and_save_DE(result, plot_title = "PreModulator_CFRD Vs PreModulator_NGT",
                         output_dir_path = "de_results_2024/transcriptomics/premod/p/",
                         plot_file_name = "PreModulator_CFRDVsPreModulator_NGT.png",
                         fc_cutoff = 1.2,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PreModulator_CFRD Vs PreModulator_NGT",
                         output_dir_path = "de_results_2024/transcriptomics/premod/padj/",
                         plot_file_name = "PreModulator_CFRDVsPreModulator_NGT.png",
                         fc_cutoff = 1.2,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         plot_width_cm = 25)


contr_matrix <- makeContrasts(contrasts = "PreModulator_IGT - PreModulator_NGT",
                              levels = colnames(model_matrix))
fit <- lmFit(data, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("rna")
result <- top.table %>%
  dplyr::select(rna, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = rna, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)

plot_volcano_and_save_DE(result, plot_title = "PreModulator_IGT Vs PreModulator_NGT",
                         output_dir_path = "de_results_2024/transcriptomics/premod/p/",
                         plot_file_name = "PreModulator_IGTVsPreModulator_NGT.png",
                         fc_cutoff = 1.2,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PreModulator_IGT Vs PreModulator_NGT",
                         output_dir_path = "de_results_2024/transcriptomics/premod/padj/",
                         plot_file_name = "PreModulator_IGTVsPreModulator_NGT.png",
                         fc_cutoff = 1.2,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
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
  rownames_to_column("rna")
result <- top.table %>%
  dplyr::select(rna, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = rna, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)

plot_volcano_and_save_DE(result, plot_title = "PostModulator_CFRD Vs PreModulator_CFRD",
                         output_dir_path = "de_results_2024/transcriptomics/postmod_premod/p/",
                         plot_file_name = "PostModulator_CFRDVsPreModulator_CFRD.png",
                         fc_cutoff = 1.2,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PostModulator_CFRD Vs PreModulator_CFRD",
                         output_dir_path = "de_results_2024/transcriptomics/postmod_premod/padj/",
                         plot_file_name = "PostModulator_CFRDVsPreModulator_CFRD.png",
                         fc_cutoff = 1.2,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         plot_width_cm = 25)


contr_matrix <- makeContrasts(contrasts = "PostModulator_NGT - PreModulator_NGT",
                              levels = colnames(model_matrix))
fit <- lmFit(data, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("rna")
result <- top.table %>%
  dplyr::select(rna, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = rna, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)

plot_volcano_and_save_DE(result, plot_title = "PostModulator_NGT Vs PreModulator_NGT",
                         output_dir_path = "de_results_2024/transcriptomics/postmod_premod/p/",
                         plot_file_name = "PostModulator_NGTVsPreModulator_NGT.png",
                         fc_cutoff = 1.2,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PostModulator_NGT Vs PreModulator_NGT",
                         output_dir_path = "de_results_2024/transcriptomics/postmod_premod/padj/",
                         plot_file_name = "PostModulator_NGTVsPreModulator_NGT.png",
                         fc_cutoff = 1.2,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         plot_width_cm = 25)


contr_matrix <- makeContrasts(contrasts = "PostModulator_IGT - PreModulator_IGT",
                              levels = colnames(model_matrix))
fit <- lmFit(data, model_matrix)
# head(coef(fit))
fit <- contrasts.fit(fit, contr_matrix)
# head(coef(fit))
efit <- eBayes(fit)
top.table <- topTable(efit, n = Inf, sort.by = "p") %>%
  rownames_to_column("rna")
result <- top.table %>%
  dplyr::select(rna, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = rna, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)

plot_volcano_and_save_DE(result, plot_title = "PostModulator_IGT Vs PreModulator_IGT",
                         output_dir_path = "de_results_2024/transcriptomics/postmod_premod/p/",
                         plot_file_name = "PostModulator_IGTVsPreModulator_IGT.png",
                         fc_cutoff = 1.2,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PostModulator_IGT Vs PreModulator_IGT",
                         output_dir_path = "de_results_2024/transcriptomics/postmod_premod/padj/",
                         plot_file_name = "PostModulator_IGTVsPreModulator_IGT.png",
                         fc_cutoff = 1.2,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
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
  rownames_to_column("rna")
result <- top.table %>%
  dplyr::select(rna, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(Molecule = rna, adjPVal = adj.P.Val, PVal = P.Value) %>%
  arrange(logFC)

plot_volcano_and_save_DE(result, plot_title = "PostModulator Vs PreModulator",
                         output_dir_path = "de_results_2024/transcriptomics/postmod_premod/p/",
                         plot_file_name = "PostModulatorVsPreModulator.png",
                         fc_cutoff = 1.2,
                         pval_cutoff = 0.05,
                         use_adj_pval = FALSE,
                         plot_width_cm = 25)
plot_volcano_and_save_DE(result, plot_title = "PostModulator Vs PreModulator",
                         output_dir_path = "de_results_2024/transcriptomics/postmod_premod/padj/",
                         plot_file_name = "PostModulatorVsPreModulator.png",
                         fc_cutoff = 1.2,
                         pval_cutoff = 0.05,
                         use_adj_pval = TRUE,
                         plot_width_cm = 25)

####################################
#create boxplots of upreg and downreg transcripts in premod cfrd vs ngt, across premod cfrd, premod ngt, postmod cfrd

create_DE_boxplot(data, phenotype, 
                  conditions_of_interest = c("PreModulator_CFRD", "PreModulator_NGT", "PostModulator_CFRD"),
                  x_lab = "Transcripts", output_dir_path = "plots_updated/post_mod/shift_cfrd_to_ngt/",
                  de_file_path = "de_results_2024/transcriptomics/premod/p/sig_PreModulator_CFRDVsPreModulator_NGT.csv", 
                  k = 5)
create_DE_boxplot(data, phenotype, 
                  conditions_of_interest = c("PreModulator_CFRD", "PreModulator_NGT", "PostModulator_CFRD"),
                  x_lab = "Transcripts", output_dir_path = "plots_updated/post_mod/shift_cfrd_to_ngt/",
                  de_file_path = "de_results_2024/transcriptomics/premod/p/sig_PreModulator_CFRDVsPreModulator_NGT.csv", 
                  k = 10)


####################################
#venn diagram of post vs pre overlap DE

de.CFRD <- read.table("de_results_2024/transcriptomics/postmod_premod/p/sig_PostModulator_CFRDVsPreModulator_CFRD.csv", 
                      sep = "\t", header = TRUE)  
de.CFRD.up <- de.CFRD %>%
  filter(logFC > 0)
de.CFRD.down <- de.CFRD %>%
  filter(logFC < 0)

de.IGT <- read.table("de_results_2024/transcriptomics/postmod_premod/p/sig_PostModulator_IGTVsPreModulator_IGT.csv", 
                     sep = "\t", header = TRUE)  
de.IGT.up <- de.IGT %>%
  filter(logFC > 0)
de.IGT.down <- de.IGT %>%
  filter(logFC < 0)

de.NGT <- read.table("de_results_2024/transcriptomics/postmod_premod/p/sig_PostModulator_NGTVsPreModulator_NGT.csv", 
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
  ggtitle("PostModulator Vs PreModulator upregulated transcript overlap") +
  theme(plot.title = element_text(vjust = 0, hjust = 0.5, size = rel(1.2), face = "bold"))
ggsave("de_results_2024/transcriptomics/postmod_premod/overlap_up.png")

intersect(de.CFRD.up$Molecule, de.IGT.up$Molecule)
intersect(de.CFRD.up$Molecule, de.NGT.up$Molecule)
intersect(de.IGT.up$Molecule, de.NGT.up$Molecule)

ggvenn(list("CFRD" = de.CFRD.down$Molecule,
            "IGT" = de.IGT.down$Molecule,
            "NGT" = de.NGT.down$Molecule),
       stroke_size = 0.1,
       set_name_size = 4,
       text_size = 3,
       fill_alpha = 0.5,
       fill_color = c("red", "orange", "yellow")) +
  ggtitle("PostModulator Vs PreModulator downregulated transcript overlap") +
  theme(plot.title = element_text(vjust = 0, hjust = 0.5, size = rel(1.2), face = "bold"))
ggsave("de_results_2024/transcriptomics/postmod_premod/overlap_down.png")

common <- intersect(intersect(de.CFRD.down$Molecule, de.IGT.down$Molecule), de.NGT.down$Molecule)
common

setdiff(intersect(de.CFRD.down$Molecule, de.IGT.down$Molecule), common)
setdiff(intersect(de.CFRD.down$Molecule, de.NGT.down$Molecule), common)
setdiff(intersect(de.IGT.down$Molecule, de.NGT.down$Molecule), common)
