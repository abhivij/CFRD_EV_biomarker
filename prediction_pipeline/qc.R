#obtain count of number of samples in different groups like CF, CFRD, IGT, ...
#create boxplots and dim red plots of transcripts to check distribution across different cohorts

library(tidyverse)
library(ggrepel)
library(umap)
library(ggvenn)
library(sva)
library(harmony)
library(xlsx)

base_dir <- "~/UNSW/VafaeeLab/CysticFibrosisGroup/ExoCF/CFRD_EV_biomarker/"
setwd(base_dir)

phenotype <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t")
summary(factor(phenotype$condition))


count_df <- data.frame(matrix(nrow = 0, ncol = 4, dimnames = list(c(),
                                                                  c("Class", "Total", "AU", "DK"))))
count_row <- c()
for(cl in c("HC", "IND", "CF_non_modulator", "CF_pre_post_modulator", 
            "CFRD", "IGT", "NGT")){
  #cl <- "CFRD"
  samples_subset_by_class <- phenotype %>%
    filter(is.na(pre_post_modulator) | pre_post_modulator == 0) %>%
    filter(condition == cl)
  count_row["total"] <- nrow(samples_subset_by_class)
  for(ct in c("AU", "DK")){
    #ct <- "AU"
    samples_subset_by_class_country <- samples_subset_by_class %>%
      filter(country == ct)
    count_row[ct] <- nrow(samples_subset_by_class_country)  
  }
  count_df[nrow(count_df)+1, ] <- c(cl, count_row)
  count_row <- c()
}

#CF
cl <- "CF"
samples_subset_by_class <- phenotype %>%
  filter(is.na(pre_post_modulator) | pre_post_modulator == 0) %>%
  filter(condition != "HC")
count_row["total"] <- nrow(samples_subset_by_class)
for(ct in c("AU", "DK")){
  samples_subset_by_class_country <- samples_subset_by_class %>%
    filter(country == ct)
  count_row[ct] <- nrow(samples_subset_by_class_country)    
}
count_df[nrow(count_df)+1, ] <- c(cl, count_row)
count_row <- c()



#premodulator and postmodulator in pre_post_modulator groups
cl <- "pre_modulator-prepostgroup"
samples_subset_by_class <- phenotype %>%
  filter(is.na(pre_post_modulator) | pre_post_modulator == 0) %>%
  filter(condition %in% c("CF_non_modulator", "CF_pre_post_modulator"))
count_row["total"] <- nrow(samples_subset_by_class)
for(ct in c("AU", "DK")){
  samples_subset_by_class_country <- samples_subset_by_class %>%
    filter(country == ct)
  count_row[ct] <- nrow(samples_subset_by_class_country)    
}
count_df[nrow(count_df)+1, ] <- c(cl, count_row)
count_row <- c()


cl <- "post_modulator-prepostgroup"
samples_subset_by_class <- phenotype %>%
  filter(!is.na(pre_post_modulator) & pre_post_modulator == 1) %>%
  filter(condition %in% c("CF_non_modulator", "CF_pre_post_modulator"))
count_row["total"] <- nrow(samples_subset_by_class)
for(ct in c("AU", "DK")){
  samples_subset_by_class_country <- samples_subset_by_class %>%
    filter(country == ct)
  count_row[ct] <- nrow(samples_subset_by_class_country)    
}
count_df[nrow(count_df)+1, ] <- c(cl, count_row)
count_row <- c()



#premodulator and postmodulator among all groups
cl <- "pre_modulator"
samples_subset_by_class <- phenotype %>%
  filter(is.na(pre_post_modulator) | pre_post_modulator == 0) %>%
  filter(condition != "HC")
count_row["total"] <- nrow(samples_subset_by_class)
for(ct in c("AU", "DK")){
  samples_subset_by_class_country <- samples_subset_by_class %>%
    filter(country == ct)
  count_row[ct] <- nrow(samples_subset_by_class_country)    
}
count_df[nrow(count_df)+1, ] <- c(cl, count_row)
count_row <- c()


cl <- "post_modulator"
samples_subset_by_class <- phenotype %>%
  filter(!is.na(pre_post_modulator) & pre_post_modulator == 1) %>%
  filter(condition != "HC")
count_row["total"] <- nrow(samples_subset_by_class)
for(ct in c("AU", "DK")){
  samples_subset_by_class_country <- samples_subset_by_class %>%
    filter(country == ct)
  count_row[ct] <- nrow(samples_subset_by_class_country)    
}
count_df[nrow(count_df)+1, ] <- c(cl, count_row)
count_row <- c()



write.table(count_df, "data/formatted/summary.csv", sep = ",", row.names = FALSE)



####################

#create dimensionality reduction plots (UMAP and tSNE) for AU and DK cohorts

comparison = "CFRDVsIGT"
classes = c("IGT", "CFRD")
perform_filter = TRUE
use_train_param = FALSE
#norm = "none"
norm = "norm_log_tmm"
use_best_transcripts = TRUE
best_features_file_path = "data/selected_features/best_features_with_is_best.csv"

#function to plot umap or tsne dimensionality reduction plots of the transcripts in the au and dk cohort
plot_data <- function(comparison, classes, 
                              best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                              use_best_transcripts = TRUE,
                              perform_filter = TRUE, norm = "norm_log_tmm", use_train_param = TRUE){
  
  data <- read.table("data/formatted/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                     nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
  phenotype <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t")
  
  output_labels.train <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes, country == "AU") %>%
    dplyr::select(Sample, Label, age, age_group, sex, FEV1) %>%
    dplyr::mutate(Label = factor(Label), age = as.numeric(age), 
                  age_group = factor(age_group), sex = factor(sex),
                  FEV1 = as.numeric(FEV1)) %>%
    arrange(Label, Sample)
  data.train <- data[, output_labels.train$Sample]
  print(summary(output_labels.train))
  
  output_labels.test <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes, country == "DK") %>%
    dplyr::select(Sample, Label, age, age_group, sex, FEV1) %>%
    dplyr::mutate(Label = factor(Label), age = as.numeric(age), 
                  age_group = factor(age_group), sex = factor(sex),
                  FEV1 = as.numeric(FEV1)) %>%    
    arrange(Label, Sample)
  data.test <- data[, output_labels.test$Sample]
  print(summary(output_labels.test))
  
  #currently data.train, data.test format : (transcripts x samples)
  
  if(perform_filter){
    keep <- edgeR::filterByExpr(data.train, group = output_labels.train$Label)
    data.train <- data.train[keep, ]
    if(use_train_param){
      data.test <- data.test[keep, ]  
    } else{
      keep <- edgeR::filterByExpr(data.test, group = output_labels.test$Label)  
      data.test <- data.test[keep, ] 
    }
  }
  
  
  if(norm == "norm_log_tmm"){
    #calculating norm log tmm
    dge <- edgeR::DGEList(counts = data.train, group = output_labels.train$Label)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    tmm <- edgeR::cpm(dge, log = TRUE)
    data.train <- tmm
    
    dge <- edgeR::DGEList(counts = data.test, group = output_labels.test$Label)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    tmm <- edgeR::cpm(dge, log = TRUE)
    data.test <- tmm
    
    data.train <- as.data.frame(t(as.matrix(data.train)))
    data.test <- as.data.frame(t(as.matrix(data.test)))  
    
    #normalizing the data
    normparam <- caret::preProcess(data.train) 
    data.train <- predict(normparam, data.train)
    
    if(use_train_param){
      data.test <- predict(normparam, data.test) #normalizing test data using params from train data       
    } else{
      normparam <- caret::preProcess(data.test)
      data.test <- predict(normparam, data.test) #normalizing test data using params from train data 
    }
  }
  
  #now data.train, data.test format : (samples x transcripts)
  
  #get best biomarkers only
  if(use_best_transcripts){
    best_features <- read.csv(best_features_file_path)  
    best_features_sub <- best_features %>%
      mutate(dataset_id = gsub("CF_EV_AU_zlogtmm_", "", dataset_id)) %>%
      filter(is_best == 1, dataset_id == comparison)
    
    biomarkers <- strsplit(best_features_sub$biomarkers, split = "|", fixed = TRUE)[[1]]  
    
    features_with_slash <- colnames(data.train)[grepl("/", colnames(data.train), fixed = TRUE)] 
    for(f in features_with_slash){
      f_replaced <- gsub("/|-", ".", f) 
      if(f_replaced %in% biomarkers){
        biomarkers[biomarkers == f_replaced] = f
      }
    }
    biomarkers <- gsub(".", "-", biomarkers, fixed = TRUE)
    
    
    data.train <- data.train[, biomarkers]
    data.test <- data.test[, biomarkers]    
  }
  
  if(perform_filter && norm == "norm_log_tmm"){
    pp = TRUE
  } else{
    pp = FALSE
  }
  
  # create_dim_red_plots(data = data.train, groups = output_labels.train$Label, 
  #                      title_prefix = paste(comparison, 
  #                                           "best_transcripts", use_best_transcripts,
  #                                           "pp", pp, "train_params", use_train_param,
  #                                           "train"), 
  #                      dim_red = "UMAP")
  # create_dim_red_plots(data = data.test, groups = output_labels.test$Label, 
  #                      title_prefix = paste(comparison, 
  #                                           "best_transcripts", use_best_transcripts,
  #                                           "pp", pp, "train_params", use_train_param,
  #                                           "test"), 
  #                      dim_red = "UMAP")  
  # 
  # 
  # create_dim_red_plots(data = data.train, groups = output_labels.train$Label, 
  #                      title_prefix = paste(comparison, 
  #                                           "best_transcripts", use_best_transcripts,
  #                                           "pp", pp, "train_params", use_train_param,
  #                                           "train"), 
  #                      dim_red = "tSNE")
  # create_dim_red_plots(data = data.test, groups = output_labels.test$Label, 
  #                      title_prefix = paste(comparison, 
  #                                           "best_transcripts", use_best_transcripts,
  #                                           "pp", pp, "train_params", use_train_param,
  #                                           "test"), 
  #                      dim_red = "tSNE")  
  # 
  
  create_box_plot(data = data.train, groups = output_labels.train$Label, 
                  title_prefix = paste(comparison, 
                                       "best_transcripts", use_best_transcripts,
                                       "pp", pp, "train_params", use_train_param,
                                       "train"))
  create_box_plot(data = data.test, groups = output_labels.test$Label, 
                  title_prefix = paste(comparison, 
                                       "best_transcripts", use_best_transcripts,
                                       "pp", pp, "train_params", use_train_param,
                                       "test"))
  
}


# does not work
# plot_data(comparison = "CFRDVsIGT", classes = c("IGT", "CFRD"), 
#                   best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
#                   use_best_transcripts = FALSE,
#                   perform_filter = FALSE, norm = "none", use_train_param = FALSE)

plot_data(comparison = "CFRDVsIGT", classes = c("IGT", "CFRD"), 
                  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                  use_best_transcripts = FALSE,
                  perform_filter = TRUE, norm = "norm_log_tmm", use_train_param = FALSE)

plot_data(comparison = "CFRDVsIGT", classes = c("IGT", "CFRD"), 
                  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                  use_best_transcripts = FALSE,
                  perform_filter = TRUE, norm = "norm_log_tmm", use_train_param = TRUE)

#does not work
# plot_data(comparison = "CFRDVsIGT", classes = c("IGT", "CFRD"), 
#                   best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
#                   use_best_transcripts = TRUE,
#                   perform_filter = FALSE, norm = "none", use_train_param = FALSE)

plot_data(comparison = "CFRDVsIGT", classes = c("IGT", "CFRD"), 
                  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                  use_best_transcripts = TRUE,
                  perform_filter = TRUE, norm = "norm_log_tmm", use_train_param = FALSE)

plot_data(comparison = "CFRDVsIGT", classes = c("IGT", "CFRD"), 
                  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                  use_best_transcripts = TRUE,
                  perform_filter = TRUE, norm = "norm_log_tmm", use_train_param = TRUE)

# data <- data.train
# groups <- output_labels.train$Label
# perplexity = 5
# shownames = TRUE
# dim_red = "tSNE"
# colour_label = "Condition"
# dim_red = "UMAP"
# data_cohort = "train (AU)"
# use_best_transcripts = FALSE
# 
# title_prefix = paste(comparison,
#                      "best_transcripts", FALSE,
#                      "pp", TRUE, "train_params", TRUE,
#                     "train")
# dim_red = "UMAP"


#THE FUNCTION BELOW WAS CREATED PREVIOUSLY
#NO LONGER USED
#data and groups expected to be matching
#should be ensured from the caller function
# create_dim_red_plots <- function(data, groups,
#                                  title_prefix = "",
#                                  dim_red = "UMAP",
#                                  perplexity = 5,
#                                  colour_label = "Condition",
#                                  shownames = FALSE){
#   
#   groups_modified <- groups
#   summa <- summary(factor(groups))
#   
#   groups_modified <- gsub(pattern = names(summa)[1], 
#                           replacement = paste0(names(summa)[1], "(", summa[1], ")"),
#                           groups_modified)
#   groups_modified <- gsub(pattern = names(summa)[2], 
#                           replacement = paste0(names(summa)[2], "(", summa[2], ")"),
#                           groups_modified)
#   groups <- factor(groups_modified)
#   
#   if(shownames){
#     text <- rownames(data)
#   } else{
#     text <- ""
#   }
#   set.seed(1)
#   if(dim_red == "tSNE"){
#     result <- Rtsne::Rtsne(data, perplexity = perplexity)
#     dim_red_df <- data.frame(x = result$Y[,1], y = result$Y[,2], 
#                              Colour = groups, 
#                              Sample = text)    
#     xlab <- "tSNE 1"
#     ylab <- "tSNE 2"
#   } else if(dim_red == "UMAP"){
#     print(dim(data)[1])
#     # n_neighbors <- max(floor(dim(data)[1] / 4), 2)
#     # print(n_neighbors)
#     # result <- umap(data, n_neighbors = n_neighbors)
#     result <- umap(data)
#     dim_red_df <- data.frame(x = result$layout[,1], y = result$layout[,2], 
#                              Colour = groups, 
#                              Sample = text)  
#     xlab <- "UMAP 1"
#     ylab <- "UMAP 2"
#   }
#   
#   title <- paste(title_prefix, dim_red)
#   
#   if (shownames) {
#     dim_red_plot <- ggplot2::ggplot(dim_red_df,
#                                     ggplot2::aes(x = x, y = y, colour = Colour)) +
#       ggplot2::geom_point() +
#       geom_text_repel(aes(label=Sample)) +
#       ggplot2::labs(title = title, colour = colour_label) +
#       ggplot2::xlab(xlab) +
#       ggplot2::ylab(ylab) +
#       labs(caption = paste("Data dimension :", paste(dim(data), collapse = "x")))
#     dim_red_plot
#   } else {
#     dim_red_plot <- ggplot2::ggplot(dim_red_df) +
#       ggplot2::geom_point(ggplot2::aes(x = x, y = y, colour = Colour)) +
#       ggplot2::labs(title = title, colour = colour_label) +
#       ggplot2::xlab(xlab) +
#       ggplot2::ylab(ylab) +
#       labs(caption = paste("Data dimension :", paste(dim(data), collapse = "x")))
#     dim_red_plot
#   }
#   
#   dir_path <- "prediction_pipeline/plots"
#   file_name <- paste0(gsub(title, pattern = " ", replacement = "-"), ".jpg")
#   file_path <- paste(dir_path, file_name, sep = "/")
#   ggplot2::ggsave(file_path, dim_red_plot, units = "cm", width = 30)
# }


#data and groups expected to be matching
#should be ensured from the caller function
create_box_plot <- function(data, groups,
                            title_prefix = ""){
  data_to_plot <- cbind(data, "label" = groups) %>%
    rownames_to_column("sample_name")
  data_to_plot <- data_to_plot %>%
    pivot_longer(!c(sample_name, label), names_to = "transcripts") %>%
    arrange(label)
  data_to_plot <- data_to_plot %>%
    mutate(sample_name = factor(sample_name, levels = unique(data_to_plot$sample_name)))
  
  ggplot(data_to_plot, aes(x = sample_name, y = value)) +
    geom_boxplot(aes(fill = label)) +
    xlab("Sample Name") +
    ylab("Expression across transcripts") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle(title_prefix)
  
  title <- paste(title_prefix, "boxplot")
  
  dir_path <- "prediction_pipeline/plots"
  file_name <- paste0(gsub(title, pattern = " ", replacement = "-"), ".jpg")
  file_path <- paste(dir_path, file_name, sep = "/")
  ggplot2::ggsave(file_path, units = "cm", width = 20)
  
}





# does not work
# plot_data(comparison = "CFRDVsNGT", classes = c("NGT", "CFRD"), 
#                   best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
#                   use_best_transcripts = FALSE,
#                   perform_filter = FALSE, norm = "none", use_train_param = FALSE)

plot_data(comparison = "CFRDVsNGT", classes = c("NGT", "CFRD"), 
                  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                  use_best_transcripts = FALSE,
                  perform_filter = TRUE, norm = "norm_log_tmm", use_train_param = FALSE)

plot_data(comparison = "CFRDVsNGT", classes = c("NGT", "CFRD"), 
                  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                  use_best_transcripts = FALSE,
                  perform_filter = TRUE, norm = "norm_log_tmm", use_train_param = TRUE)

#does not work
# plot_data(comparison = "CFRDVsNGT", classes = c("NGT", "CFRD"), 
#                   best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
#                   use_best_transcripts = TRUE,
#                   perform_filter = FALSE, norm = "none", use_train_param = FALSE)

plot_data(comparison = "CFRDVsNGT", classes = c("NGT", "CFRD"), 
                  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                  use_best_transcripts = TRUE,
                  perform_filter = TRUE, norm = "norm_log_tmm", use_train_param = FALSE)

plot_data(comparison = "CFRDVsNGT", classes = c("NGT", "CFRD"), 
                  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                  use_best_transcripts = TRUE,
                  perform_filter = TRUE, norm = "norm_log_tmm", use_train_param = TRUE)




# does not work
# plot_data(comparison = "IGTVsNGT", classes = c("NGT", "IGT"), 
#                   best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
#                   use_best_transcripts = FALSE,
#                   perform_filter = FALSE, norm = "none", use_train_param = FALSE)

plot_data(comparison = "IGTVsNGT", classes = c("NGT", "IGT"), 
                  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                  use_best_transcripts = FALSE,
                  perform_filter = TRUE, norm = "norm_log_tmm", use_train_param = FALSE)

plot_data(comparison = "IGTVsNGT", classes = c("NGT", "IGT"), 
                  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                  use_best_transcripts = FALSE,
                  perform_filter = TRUE, norm = "norm_log_tmm", use_train_param = TRUE)

#does not work
# plot_data(comparison = "IGTVsNGT", classes = c("NGT", "IGT"), 
#                   best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
#                   use_best_transcripts = TRUE,
#                   perform_filter = FALSE, norm = "none", use_train_param = FALSE)

plot_data(comparison = "IGTVsNGT", classes = c("NGT", "IGT"), 
                  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                  use_best_transcripts = TRUE,
                  perform_filter = TRUE, norm = "norm_log_tmm", use_train_param = FALSE)

plot_data(comparison = "IGTVsNGT", classes = c("NGT", "IGT"), 
                  best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                  use_best_transcripts = TRUE,
                  perform_filter = TRUE, norm = "norm_log_tmm", use_train_param = TRUE)



###############


#plots for FEV1, Age, Agegroup, Sex, PatientRecruitmentYear

data_of_interest <- phenotype %>%
  filter(is.na(pre_post_modulator) | pre_post_modulator == 0) %>%
  filter(condition %in% c("CFRD", "IGT", "NGT")) %>%
  select(Sample, condition, country, patient_recruitment_year, age, age_group, sex, FEV1) %>%
  mutate(age = as.numeric(str_trim(age)), FEV1 = as.numeric(str_trim(FEV1)))

#FEV1 and Sex missing for child samples


#plot FEV1

data_to_plot <- data_of_interest %>%
  select(Sample, condition, country, FEV1)
ggplot(data_to_plot, aes(x = condition, y = FEV1)) +
  geom_boxplot(aes(fill = country)) +
  xlab("Condition") +
  ylab("FEV1") +
  #specifying 10 here causes nearest good value to 10 to be chosen for number of breaks in y-axis
  scale_y_continuous(n.breaks = 10)
ggsave("prediction_pipeline/plots/other_features/FEV1_original.jpg")


#FEV1 replace missing with average
mean_values <- data_of_interest %>% filter(country == "AU") %>%
  group_by(condition) %>%
  summarise(mean_FEV1 = mean(FEV1, na.rm = TRUE))

#manually verified that DK samples have no missing
#so replacing NAs does not affect DK samples
data_of_interest <- data_of_interest %>%
  inner_join(mean_values, by = "condition")
data_of_interest <- data_of_interest %>%
  mutate(FEV1 = case_when(is.na(FEV1) ~ mean_FEV1,
                          TRUE ~ FEV1))

data_to_plot <- data_of_interest %>%
  select(Sample, condition, country, FEV1)
ggplot(data_to_plot, aes(x = condition, y = FEV1)) +
  geom_boxplot(aes(fill = country)) +
  xlab("Condition") +
  ylab("FEV1") +
  #specifying 10 here causes nearest good value to 10 to be chosen for number of breaks in y-axis
  scale_y_continuous(n.breaks = 10)
ggsave("prediction_pipeline/plots/other_features/FEV1_afterreplacebymean.jpg")


data_of_interest <- data_of_interest %>%
  select(-c(mean_FEV1))
data_to_plot <- data_of_interest %>%
  select(Sample, condition, country, age)
ggplot(data_to_plot, aes(x = condition, y = age)) +
  geom_boxplot(aes(fill = country)) +
  xlab("Condition") +
  ylab("Age") +
  #specifying 10 here causes nearest good value to 10 to be chosen for number of breaks in y-axis
  scale_y_continuous(n.breaks = 10)
ggsave("prediction_pipeline/plots/other_features/age.jpg")


#NA values in Sex column not replaced
data_to_plot <- data_of_interest %>%
  select(Sample, condition, country, sex)
ggplot(data_to_plot, aes(x = condition)) +
  geom_bar(position = "dodge", aes(fill = country)) +
  xlab("Condition") +
  facet_wrap(~sex) +
  scale_y_continuous(n.breaks = 10)
ggsave("prediction_pipeline/plots/other_features/sex.jpg")


data_to_plot <- data_of_interest %>%
  select(Sample, condition, country, patient_recruitment_year)
ggplot(data_to_plot, aes(x = condition)) +
  geom_bar(position = "dodge", aes(fill = country)) +
  xlab("Condition") +
  facet_wrap(~patient_recruitment_year) +
  scale_y_continuous(n.breaks = 10)
ggsave("prediction_pipeline/plots/other_features/patient_recruitment_year.jpg")


data_to_plot <- data_of_interest %>%
  select(Sample, condition, country, age_group)
ggplot(data_to_plot, aes(x = condition)) +
  geom_bar(position = "dodge", aes(fill = country)) +
  xlab("Condition") +
  facet_wrap(~age_group) +
  scale_y_continuous(n.breaks = 5)
ggsave("prediction_pipeline/plots/other_features/age_group.jpg")


###############
#compare best biomarker set from AU and DK

best_features_file_path = "data/selected_features/best_features_with_is_best.csv"
comparison = "CFRDVsIGT"

plot_best_biomarker_venn <- function(comparison,
                                     plot_dir_path = "plots/best_biomarkers",
                                     best_features_file_path = "data/selected_features/best_features_with_is_best.csv"){
  best_features <- read.csv(best_features_file_path)  
  
  #from AU
  best_features_sub <- best_features %>%
    filter(is_best == 1, dataset_id == paste0("CF_EV_AU_adult_logtmm_", comparison))
  biomarkers_AU <- strsplit(best_features_sub$biomarkers, split = "|", fixed = TRUE)[[1]]  
  
  #from DK
  best_features_sub <- best_features %>%
    filter(is_best == 1, dataset_id == paste0("CF_EV_DK_adult_logtmm_", comparison))
  biomarkers_DK <- strsplit(best_features_sub$biomarkers, split = "|", fixed = TRUE)[[1]]  
  
  ggvenn(list("AU" = biomarkers_AU,
              "DK" = biomarkers_DK),
         stroke_size = 0.1,
         set_name_size = 4,
         text_size = 3,
         fill_color = c("#F8766D", "#00BFC4")) +
    ggtitle(paste("Best biomarkers in ", sub("Vs", " Vs ", comparison))) +
    theme(plot.title = element_text(vjust = - 20, hjust = 0.5))
  file_name <- paste0("plots/best_biomarkers/", comparison, ".png")
  ggsave(file_name)
}

plot_best_biomarker_venn("CFRDVsIGT")
plot_best_biomarker_venn("CFRDVsNGT")
plot_best_biomarker_venn("IGTVsNGT")


###############
#create box plots of best biomarker sets and all transcripts
train_cohort_country = "AU"

data_file_path = "data/formatted/umi_counts_filtered_seurat3_with_norm_and_find_var_feat.csv"
comparison = "CFRDVsIGT"
classes  = c("IGT", "CFRD")
use_best_features = TRUE
norm = "none"
best_features_file_path = "data/selected_features/best_features_with_is_best.csv"
perform_filter = FALSE
combined_plot = TRUE
dataset_replace_str = "CF_EV_adult_filtered_seurat3_norm_find_var_none_"
plot_width_cm = 25
dir_path = "prediction_pipeline/plots/box_plots/filtered_then_seurat3_norm_and_find_var_feat_custom"
plot_title_input = "filter and then seurat3 best biomarkers CFRD Vs IGT"

data_file_path = "data/formatted/umi_counts_filtered_seurat3_with_norm_and_find_var_feat.csv"
comparison = "CFRDVsNGT"
classes  = c("NGT", "CFRD")
use_best_features = TRUE
norm = "none"
best_features_file_path = "data/selected_features/best_features_with_is_best.csv"
perform_filter = FALSE
combined_plot = TRUE
dataset_replace_str = "CF_EV_adult_filtered_seurat3_norm_find_var_none_"
plot_width_cm = 25
dir_path = "prediction_pipeline/plots/box_plots/filtered_then_seurat3_norm_and_find_var_feat_custom"
plot_title_input = "filter and then seurat3 best biomarkers CFRD Vs NGT"
y_lim = NA

#filter_type - possible values : combined - apply filter taking AU+DK together
#                                regular  - filter on train, apply on test
create_transcript_box_plots <- function(comparison,
                                        classes,
                                        use_best_features = TRUE,
                                        train_cohort_country = "AU",
                                        norm = "none",
                                        best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                                        perform_filter = TRUE,
                                        combat_seq = FALSE,
                                        filter_type = "regular",
                                        dataset_replace_str = "",
                                        plot_width_cm = 20,
                                        dir_path = "prediction_pipeline/plots/box_plots",
                                        data_file_path = NA,
                                        combined_plot = FALSE,
                                        plot_title_input = NA,
                                        y_lim = c(),
                                        biomarkers_specific = c()){
  
  if(is.na(data_file_path)){
    if(combat_seq){
      data <- read.table("data/formatted/umi_counts_combat_seq.csv", header=TRUE, sep=",", row.names=1, skip=0,
                         nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
    } else{
      data <- read.table("data/formatted/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                         nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
    }  
  } else{
    data <- read.table(data_file_path, header=TRUE, sep=",", row.names=1, skip=0,
                       nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
  }
  
  phenotype <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t")
  phenotype <- phenotype %>% filter(age_group == "adult")
  
  test_cohort_country <- ifelse(train_cohort_country == "AU", "DK", "AU") 
  
  train.output_labels <- phenotype %>%
    dplyr::rename("Label" = comparison) %>%
    filter(Label %in% classes, country == train_cohort_country) %>%
    dplyr::select(Sample, Label) %>%
    dplyr::mutate(Label = factor(Label)) %>%
    arrange(Label, Sample)
  train.tra_data <- data[, train.output_labels$Sample]
  
  test.output_labels <- phenotype %>%
    dplyr::rename("Label" = comparison) %>%
    filter(Label %in% classes, country == test_cohort_country) %>%
    dplyr::select(Sample, Label) %>%
    dplyr::mutate(Label = factor(Label)) %>%    
    arrange(Label, Sample)
  test.tra_data <- data[, test.output_labels$Sample]

  #currently data format : (transcripts x samples)
  if(perform_filter){
    if(filter_type == "regular"){
      keep <- edgeR::filterByExpr(train.tra_data, group = train.output_labels$Label)
      train.tra_data <- train.tra_data[keep, ]
      test.tra_data <- test.tra_data[keep, ]  
      
    } else if(filter_type == "combined"){
      output_labels <- phenotype %>%
        rename("Label" = comparison) %>%
        filter(Label %in% classes) %>%
        dplyr::mutate(Label = factor(Label)) %>%
        arrange(Label, Sample)
      data_of_interest <- data[, output_labels$Sample]
      
      keep <- edgeR::filterByExpr(data_of_interest, group = output_labels$Label)
      data_of_interest <- data_of_interest[keep, ]

      train.tra_data <- data_of_interest[, train.output_labels$Sample]
      test.tra_data <- data_of_interest[, test.output_labels$Sample]
    }
  } else{
    filter_type <- "none"
  }

  if(norm == "log_tmm"){
    dge <- edgeR::DGEList(counts = train.tra_data, group = train.output_labels$Label)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    tmm <- edgeR::cpm(dge, log = TRUE)
    train.tra_data <- tmm
    
    dge <- edgeR::DGEList(counts = test.tra_data, group = test.output_labels$Label)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    tmm <- edgeR::cpm(dge, log = TRUE)
    test.tra_data <- tmm
  } else if(norm == "log_cpm"){
    train.tra_data <- edgeR::cpm(train.tra_data, log = TRUE)
    test.tra_data <- edgeR::cpm(test.tra_data, log = TRUE)
  }
  train.tra_data <- as.data.frame(t(as.matrix(train.tra_data)))
  test.tra_data <- as.data.frame(t(as.matrix(test.tra_data)))  
  #now data : (samples x transcripts)

  if(use_best_features){
    best_features <- read.csv(best_features_file_path)  
    
    if(dataset_replace_str == ""){
      dataset_replace_str <- paste0("CF_EV_", train_cohort_country, 
                                    "_adult_logtmm_")      
    }
    
    best_features_sub <- best_features %>%
      mutate(dataset_id = gsub(dataset_replace_str, "", dataset_id)) %>%
      filter(is_best == 1, dataset_id == comparison)
    
    biomarkers <- strsplit(best_features_sub$biomarkers, split = "|", fixed = TRUE)[[1]]     
    
    features_with_slash <- colnames(train.tra_data)[grepl("/", colnames(train.tra_data), fixed = TRUE)] 
    for(f in features_with_slash){
      f_replaced <- gsub("/|-", ".", f) 
      if(f_replaced %in% biomarkers){
        biomarkers[biomarkers == f_replaced] = f
      }
    }
    biomarkers <- gsub(".", "-", biomarkers, fixed = TRUE)
    
    print(biomarkers)
    print(length(biomarkers))
    biomarkers_present <- biomarkers[biomarkers %in% colnames(train.tra_data)]
    
    print(biomarkers_present)
    print(length(biomarkers_present))
    
    train.tra_data <- train.tra_data[, biomarkers_present]
    test.tra_data <- test.tra_data[, biomarkers_present]
  }
  num_transcripts <- dim(train.tra_data)[2]
  
  data_to_plot <- rbind(
    cbind(train.tra_data, 
          "label" = train.output_labels$Label,
          "country" = train_cohort_country),
    cbind(test.tra_data,
          "label" = test.output_labels$Label,
          "country" = test_cohort_country)
  )
  if(!combined_plot){
    data_to_plot <- data_to_plot %>%
      mutate(label = paste(country, label, sep = "_"))
  }
  data_to_plot <- data_to_plot %>%
    select(-c(country)) %>%
    rownames_to_column("sample_name")  
  
  group_counts <- data_to_plot %>%
    select(c(sample_name, label)) %>%
    group_by(label) %>%
    summarise(n = n())
  
  data_to_plot <- data_to_plot %>%
    inner_join(group_counts) %>%
    mutate(label = paste0(label, " (", n, ") ")) %>%
    select(-c(n))
  
  data_to_plot <- data_to_plot %>%
    pivot_longer(!c(sample_name, label), names_to = "transcripts") %>%
    arrange(transcripts, label)
  
  if(norm == "log_tmm"){
    y_lab <- "Log TMM of expression level across samples"
    plot_title <- "Log TMM expression of "
  } else if(norm == "log_cpm"){
    y_lab <- "Log CPM of expression level across samples"
    plot_title <- "Log CPM expression of "    
  } else{
    y_lab <- "Expression level across samples"
    plot_title <- "Expression of "
  }
  
  if(use_best_features){
    plot_title <- paste0(plot_title, "best transcripts from ", train_cohort_country)
  } else{
    plot_title <- paste0(plot_title, "all transcripts with params from ", train_cohort_country)
  }
  plot_title <- paste0(plot_title, " in classes ", paste0(rev(classes), collapse = ","))
  plot_title <- paste0(plot_title, " filter_type ", filter_type)

  if(!is.na(plot_title_input)){
    plot_title = plot_title_input
  }
  
  # biomarkers_specific <- c("hsa-miR-122-5p", "hsa-miR-342-3p", "hsa-miR-486-3p", "hsa-miR-182-5p",
  #                          "hsa-piR-020497", "hsa-miR-17-5p", "hsa-piR-016926")
  data_to_plot_sub <- data_to_plot %>%
    filter(!transcripts %in% biomarkers_specific) %>%
    group_by(transcripts) %>%
    summarize(median = median(value), mean = mean(value)) %>%
    arrange(desc(median), desc(mean))
  biomarkers_rem <- data_to_plot_sub$transcripts
  
  biomarkers_ordered <- c(biomarkers_specific, biomarkers_rem)

  data_to_write <- data.frame(transcripts = biomarkers_ordered)
  write.xlsx(data_to_write,
             "data/formatted/identified_biomarkers.xlsx",
             sheetName = comparison,
             col.names = TRUE, row.names = FALSE, append = TRUE)
  
  data_to_plot <- data_to_plot %>%
    mutate(transcripts = factor(transcripts, levels = biomarkers_ordered)) %>%
    mutate(plot_num = ifelse(transcripts %in% biomarkers_specific,
                               1, 2))
  plot <- ggplot() +
    aes(x = transcripts, y = value, fill = label) +
    geom_boxplot(data = (data_to_plot %>% filter(plot_num == 1))) +
    geom_boxplot(data = (data_to_plot %>% filter(plot_num == 2))) +
    xlab(paste0("Transcripts (", num_transcripts, ")")) +
    ylab(y_lab) +
    ggtitle(plot_title) +
    labs(fill = "") +
    facet_wrap(~plot_num, scales = "free") + 
    theme(
      strip.text.x = element_blank()
    )

  if(length(y_lim) != 0){
    plot <- plot + ylim(y_lim)
  }
  if(use_best_features){
    plot <- plot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  } else {
    plot <- plot + theme(axis.text.x = element_blank())
  }
  plot
  
  if(!dir.exists(dir_path)){
    dir.create(dir_path, recursive = TRUE)
  }
  file_name <- paste0(gsub(plot_title, pattern = " ", replacement = "-"), ".png")
  file_path <- paste(dir_path, file_name, sep = "/")
  ggplot2::ggsave(file_path, units = "cm", width = plot_width_cm)
}


#############
#create dim_red plots with both AU and DK cohort in one figure
#this is to check if there are separate clusters for AU and DK within each of the pairs of CFRD, IGT, NGT
# comparison = "CFRDVsIGT" #this will be NA if 3 classes are involved
#                          #in that case non-modulator sample filtering will be separately done too
# classes = c("NGT", "CFRD")
# 
# comparison <- NA
# classes = c("CFRD", "IGT", "NGT")
# dim_red = "UMAP"
# shownames = FALSE
# norm = "log_tmm"

fill_column = NA
combat_seq = FALSE
combat_seq_specify_group = FALSE
combatref = FALSE
shownames = FALSE
perform_filter = TRUE
best_features_file_path = NA
dataset_replace_str = NA
dir_path = NA
data_file_path = "data/formatted/umi_counts.csv"
filter_out_child = FALSE
plot_title_prefix = ""
plot_title_suffix = ""
box_plot_dir_path <- "plots_updated/boxplots"

comparison = "CFRDVsIGT"
classes = c("CFRD", "IGT")
class_colours = c("red", "orange")
dim_red = "UMAP"
norm = "log_tmm"
combat = TRUE
dir_path = "plots_updated/dim_red"
plot_width_cm = 21

best_features_file_path  = "data/selected_features/best_features_with_is_best.csv"
dataset_replace_str = "CF_EV_tra_combat_"
colour_column = "age_group"


combattwice = FALSE

comparison = "CFRDVsIGT"
classes = c("CFRD", "IGT")
class_colours = c("red", "orange")
dim_red = "UMAP"
norm = "non-normalized"
dir_path = "plots_updated/tra_334"
plot_width = 21
perform_filter = TRUE
colour_column = "batch_name"
data_file_path = "data/formatted/rna_all/umi_counts_filter90.csv"
phenotype_file_path = "data/formatted/rna_all/phenotype.txt"
plot_title_prefix = "1 "


comparison = NA
classes = c("PostModulator_CFRD", "PostModulator_IGT", "PostModulator_NGT",
            "PreModulator_CFRD", "PreModulator_IGT", "PreModulator_NGT")
class_colours = c("indianred", "coral", "khaki", 
                  "red", "orange", "yellow")
dim_red = "UMAP"
norm = "non-normalized"
dir_path = "plots_updated/post_mod/proteomics"
plot_width = 21
perform_filter = FALSE
colour_column = "batch_name"
point_border_colours = c("black", "green", "purple")
data = data
phenotype = phenotype
filter_post_modulator = FALSE
custom_title = "proteomics pre and post modulator samples with Quantile norm + ComBat"

create_dim_red_plots <- function(comparison, classes,
                                 class_colours,
                                 dim_red, norm,
                                 fill_column = NA,
                                 simplified = FALSE,
                                 colour_column = "age_group",
                                 point_border_colours = c("black", "green"),
                                 combat_seq = FALSE, combat_seq_specify_group = FALSE,
                                 combatref = FALSE,
                                 combat = FALSE,
                                 combattwice = FALSE,
                                 shownames = FALSE,
                                 perform_filter = TRUE,
                                 best_features_file_path = NA, 
                                 dataset_replace_str = NA,
                                 dir_path = NA,
                                 data_file_path = "data/formatted/umi_counts.csv",
                                 phenotype_file_path = "data/formatted/phenotype.txt",
                                 data = NA,
                                 phenotype = NA,
                                 filter_out_child = FALSE,
                                 dimred_plot_width_cm = 30,
                                 plot_width_cm = 21,
                                 plot_title_prefix = "",
                                 plot_title_suffix = "",
                                 omics_type = "tra",
                                 box_plot_dir_path = "plots_updated/boxplots",
                                 biomarker_expr_file_path = "data/selected_features/biomarkers_expr.xlsx",
                                 filter_post_modulator = TRUE,
                                 custom_title = NA){
  if(is.null(dim(data))){
    data <- read.table(data_file_path, header=TRUE, sep=",", row.names=1, skip=0,
                       nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")    
  }
  if(is.null(dim(phenotype))){
    phenotype <- read.table(phenotype_file_path, header=TRUE, sep="\t")    
  }
  
  title <- paste0(dim_red, " plot of ", paste(classes, collapse = ", "), " samples from ", norm, " data")
  
  if(!is.na(comparison)){
    output_labels <- phenotype %>%
      dplyr::rename("Label" = comparison) %>%
      dplyr::rename("colour_column" = colour_column)
  } else{
    #this case is to include greater than 2 conditions in the plot
    output_labels <- phenotype %>%
      dplyr::rename("Label" = "condition") %>%
      dplyr::rename("colour_column" = colour_column)
    
    #ensure that samples on modulator are not chosen
    #using comparison field does not require the below, because comparison field creation 
    #                 incorporates this filter
    if(filter_post_modulator){
      output_labels <- output_labels %>%
        filter(is.na(pre_post_modulator) | pre_post_modulator == 0)      
    }  

  }
  output_labels <- output_labels %>%
    filter(Label %in% classes) %>%
    dplyr::select(Sample, Label, country, colour_column, mutation, 
                  sex, cohort, patient_recruitment_year,
                  seq_plate, seq_miR_library_quality) %>%
    dplyr::mutate(Label = factor(Label)) %>%
    dplyr::mutate(seq_plate = factor(seq_plate),
                  patient_recruitment_year = factor(patient_recruitment_year)) %>%
    arrange(Label, Sample)
  
  if(filter_out_child){
    output_labels <- output_labels %>%
      filter(age_group != "child")
  }
  
  if(!is.na(fill_column)){
    output_labels <- output_labels %>%
      rename("fill_column" = fill_column)
    title <- paste0(title, " based on ", fill_column)
  }

  group_counts <- output_labels %>%
    dplyr::mutate(Label = paste(country, Label, sep = "_")) %>%
    group_by(Label) %>%
    summarise(n = n())
  
  group_counts_text <- paste(apply(group_counts, MARGIN = 1, FUN = function(x){paste(x[1], x[2], sep = ":")}),
                             collapse = "  ")
  data <- data[, output_labels$Sample]
  
  #currently data format : (transcripts x samples)
  if(perform_filter){
    keep <- edgeR::filterByExpr(data, group = output_labels$Label)
    data <- data[keep, ]
  }
  
  if(combat_seq){
    if(combat_seq_specify_group){
      data.combat_seq <- ComBat_seq(counts = as.matrix(data),
                                    batch = factor(output_labels$country),
                                    group = factor(output_labels$Label))  
      title <- paste0(title, " after combat seq with group")
    } else{
      data.combat_seq <- ComBat_seq(counts = as.matrix(data),
                                    batch = factor(output_labels$country))
      title <- paste0(title, " after combat seq")
    }
    # dim(data)
    # dim(data.combat_seq)
    # all.equal(rownames(data), rownames(data.combat_seq))
    # all.equal(colnames(data), colnames(data.combat_seq))
    data <- data.combat_seq
  }
  
  if(norm == "log_tmm"){
    dge <- edgeR::DGEList(counts = data, group = output_labels$Label)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    tmm <- edgeR::cpm(dge, log = TRUE)
    data <- tmm
  }else if(norm == "log"){
    #taking log of 0, causes UMAP/PCA computation to fail
    #so replace 0 with aribitrary small number
    #min value in this data other than 0 is 10
    data[data == 0] <- 2^-30
    data <- log2(data)
  } else if(norm == "log_cpm"){
    data <- edgeR::cpm(data, log = TRUE)
  } else if(norm == "quantile"){
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
  }
  data <- as.data.frame(t(as.matrix(data)))
  
  
  # if(harmony){
  #   harmony_embeddings <- HarmonyMatrix(
  #     data_mat  = data,
  #     meta_data = output_labels,
  #     vars_use  = "country"
  #   )
  # }
  
  #perform batch correction on test data using train data as reference
  if(combatref){
    ref_batch = 'AU'
    data_of_interest <- as.data.frame(t(as.matrix(data)))

    all.equal(colnames(data_of_interest), output_labels$Sample)
    
    data_of_interest.combat = ComBat(dat=data_of_interest, 
                                     batch=output_labels$country, 
                                     ref.batch = ref_batch)
    data_of_interest.combat <- as.data.frame(t(as.matrix(data_of_interest.combat)))
    
    output_labels.au <- output_labels %>% filter(country == 'AU')
    data.au <- data[output_labels.au$Sample, ]
    output_labels.dk <- output_labels %>% filter(country == 'DK')
    data.dk <- data[output_labels.dk$Sample, ]
    
    data_of_interest.combat.au <- data_of_interest.combat[output_labels.au$Sample, ]
    data_of_interest.combat.dk <- data_of_interest.combat[output_labels.dk$Sample, ]
    
    #all.equal(data.au, data_of_interest.combat.au)
    #TRUE
    #all.equal(data.dk, data_of_interest.combat.dk)
    #FALSE
    data <- data_of_interest.combat
    
    title <- paste(title, '+ combat with ref as', ref_batch)
  } else if(combat){
    data_of_interest <- as.data.frame(t(as.matrix(data)))
    
    all.equal(colnames(data_of_interest), output_labels$Sample)
    
    data_of_interest.combat = ComBat(dat=data_of_interest, 
                                     batch=output_labels$country)
    if(combattwice){
      data_of_interest.combat = ComBat(dat=data_of_interest.combat, 
                                       batch=output_labels$colour_column)
    }
    data_of_interest.combat <- as.data.frame(t(as.matrix(data_of_interest.combat)))
    data <- data_of_interest.combat
    
    title <- paste(title, '+ combat')
  }
  
  if(!is.na(best_features_file_path) & !is.na(dataset_replace_str)){
    print('selecting best biomarkers')
    best_features <- read.csv(best_features_file_path)  
    best_features_sub <- best_features %>%
      mutate(dataset_id = gsub(dataset_replace_str, "", dataset_id)) %>%
      filter(is_best == 1, dataset_id == comparison)
    
    biomarkers <- strsplit(best_features_sub$biomarkers, split = "|", fixed = TRUE)[[1]]  
    
    data <- as.data.frame(t(as.matrix(data)))
    features_with_slash <- colnames(data)[grepl("/", colnames(data), fixed = TRUE)] 
    for(f in features_with_slash){
      f_replaced <- gsub("/|-", ".", f) 
      if(f_replaced %in% biomarkers){
        biomarkers[biomarkers == f_replaced] = f
      }
    }
    biomarkers <- gsub(".", "-", biomarkers, fixed = TRUE)
    data <- data[biomarkers, ]
    data <- as.data.frame(t(as.matrix(data)))
    
    title <- paste0(title, " best biomarkers")
  }
  
  
  if(shownames){
    text <- rownames(data)
  } else{
    text <- ""
  }
  set.seed(1)
  if(dim_red == "PCA"){
    result <- prcomp(data)
    dim_red_df <- data.frame(x = result$x[,1], y = result$x[,2])    
    xlab <- "PCA 1"
    ylab <- "PCA 2"  
  } else if(dim_red == "UMAP"){
    result <- umap(data)
    dim_red_df <- data.frame(x = result$layout[,1], y = result$layout[,2])  
    xlab <- "UMAP 1"
    ylab <- "UMAP 2"
  }
  title <- paste0(plot_title_prefix, title, plot_title_suffix)
  if(!is.na(custom_title)){
    title <- custom_title
  }
  
  if(simplified){
    # ggplot2::ggplot(dim_red_df, ggplot2::aes(x = x, y = y)) +
    #   ggplot2::geom_point(ggplot2::aes(colour = output_labels$colour_column), size = 3) +
    #   geom_text_repel(aes(label = text)) +
    #   ggplot2::scale_fill_manual(name = "Condition", values = class_colours) +
    #   ggplot2::guides(fill = guide_legend(override.aes = list(shape = 21,
    #                                                           colour = class_colours))) +
    #   ggplot2::scale_colour_manual(name = sub("_", " ", colour_column), values = point_border_colours) +
    #   ggplot2::guides(colour = guide_legend(override.aes = list(shape = 1))) +
    #   ggplot2::labs(title = title) +
    #   ggplot2::xlab(xlab) +
    #   ggplot2::ylab(ylab) +
    #   labs(caption = paste("Data dimension :", paste(dim(data), collapse = "x")))  
    
    ggplot2::ggplot(dim_red_df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(ggplot2::aes(fill = output_labels$Label), size = 3, shape = 21) +
      geom_text_repel(aes(label = text)) +
      ggplot2::scale_fill_manual(name = "Condition", values = class_colours) +
      ggplot2::guides(fill = guide_legend(override.aes = list(shape = 21,
                                                              colour = class_colours))) +
      ggplot2::labs(title = title) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab) +
      labs(caption = paste(paste("Data dimension :", paste(dim(data), collapse = "x")), "\n",
                           group_counts_text),
           fill = fill_column)  
  }
  else if(is.na(fill_column)){
    ggplot2::ggplot(dim_red_df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(ggplot2::aes(fill = output_labels$Label,
                                       shape = output_labels$country,
                                       colour = output_labels$colour_column), size = 3) +
      geom_text_repel(aes(label = text)) +
      ggplot2::scale_fill_manual(name = "Condition", values = class_colours) +
      ggplot2::guides(fill = guide_legend(override.aes = list(shape = 21,
                                                              colour = class_colours))) +
      ggplot2::scale_shape_manual(name = "Country", values = c(21, 22)) +
      ggplot2::guides(shape = guide_legend(override.aes = list(fill = c("black", "black")))) +
      ggplot2::scale_colour_manual(name = sub("_", " ", colour_column), values = point_border_colours) +
      ggplot2::guides(colour = guide_legend(override.aes = list(shape = 1))) +
      ggplot2::labs(title = title) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab) +
      labs(caption = paste(paste("Data dimension :", paste(dim(data), collapse = "x")), "\n",
                           group_counts_text),
           fill = fill_column)    
  } else{
    ggplot2::ggplot(dim_red_df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(ggplot2::aes(fill = output_labels$fill_column,
                                       shape = output_labels$country,
                                       colour = output_labels$Label), size = 3) +
      geom_text_repel(aes(label = text)) +
      ggplot2::guides(fill = guide_legend(override.aes = list(shape = 21))) +
      ggplot2::scale_shape_manual(name = "Country", values = c(21, 22)) +
      ggplot2::guides(shape = guide_legend(override.aes = list(fill = c("black", "black")))) +
      ggplot2::scale_colour_manual(name = "Condition", values = class_colours) +
      ggplot2::guides(colour = guide_legend(override.aes = list(shape = 1))) +
      ggplot2::labs(title = title) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab) +
      labs(caption = paste(paste("Data dimension :", paste(dim(data), collapse = "x")), "\n",
                           group_counts_text),
           fill = fill_column)    
  }

  if(is.na(dir_path)){
    dir_path <- paste0("prediction_pipeline/plots/dim_red/", fill_column)  
  }
  if(!dir.exists(dir_path)){
    dir.create(dir_path, recursive = TRUE)
  }
  file_name <- paste0(gsub(title, pattern = " |,", replacement = "-"), ".jpg")
  file_path <- paste(dir_path, file_name, sep = "/")
  ggplot2::ggsave(file_path, units = "cm", width = dimred_plot_width_cm)
  
  ###############
  #create biomarker box plots if best biomarkers only
  
  if(!is.na(best_features_file_path) & !is.na(dataset_replace_str)){
    data_to_plot <- data %>%
      rownames_to_column("sample_name") %>%
      inner_join(output_labels %>%
                   dplyr::select(c(Sample, Label)) %>%
                   rename(c("sample_name" = "Sample")))
    group_counts <- data_to_plot %>%
      group_by(Label) %>%
      summarize(n = n())
    
    data_to_plot <- data_to_plot %>%
      inner_join(group_counts) %>%
      mutate(Label = paste0(Label, " (", n, ") ")) %>%
      select(-c(n))
    
    data_to_plot <- data_to_plot %>%
      pivot_longer(!c(sample_name, Label), names_to = "biomarkers") %>%
      arrange(biomarkers, Label)
    
    data_to_plot_sub <- data_to_plot %>%
      group_by(biomarkers) %>%
      summarize(median = median(value), mean = mean(value)) %>%
      arrange(desc(median), desc(mean))
    
    # biomarker is already gene-name
    # adding this below mapping causes issues when multiple names present in protein names file
    
    # if(omics_type == "prot"){
    #   protein_names <- read.csv("data/proteomics/protein_names.csv") %>%
    #     dplyr::select(c(gene_name, protein_name, protein_id))
    #   
    #   data_to_plot_sub <- protein_names %>%
    #     right_join(data_to_plot_sub, by = c("gene_name" = "biomarkers")) %>%
    #     dplyr::rename(c("biomarkers" = "gene_name"))
    # }
    
    # biomarkers_rem <- data_to_plot_sub$biomarkers
    # 
    # biomarkers_ordered <- c(biomarkers_specific, biomarkers_rem)
    # 
    data_to_write <- as.data.frame(data_to_plot_sub)
    write.xlsx(data_to_write, biomarker_expr_file_path,
               sheetName = paste(omics_type, comparison, sep = "_"),
               col.names = TRUE, row.names = FALSE, append = TRUE)
    
    data_to_plot <- data_to_plot %>%
      mutate(biomarkers = factor(biomarkers, levels = data_to_plot_sub$biomarkers))
    plot_title <- sub("UMAP plot of ", "", title)
    if(norm == "log_tmm"){
      y_lab <- "Log TMM of expression level across samples"
    } else if(norm == "log_cpm"){
      y_lab <- "Log CPM of expression level across samples"
    } else if(norm == "quantile"){
      y_lab <- "Quantile normed expression level across samples"
    }
    
    num_biomarkers <- nrow(data_to_plot_sub)
    if(omics_type == "tra"){
      x_lab <- paste0("Transcripts (", num_biomarkers, ")")
    } else{
      x_lab <- paste0("Proteins (", num_biomarkers, ")")
    }
    
    plot <- ggplot() +
      aes(x = biomarkers, y = value, fill = Label) +
      geom_boxplot(data = data_to_plot) +
      xlab(x_lab) +
      ylab(y_lab) +
      ggtitle(plot_title) +
      labs(fill = "") +
      theme(
        strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
      )
    
    plot
    if(!dir.exists(box_plot_dir_path)){
      dir.create(box_plot_dir_path, recursive = TRUE)
    }
    file_name <- paste0(gsub(plot_title, pattern = " ", replacement = "-"), ".png")
    file_path <- paste(box_plot_dir_path, file_name, sep = "/")
    ggplot2::ggsave(file_path, units = "cm", width = plot_width_cm)
  }

}



#plot mutation
phenotype <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t")
summary(factor(phenotype$condition))

data_of_interest <- phenotype %>%
  filter(is.na(pre_post_modulator) | pre_post_modulator == 0) %>%
  filter(condition %in% c("CFRD", "IGT", "NGT")) %>%
  select(Sample, condition, country, patient_recruitment_year, age, age_group, sex, FEV1,
         mutation) %>%
  mutate(age = as.numeric(str_trim(age)), FEV1 = as.numeric(str_trim(FEV1)))
ggplot(data_of_interest, aes(x = mutation)) +
  geom_bar(position = "dodge", aes(fill = country)) +
  facet_wrap(~condition) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = rel(0.8)))
ggsave("prediction_pipeline/plots/other_features/mutation_cfrd_igt_ngt.jpg")


data_of_interest <- phenotype %>%
  select(Sample, condition, country, patient_recruitment_year, age, age_group, sex, FEV1,
         mutation) %>%
  mutate(age = as.numeric(str_trim(age)), FEV1 = as.numeric(str_trim(FEV1)))

ggplot(data_of_interest, aes(x = mutation)) +
  geom_bar(position = "dodge", aes(fill = country)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("prediction_pipeline/plots/other_features/mutation_country.jpg")


###compare transcripts after filter_by_expr b/w AU and DK cohorts
compare_transcripts_from_filter_by_expr <- function(comparison, classes, combat_seq = FALSE){
  if(combat_seq){
    data <- read.table("data/formatted/umi_counts_combat_seq.csv", header=TRUE, sep=",", row.names=1, skip=0,
                       nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
  } else{
    data <- read.table("data/formatted/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                       nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
  }
  phenotype <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t")
  output_labels <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes) %>%
    select(Sample, Label, country)
  
  au.output_labels <- output_labels %>%
    filter(country == "AU") %>%
    select(-c(country))
  dk.output_labels <- output_labels %>%
    filter(country == "DK") %>%
    select(-c(country))
  
  au.data <- data[, au.output_labels$Sample]
  dk.data <- data[, dk.output_labels$Sample]
  
  keep <- edgeR::filterByExpr(au.data, group = au.output_labels$Label)
  au.data <- au.data[keep, ]
  
  keep <- edgeR::filterByExpr(dk.data, group = dk.output_labels$Label)
  dk.data <- dk.data[keep, ]
  
  title <- paste("Filtered transcripts in samples from", sub("Vs", " Vs ", comparison))
  if(combat_seq){
    title <- paste(title, "after combat_seq")
  }
  
  ggvenn(list("AU" = rownames(au.data),
              "DK" = rownames(dk.data)),
         stroke_size = 0.1,
         set_name_size = 4,
         text_size = 3,
         fill_color = c("#F8766D", "#00BFC4")) +
    ggtitle(title) +
    theme(plot.title = element_text(vjust = - 20, hjust = 0.5))
  file_name <- paste0("prediction_pipeline/plots/filtered_transcripts/", 
                      comparison, "combat_seq", combat_seq, ".png")
  ggsave(file_name)
}

compare_transcripts_from_filter_by_expr(comparison = "CFRDVsIGT", 
                                        classes = c("CFRD", "IGT"), 
                                        combat_seq = FALSE)
compare_transcripts_from_filter_by_expr(comparison = "CFRDVsNGT", 
                                        classes = c("CFRD", "NGT"), 
                                        combat_seq = FALSE)
compare_transcripts_from_filter_by_expr(comparison = "IGTVsNGT", 
                                        classes = c("IGT", "NGT"), 
                                        combat_seq = FALSE)
compare_transcripts_from_filter_by_expr(comparison = "CFRDVsIGT", 
                                        classes = c("CFRD", "IGT"), 
                                        combat_seq = TRUE)
compare_transcripts_from_filter_by_expr(comparison = "CFRDVsNGT", 
                                        classes = c("CFRD", "NGT"), 
                                        combat_seq = TRUE)
compare_transcripts_from_filter_by_expr(comparison = "IGTVsNGT", 
                                        classes = c("IGT", "NGT"), 
                                        combat_seq = TRUE)



comparison = "CFRDVsIGT"
classes = c("IGT", "CFRD")
dataset_replace_str = "CF_EV_AU_adult_combat_logtmm_"
combat_seq = TRUE
best_features_file_path  = "data/selected_features/best_features_with_is_best.csv"

#compare best transcripts identified from AU with AU+DK filtered transcripts
compare_filtered_transcripts_with_best_biomarkers <- function(comparison, 
                                                              classes, 
                                                              dataset_replace_str,
                                                              combat_seq = FALSE,
                                                              best_features_file_path  = "data/selected_features/best_features_with_is_best.csv"){
  
  best_features <- read.csv(best_features_file_path)  
  best_features_sub <- best_features %>%
    mutate(dataset_id = gsub(dataset_replace_str, "", dataset_id)) %>%
    filter(is_best == 1, dataset_id == comparison)
  biomarkers <- strsplit(best_features_sub$biomarkers, split = "|", fixed = TRUE)[[1]]  
  
  if(combat_seq){
    data <- read.table("data/formatted/umi_counts_combat_seq.csv", header=TRUE, sep=",", row.names=1, skip=0,
                       nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
  } else{
    data <- read.table("data/formatted/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                       nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
  }
  phenotype <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t")
  output_labels <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes) %>%
    select(Sample, Label)
  
  data <- data[, output_labels$Sample]
  
  keep <- edgeR::filterByExpr(data, group = output_labels$Label)
  data <- data[keep, ]
  
  features_with_slash <- rownames(data)[grepl("/", rownames(data), fixed = TRUE)] 
  for(f in features_with_slash){
    f_replaced <- gsub("/|-", ".", f) 
    if(f_replaced %in% biomarkers){
      biomarkers[biomarkers == f_replaced] = f
    }
  }
  biomarkers <- gsub(".", "-", biomarkers, fixed = TRUE)

  title <- paste("Filtered transcripts against best biomarkers",
                 gsub("_$", "", gsub("CF_EV_", "", dataset_replace_str)),
                 " in samples from", sub("Vs", " Vs ", comparison))
  
  ggvenn(list("AU + DK filtered" = rownames(data),
              "Best biomarkers" = biomarkers),
         stroke_size = 0.1,
         set_name_size = 4,
         text_size = 3,
         fill_color = c("#F8766D", "#00BFC4")) +
    ggtitle(title) +
    theme(plot.title = element_text(vjust = - 20, hjust = 0.5))
  file_name <- paste0("prediction_pipeline/plots/filtered_transcripts_with_best_biomarkers/", 
                      comparison, "_",
                      gsub("_$", "", gsub("CF_EV_", "", dataset_replace_str)), 
                      ".png")
  ggsave(file_name, units = "cm", height = 20, width = 30)
}


compare_filtered_transcripts_with_best_biomarkers(comparison = "CFRDVsIGT",
                                                  classes = c("IGT", "CFRD"),
                                                  dataset_replace_str = "CF_EV_AU_adult_combat_logtmm_",
                                                  combat_seq = TRUE,
                                                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv")
compare_filtered_transcripts_with_best_biomarkers(comparison = "CFRDVsNGT",
                                                  classes = c("NGT", "CFRD"),
                                                  dataset_replace_str = "CF_EV_AU_adult_combat_logtmm_",
                                                  combat_seq = TRUE,
                                                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv")
compare_filtered_transcripts_with_best_biomarkers(comparison = "IGTVsNGT",
                                                  classes = c("NGT", "IGT"),
                                                  dataset_replace_str = "CF_EV_AU_adult_combat_logtmm_",
                                                  combat_seq = TRUE,
                                                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv")

compare_filtered_transcripts_with_best_biomarkers(comparison = "CFRDVsIGT",
                                                  classes = c("IGT", "CFRD"),
                                                  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
                                                  combat_seq = FALSE,
                                                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv")
compare_filtered_transcripts_with_best_biomarkers(comparison = "CFRDVsNGT",
                                                  classes = c("NGT", "CFRD"),
                                                  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
                                                  combat_seq = FALSE,
                                                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv")
compare_filtered_transcripts_with_best_biomarkers(comparison = "IGTVsNGT",
                                                  classes = c("NGT", "IGT"),
                                                  dataset_replace_str = "CF_EV_AU_adult_logtmm_",
                                                  combat_seq = FALSE,
                                                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv")



compare_filtered_transcripts_with_best_biomarkers(comparison = "CFRDVsIGT",
                                                  classes = c("IGT", "CFRD"),
                                                  dataset_replace_str = "CF_EV_AU_adult_combat_logcpm_",
                                                  combat_seq = TRUE,
                                                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv")
compare_filtered_transcripts_with_best_biomarkers(comparison = "CFRDVsNGT",
                                                  classes = c("NGT", "CFRD"),
                                                  dataset_replace_str = "CF_EV_AU_adult_combat_logcpm_",
                                                  combat_seq = TRUE,
                                                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv")
compare_filtered_transcripts_with_best_biomarkers(comparison = "IGTVsNGT",
                                                  classes = c("NGT", "IGT"),
                                                  dataset_replace_str = "CF_EV_AU_adult_combat_logcpm_",
                                                  combat_seq = TRUE,
                                                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv")

compare_filtered_transcripts_with_best_biomarkers(comparison = "CFRDVsIGT",
                                                  classes = c("IGT", "CFRD"),
                                                  dataset_replace_str = "CF_EV_AU_adult_logcpm_",
                                                  combat_seq = FALSE,
                                                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv")
compare_filtered_transcripts_with_best_biomarkers(comparison = "CFRDVsNGT",
                                                  classes = c("NGT", "CFRD"),
                                                  dataset_replace_str = "CF_EV_AU_adult_logcpm_",
                                                  combat_seq = FALSE,
                                                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv")
compare_filtered_transcripts_with_best_biomarkers(comparison = "IGTVsNGT",
                                                  classes = c("NGT", "IGT"),
                                                  dataset_replace_str = "CF_EV_AU_adult_logcpm_",
                                                  combat_seq = FALSE,
                                                  best_features_file_path  = "data/selected_features/best_features_with_is_best.csv")

