#obtain count of number of samples in different groups like CF, CFRD, IGT, ...
#create boxplots and dim red plots of transcripts to check distribution across different cohorts

library(tidyverse)
library(ggrepel)
library(umap)
library(ggvenn)
library(sva)


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

comparison = "CFRDVsNGT"
classes = c("NGT", "CFRD")
use_best_features = TRUE
norm = "log_tmm"
best_features_file_path = "data/selected_features/best_features_with_is_best.csv"
perform_filter = TRUE
train_cohort_country = "AU"
plot_width_cm = 20

create_transcript_box_plots <- function(comparison,
                                        classes,
                                        use_best_features = TRUE,
                                        train_cohort_country = "AU",
                                        norm = "none",
                                        best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                                        perform_filter = TRUE,
                                        plot_width_cm = 20){
  
  data <- read.table("data/formatted/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                     nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
  phenotype <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t")
  phenotype <- phenotype %>% filter(age_group == "adult")
  
  test_cohort_country <- ifelse(train_cohort_country == "AU", "DK", "AU") 
  
  train.output_labels <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes, country == train_cohort_country) %>%
    dplyr::select(Sample, Label) %>%
    dplyr::mutate(Label = factor(Label)) %>%
    arrange(Label, Sample)
  train.tra_data <- data[, train.output_labels$Sample]
  
  test.output_labels <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes, country == test_cohort_country) %>%
    dplyr::select(Sample, Label) %>%
    dplyr::mutate(Label = factor(Label)) %>%    
    arrange(Label, Sample)
  test.tra_data <- data[, test.output_labels$Sample]

  #currently data format : (transcripts x samples)
  if(perform_filter){
    keep <- edgeR::filterByExpr(train.tra_data, group = train.output_labels$Label)
    train.tra_data <- train.tra_data[keep, ]
    test.tra_data <- test.tra_data[keep, ]  
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
  }
  train.tra_data <- as.data.frame(t(as.matrix(train.tra_data)))
  test.tra_data <- as.data.frame(t(as.matrix(test.tra_data)))  
  #now data : (samples x transcripts)

  if(use_best_features){
    best_features <- read.csv(best_features_file_path)  
    
    dataset_replace_str <- paste0("CF_EV_", train_cohort_country, 
                                  "_adult_logtmm_")
    
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
    
    train.tra_data <- train.tra_data[, biomarkers]
    test.tra_data <- test.tra_data[, biomarkers]
  }
  num_transcripts <- dim(train.tra_data)[2]
  
  data_to_plot <- rbind(
    cbind(train.tra_data, 
          "label" = train.output_labels$Label,
          "country" = train_cohort_country),
    cbind(test.tra_data,
          "label" = test.output_labels$Label,
          "country" = test_cohort_country)
  ) %>%
    mutate(label = paste(country, label, sep = "_")) %>%
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
  
  y_lab <- ifelse(norm == "log_tmm",
                  "Log TMM of expression level across samples",
                  "Expression level across samples")
  
  plot_title <- ifelse(norm == "log_tmm",
                         "Log tmm expression of ",
                         "Non-normalized expression of ")
  if(use_best_features){
    plot_title <- paste0(plot_title, "best transcripts from ", train_cohort_country)
  } else{
    plot_title <- paste0(plot_title, "all transcripts with params from ", train_cohort_country)
  }
  plot_title <- paste0(plot_title, " in classes ", paste0(rev(classes), collapse = ","))

  plot <- ggplot(data_to_plot, aes(x = transcripts, y = value)) +
    geom_boxplot(aes(fill = label)) +
    xlab(paste0("Transcripts (", num_transcripts, ")")) +
    ylab(y_lab) +
    ggtitle(plot_title) +
    labs(fill = "") 
  if(use_best_features){
    plot <- plot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  } else {
    plot <- plot + theme(axis.text.x = element_blank())
  }
  plot
  
  dir_path <- "prediction_pipeline/plots/box_plots"
  if(!dir.exists(dir_path)){
    dir.create(dir_path, recursive = TRUE)
  }
  file_name <- paste0(gsub(plot_title, pattern = " ", replacement = "-"), ".png")
  file_path <- paste(dir_path, file_name, sep = "/")
  ggplot2::ggsave(file_path, units = "cm", width = plot_width_cm)
}

create_transcript_box_plots(comparison = "CFRDVsIGT",
                            classes = c("IGT", "CFRD"),
                            use_best_features = TRUE,
                            train_cohort_country = "AU",
                            norm = "log_tmm",
                            best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                            perform_filter = TRUE)
create_transcript_box_plots(comparison = "CFRDVsIGT",
                            classes = c("IGT", "CFRD"),
                            use_best_features = TRUE,
                            train_cohort_country = "DK",
                            norm = "log_tmm",
                            best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                            perform_filter = TRUE)
create_transcript_box_plots(comparison = "CFRDVsIGT",
                            classes = c("IGT", "CFRD"),
                            use_best_features = TRUE,
                            train_cohort_country = "AU",
                            norm = "none",
                            best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                            perform_filter = TRUE)
create_transcript_box_plots(comparison = "CFRDVsIGT",
                            classes = c("IGT", "CFRD"),
                            use_best_features = TRUE,
                            train_cohort_country = "DK",
                            norm = "none",
                            best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                            perform_filter = TRUE)
create_transcript_box_plots(comparison = "CFRDVsIGT",
                            classes = c("IGT", "CFRD"),
                            use_best_features = FALSE, plot_width = 80,
                            train_cohort_country = "AU",
                            norm = "log_tmm",
                            best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                            perform_filter = TRUE)
create_transcript_box_plots(comparison = "CFRDVsIGT",
                            classes = c("IGT", "CFRD"),
                            use_best_features = FALSE, plot_width = 80,
                            train_cohort_country = "DK",
                            norm = "log_tmm",
                            best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                            perform_filter = TRUE)
create_transcript_box_plots(comparison = "CFRDVsIGT",
                            classes = c("IGT", "CFRD"),
                            use_best_features = FALSE, plot_width = 80,
                            train_cohort_country = "AU",
                            norm = "none",
                            best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                            perform_filter = TRUE)
create_transcript_box_plots(comparison = "CFRDVsIGT",
                            classes = c("IGT", "CFRD"),
                            use_best_features = FALSE, plot_width = 80,
                            train_cohort_country = "DK",
                            norm = "none",
                            best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                            perform_filter = TRUE)

create_transcript_box_plots(comparison = "CFRDVsNGT",
                            classes = c("NGT", "CFRD"),
                            use_best_features = TRUE,
                            train_cohort_country = "AU",
                            norm = "log_tmm",
                            best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                            perform_filter = TRUE)
create_transcript_box_plots(comparison = "CFRDVsNGT",
                            classes = c("NGT", "CFRD"),
                            use_best_features = TRUE,
                            train_cohort_country = "DK",
                            norm = "log_tmm",
                            best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                            perform_filter = TRUE)
create_transcript_box_plots(comparison = "CFRDVsNGT",
                            classes = c("NGT", "CFRD"),
                            use_best_features = TRUE,
                            train_cohort_country = "AU",
                            norm = "none",
                            best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                            perform_filter = TRUE)
create_transcript_box_plots(comparison = "CFRDVsNGT",
                            classes = c("NGT", "CFRD"),
                            use_best_features = TRUE,
                            train_cohort_country = "DK",
                            norm = "none",
                            best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                            perform_filter = TRUE)
create_transcript_box_plots(comparison = "CFRDVsNGT",
                            classes = c("NGT", "CFRD"),
                            use_best_features = FALSE, plot_width = 80,
                            train_cohort_country = "AU",
                            norm = "log_tmm",
                            best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                            perform_filter = TRUE)
create_transcript_box_plots(comparison = "CFRDVsNGT",
                            classes = c("NGT", "CFRD"),
                            use_best_features = FALSE, plot_width = 80,
                            train_cohort_country = "DK",
                            norm = "log_tmm",
                            best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                            perform_filter = TRUE)
create_transcript_box_plots(comparison = "CFRDVsNGT",
                            classes = c("NGT", "CFRD"),
                            use_best_features = FALSE, plot_width = 80,
                            train_cohort_country = "AU",
                            norm = "none",
                            best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                            perform_filter = TRUE)
create_transcript_box_plots(comparison = "CFRDVsNGT",
                            classes = c("NGT", "CFRD"),
                            use_best_features = FALSE, plot_width = 80,
                            train_cohort_country = "DK",
                            norm = "none",
                            best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                            perform_filter = TRUE)

create_transcript_box_plots(comparison = "IGTVsNGT",
                            classes = c("NGT", "IGT"),
                            use_best_features = TRUE,
                            train_cohort_country = "AU",
                            norm = "log_tmm",
                            best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                            perform_filter = TRUE)
create_transcript_box_plots(comparison = "IGTVsNGT",
                            classes = c("NGT", "IGT"),
                            use_best_features = TRUE,
                            train_cohort_country = "DK",
                            norm = "log_tmm",
                            best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                            perform_filter = TRUE)
create_transcript_box_plots(comparison = "IGTVsNGT",
                            classes = c("NGT", "IGT"),
                            use_best_features = TRUE,
                            train_cohort_country = "AU",
                            norm = "none",
                            best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                            perform_filter = TRUE)
create_transcript_box_plots(comparison = "IGTVsNGT",
                            classes = c("NGT", "IGT"),
                            use_best_features = TRUE,
                            train_cohort_country = "DK",
                            norm = "none",
                            best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                            perform_filter = TRUE)
create_transcript_box_plots(comparison = "IGTVsNGT",
                            classes = c("NGT", "IGT"),
                            use_best_features = FALSE, plot_width = 80,
                            train_cohort_country = "AU",
                            norm = "log_tmm",
                            best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                            perform_filter = TRUE)
create_transcript_box_plots(comparison = "IGTVsNGT",
                            classes = c("NGT", "IGT"),
                            use_best_features = FALSE, plot_width = 80,
                            train_cohort_country = "DK",
                            norm = "log_tmm",
                            best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                            perform_filter = TRUE)
create_transcript_box_plots(comparison = "IGTVsNGT",
                            classes = c("NGT", "IGT"),
                            use_best_features = FALSE, plot_width = 80,
                            train_cohort_country = "AU",
                            norm = "none",
                            best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                            perform_filter = TRUE)
create_transcript_box_plots(comparison = "IGTVsNGT",
                            classes = c("NGT", "IGT"),
                            use_best_features = FALSE, plot_width = 80,
                            train_cohort_country = "DK",
                            norm = "none",
                            best_features_file_path = "data/selected_features/best_features_with_is_best.csv",
                            perform_filter = TRUE)



#############
#create dim_red plots with both AU and DK cohort in one figure
#this is to check if there are separate clusters for AU and DK within each of the pairs of CFRD, IGT, NGT
# comparison = "CFRDVsIGT" #this will be NA if 3 classes are involved
#                          #in that case non-modulator sample filtering will be seprately done too
# classes = c("NGT", "CFRD")
# 
# comparison <- NA
# classes = c("CFRD", "IGT", "NGT")
# dim_red = "UMAP"
# shownames = FALSE
# norm = "log_tmm"

comparison = "CFRDVsIGT"
classes = c("CFRD", "IGT")
class_colours = c("red", "orange")
dim_red = "UMAP"
norm = "log"
shownames = FALSE
fill_column = "mutation"
fill_column = NA

combat_seq = FALSE
combat_seq_specify_group = FALSE

# comparison = "PreModulatorVsPostModulator"
# classes = c("PreModulator", "PostModulator")
# class_colours = c("brown", "skyblue")
# dim_red = "UMAP"
# norm = "log_tmm"
# shownames = FALSE

create_dim_red_plots <- function(comparison, classes,
                                 class_colours,
                                 dim_red, norm,
                                 fill_column = NA,
                                 combat_seq = FALSE, combat_seq_specify_group = FALSE,
                                 shownames = FALSE,
                                 perform_filter = TRUE){
  data <- read.table("data/formatted/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                     nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
  phenotype <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t")
  
  title <- paste0(dim_red, " plot of ", paste(classes, collapse = ", "), " samples from ", norm, " data")
  
  if(!is.na(comparison)){
    output_labels <- phenotype %>%
      rename("Label" = comparison)
  } else{
    #this case is to include greater than 2 conditions in the plot
    output_labels <- phenotype %>%
      rename("Label" = "condition") %>%
      #ensure that samples on modulator are not chosen
      #using comparison field does not require the below, because comparison field creation 
      #                 incorporates this filter
      filter(is.na(pre_post_modulator) | pre_post_modulator == 0)
  }
  output_labels <- output_labels %>%
    filter(Label %in% classes) %>%
    dplyr::select(Sample, Label, country, age_group, mutation, 
                  sex, cohort, patient_recruitment_year,
                  seq_plate, seq_miR_library_quality,
                  quant_batch) %>%
    dplyr::mutate(Label = factor(Label)) %>%
    dplyr::mutate(seq_plate = factor(seq_plate), quant_batch = factor(quant_batch)) %>%
    arrange(Label, Sample)
  
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
  }
  data <- as.data.frame(t(as.matrix(data)))
  
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
  
  if(is.na(fill_column)){
    ggplot2::ggplot(dim_red_df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(ggplot2::aes(fill = output_labels$Label,
                                       shape = output_labels$country,
                                       colour = output_labels$age_group), size = 3) +
      geom_text_repel(aes(label = text)) +
      ggplot2::scale_fill_manual(name = "Condition", values = class_colours) +
      ggplot2::guides(fill = guide_legend(override.aes = list(shape = 21,
                                                              colour = class_colours))) +
      ggplot2::scale_shape_manual(name = "Country", values = c(21, 22)) +
      ggplot2::guides(shape = guide_legend(override.aes = list(fill = c("black", "black")))) +
      ggplot2::scale_colour_manual(name = "Age group", values = c("black", "green")) +
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

  
  dir_path <- paste0("prediction_pipeline/plots/dim_red/", fill_column)
  if(!dir.exists(dir_path)){
    dir.create(dir_path, recursive = TRUE)
  }
  file_name <- paste0(gsub(title, pattern = " |,", replacement = "-"), ".jpg")
  file_path <- paste(dir_path, file_name, sep = "/")
  ggplot2::ggsave(file_path, units = "cm", width = 30)
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
