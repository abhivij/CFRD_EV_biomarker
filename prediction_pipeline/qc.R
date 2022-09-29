#obtain count of number of samples in different groups like CF, CFRD, IGT, ...
#create boxplots and dim red plots of transcripts to check distribution across different cohorts

library(tidyverse)
library(ggrepel)
library(umap)

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


#data and groups expected to be matching
#should be ensured from the caller function
create_dim_red_plots <- function(data, groups,
                                 title_prefix = "",
                                 dim_red = "UMAP",
                                 perplexity = 5,
                                 colour_label = "Condition",
                                 shownames = FALSE){
  
  groups_modified <- groups
  summa <- summary(factor(groups))
  
  groups_modified <- gsub(pattern = names(summa)[1], 
                          replacement = paste0(names(summa)[1], "(", summa[1], ")"),
                          groups_modified)
  groups_modified <- gsub(pattern = names(summa)[2], 
                          replacement = paste0(names(summa)[2], "(", summa[2], ")"),
                          groups_modified)
  groups <- factor(groups_modified)
  
  if(shownames){
    text <- rownames(data)
  } else{
    text <- ""
  }
  set.seed(1)
  if(dim_red == "tSNE"){
    result <- Rtsne::Rtsne(data, perplexity = perplexity)
    dim_red_df <- data.frame(x = result$Y[,1], y = result$Y[,2], 
                             Colour = groups, 
                             Sample = text)    
    xlab <- "tSNE 1"
    ylab <- "tSNE 2"
  } else if(dim_red == "UMAP"){
    print(dim(data)[1])
    # n_neighbors <- max(floor(dim(data)[1] / 4), 2)
    # print(n_neighbors)
    # result <- umap(data, n_neighbors = n_neighbors)
    result <- umap(data)
    dim_red_df <- data.frame(x = result$layout[,1], y = result$layout[,2], 
                             Colour = groups, 
                             Sample = text)  
    xlab <- "UMAP 1"
    ylab <- "UMAP 2"
  }
  
  title <- paste(title_prefix, dim_red)
  
  if (shownames) {
    dim_red_plot <- ggplot2::ggplot(dim_red_df,
                                    ggplot2::aes(x = x, y = y, colour = Colour)) +
      ggplot2::geom_point() +
      geom_text_repel(aes(label=Sample)) +
      ggplot2::labs(title = title, colour = colour_label) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab) +
      labs(caption = paste("Data dimension :", paste(dim(data), collapse = "x")))
    dim_red_plot
  } else {
    dim_red_plot <- ggplot2::ggplot(dim_red_df) +
      ggplot2::geom_point(ggplot2::aes(x = x, y = y, colour = Colour)) +
      ggplot2::labs(title = title, colour = colour_label) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab) +
      labs(caption = paste("Data dimension :", paste(dim(data), collapse = "x")))
    dim_red_plot
  }
  
  dir_path <- "prediction_pipeline/plots"
  file_name <- paste0(gsub(title, pattern = " ", replacement = "-"), ".jpg")
  file_path <- paste(dir_path, file_name, sep = "/")
  ggplot2::ggsave(file_path, dim_red_plot, units = "cm", width = 30)
}


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
