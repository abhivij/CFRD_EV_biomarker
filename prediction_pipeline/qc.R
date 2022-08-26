#obtain count of number of samples in different groups like CF, CFRD, IGT, ...
#create boxplots and dim red plots of transcripts to check distribution across different cohorts

library(tidyverse)

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