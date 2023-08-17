#to create data files for the analysis after including proteomics

library(tidyverse)
library(readxl)

library(sva)

library(ggvenn)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/CysticFibrosisGroup/ExoCF/CFRD_EV_biomarker/"
setwd(base_dir)

#create combat corrected datasets

comparison = "CFRDVsIGT"
classes = c("CFRD", "IGT")
norm = "log_tmm"

data_file_path = "data/proteomics/imputed_combined.csv"
phenotype_file_path = "data/formatted/prot_phenotype.txt"
combattwice = TRUE
norm = "quantile"
perform_filter = FALSE

file_dir_path = "data/formatted/proteomics/"
file_suffix = "_imputed_combined_quantile_combattwice.csv"




create_combat_files <- function(comparison, classes, 
                                norm = "log_tmm",
                                data_file_path = "data/formatted/umi_counts.csv",
                                phenotype_file_path = "data/formatted/phenotype.txt",
                                combattwice = FALSE,
                                perform_filter = TRUE,
                                file_dir_path = "data/formatted/",
                                file_suffix = "_umi_counts_combat_processed.csv"){
  
  data <- read.table(data_file_path, header=TRUE, sep=",", row.names=1, skip=0,
                     nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
  phenotype <- read.table(phenotype_file_path, header=TRUE, sep="\t")  
  
  output_labels <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes)
  if(combattwice){
    output_labels <- output_labels %>%
      dplyr::select(Sample, Label, country, age_group, mq_batch)
  } else{
    output_labels <- output_labels %>%
      dplyr::select(Sample, Label, country, age_group)
  }

    
  
  data <- data[, output_labels$Sample]
  
  #currently data format : (transcripts x samples)
  if(perform_filter){
    keep <- edgeR::filterByExpr(data, group = output_labels$Label)
    data <- data[keep, ]
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

  #perform combat
  data_of_interest <- as.data.frame(as.matrix(data))
  
  all.equal(colnames(data_of_interest), output_labels$Sample)
  
  data_of_interest.combat = ComBat(dat=data_of_interest, 
                                   batch=output_labels$country)
  if(combattwice){
    data_of_interest.combat = ComBat(dat=data_of_interest.combat, 
                                     batch=output_labels$mq_batch)
  }
  data_of_interest.combat <- as.data.frame(as.matrix(data_of_interest.combat))
  data <- data_of_interest.combat

  file_path <- paste0(file_dir_path, comparison, file_suffix)
  write.csv(data, file_path)  
}

create_combat_files(comparison = "CFRDVsIGT", 
                    classes = c("CFRD", "IGT"))
create_combat_files(comparison = "CFRDVsNGT", 
                    classes = c("CFRD", "NGT"))
create_combat_files(comparison = "IGTVsNGT", 
                    classes = c("IGT", "NGT"))


create_combat_files(comparison = "CFRDVsIGT", 
                    classes = c("CFRD", "IGT"), 
                    norm = "quantile",
                    data_file_path = "data/proteomics/imputed_combined.csv",
                    phenotype_file_path = "data/formatted/prot_phenotype.txt",
                    combattwice = FALSE,
                    perform_filter = FALSE,
                    file_dir_path = "data/formatted/proteomics/",
                    file_suffix = "_imputed_combined_quantile_combat.csv")
create_combat_files(comparison = "CFRDVsNGT", 
                    classes = c("CFRD", "NGT"), 
                    norm = "quantile",
                    data_file_path = "data/proteomics/imputed_combined.csv",
                    phenotype_file_path = "data/formatted/prot_phenotype.txt",
                    combattwice = FALSE,
                    perform_filter = FALSE,
                    file_dir_path = "data/formatted/proteomics/",
                    file_suffix = "_imputed_combined_quantile_combat.csv")
create_combat_files(comparison = "IGTVsNGT", 
                    classes = c("IGT", "NGT"), 
                    norm = "quantile",
                    data_file_path = "data/proteomics/imputed_combined.csv",
                    phenotype_file_path = "data/formatted/prot_phenotype.txt",
                    combattwice = FALSE,
                    perform_filter = FALSE,
                    file_dir_path = "data/formatted/proteomics/",
                    file_suffix = "_imputed_combined_quantile_combat.csv")


create_combat_files(comparison = "CFRDVsIGT", 
                    classes = c("CFRD", "IGT"), 
                    norm = "quantile",
                    data_file_path = "data/proteomics/imputed_combined.csv",
                    phenotype_file_path = "data/formatted/prot_phenotype.txt",
                    combattwice = TRUE,
                    perform_filter = FALSE,
                    file_dir_path = "data/formatted/proteomics/",
                    file_suffix = "_imputed_combined_quantile_combattwice.csv")
create_combat_files(comparison = "CFRDVsNGT", 
                    classes = c("CFRD", "NGT"), 
                    norm = "quantile",
                    data_file_path = "data/proteomics/imputed_combined.csv",
                    phenotype_file_path = "data/formatted/prot_phenotype.txt",
                    combattwice = TRUE,
                    perform_filter = FALSE,
                    file_dir_path = "data/formatted/proteomics/",
                    file_suffix = "_imputed_combined_quantile_combattwice.csv")
create_combat_files(comparison = "IGTVsNGT", 
                    classes = c("IGT", "NGT"), 
                    norm = "quantile",
                    data_file_path = "data/proteomics/imputed_combined.csv",
                    phenotype_file_path = "data/formatted/prot_phenotype.txt",
                    combattwice = TRUE,
                    perform_filter = FALSE,
                    file_dir_path = "data/formatted/proteomics/",
                    file_suffix = "_imputed_combined_quantile_combattwice.csv")


create_combat_files(comparison = "CFRDVsIGT", 
                    classes = c("CFRD", "IGT"), 
                    norm = "quantile",
                    data_file_path = "data/proteomics/imputed_main.csv",
                    phenotype_file_path = "data/formatted/prot_phenotype_only_main.txt",
                    combattwice = FALSE,
                    perform_filter = FALSE,
                    file_dir_path = "data/formatted/proteomics/",
                    file_suffix = "_imputed_main_quantile_combat.csv")
create_combat_files(comparison = "CFRDVsNGT", 
                    classes = c("CFRD", "NGT"), 
                    norm = "quantile",
                    data_file_path = "data/proteomics/imputed_main.csv",
                    phenotype_file_path = "data/formatted/prot_phenotype_only_main.txt",
                    combattwice = FALSE,
                    perform_filter = FALSE,
                    file_dir_path = "data/formatted/proteomics/",
                    file_suffix = "_imputed_main_quantile_combat.csv")
create_combat_files(comparison = "IGTVsNGT", 
                    classes = c("IGT", "NGT"), 
                    norm = "quantile",
                    data_file_path = "data/proteomics/imputed_main.csv",
                    phenotype_file_path = "data/formatted/prot_phenotype_only_main.txt",
                    combattwice = FALSE,
                    perform_filter = FALSE,
                    file_dir_path = "data/formatted/proteomics/",
                    file_suffix = "_imputed_main_quantile_combat.csv")



create_combat_files(comparison = "CFRDVsIGT", 
                    classes = c("CFRD", "IGT"), 
                    norm = "quantile",
                    data_file_path = "data/proteomics/data_333samples_imputed_mf.csv",
                    phenotype_file_path = "data/formatted/prot_phenotype_333.txt",
                    combattwice = FALSE,
                    perform_filter = FALSE,
                    file_dir_path = "data/formatted/proteomics/",
                    file_suffix = "_imputed333_mf_quantile_combat.csv")
create_combat_files(comparison = "CFRDVsNGT", 
                    classes = c("CFRD", "NGT"), 
                    norm = "quantile",
                    data_file_path = "data/proteomics/data_333samples_imputed_mf.csv",
                    phenotype_file_path = "data/formatted/prot_phenotype_333.txt",
                    combattwice = FALSE,
                    perform_filter = FALSE,
                    file_dir_path = "data/formatted/proteomics/",
                    file_suffix = "_imputed333_mf_quantile_combat.csv")
create_combat_files(comparison = "IGTVsNGT", 
                    classes = c("IGT", "NGT"), 
                    norm = "quantile",
                    data_file_path = "data/proteomics/data_333samples_imputed_mf.csv",
                    phenotype_file_path = "data/formatted/prot_phenotype_333.txt",
                    combattwice = FALSE,
                    perform_filter = FALSE,
                    file_dir_path = "data/formatted/proteomics/",
                    file_suffix = "_imputed333_mf_quantile_combat.csv")
