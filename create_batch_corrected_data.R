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

create_combat_files <- function(comparison, classes, norm = "log_tmm"){
  
  data <- read.table("data/formatted/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                     nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
  phenotype <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t")  
  
  output_labels <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes) %>%
    dplyr::select(Sample, Label, country, age_group)
  
  data <- data[, output_labels$Sample]
  
  #currently data format : (transcripts x samples)
  keep <- edgeR::filterByExpr(data, group = output_labels$Label)
  data <- data[keep, ]
  
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
  }

  #perform combat
  data_of_interest <- as.data.frame(as.matrix(data))
  
  all.equal(colnames(data_of_interest), output_labels$Sample)
  
  data_of_interest.combat = ComBat(dat=data_of_interest, 
                                   batch=output_labels$country)
  data_of_interest.combat <- as.data.frame(as.matrix(data_of_interest.combat))
  
  data <- data_of_interest.combat

  file_path <- paste0("data/formatted/", comparison, "_umi_counts_combat_processed.csv")
  write.csv(data, file_path)  
}

create_combat_files(comparison = "CFRDVsIGT", 
                    classes = c("CFRD", "IGT"))
create_combat_files(comparison = "CFRDVsNGT", 
                    classes = c("CFRD", "NGT"))
create_combat_files(comparison = "IGTVsNGT", 
                    classes = c("IGT", "NGT"))
