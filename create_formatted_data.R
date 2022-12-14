#author : Abhishek Vijayan

library(tidyverse)
library(readxl)

library(sva)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/CysticFibrosisGroup/ExoCF/CFRD_EV_biomarker/"
setwd(base_dir)

#Each of CFRD_100_1, CFRD_100_2, ... corresponds to results of different quantification jobs
#   from Qiagen Geneglobe
for(i in c(1:11)){
  data_file_name <- paste0("CFRD_100_", i, ".xlsx")
  print(paste("data", data_file_name, sep = "/"))
  
  data <- read_excel(paste("data", data_file_name, sep = "/"), sheet = "miRNA_piRNA") %>%
    column_to_rownames("miRNA") %>%
    select(ends_with("UMIs"))
  colnames(data) <- gsub("-UMIs", "", colnames(data))
  
  pirna_header <- which(rownames(data) == "piRNA")
  
  mirna_data <- data[1:pirna_header-1, ]
  pirna_data <- data[(pirna_header+1):dim(data)[1], ]
  rownames(pirna_data) <- gsub("/gb", "", rownames(pirna_data), fixed = TRUE)
  rownames(pirna_data) <- gsub("/Homo", "", rownames(pirna_data), fixed = TRUE)

  pirna_mapping <- data.frame(combined_id = rownames(pirna_data)) %>%
    separate(col = combined_id, into = c("piRNA_name", "accession"), sep = "/")
  
  rownames(pirna_data) <- sapply(rownames(pirna_data), FUN = 
                                   function(x){
                                     strsplit(x, split = "/", fixed = TRUE)[[1]][1]
                                   })
  
  mirna_data <- mirna_data %>% rownames_to_column(var = "transcript")
  pirna_data <- pirna_data %>% rownames_to_column(var = "transcript")
  if(i == 1){
    mirna_data_all <- mirna_data
    pirna_data_all <- pirna_data
    pirna_mapping_all <- pirna_mapping
    quantification_batch <- data.frame("sample_long_name" = colnames(data),
                                       "quant_batch" = i)
  } else{
    mirna_data_all <- mirna_data_all %>%
      full_join(mirna_data)
    pirna_data_all <- pirna_data_all %>%
      full_join(pirna_data)
    pirna_mapping_all <- rbind(pirna_mapping_all, pirna_mapping) %>%
      distinct()
    quantification_batch <- rbind(quantification_batch,
                                  data.frame("sample_long_name" = colnames(data),
                                             "quant_batch" = i))
    
  }

}      

dim(mirna_data_all)
# [1] 747 273
dim(pirna_data_all)
# [1] 655 273

combined_data <- rbind(mirna_data_all, pirna_data_all)
dim(combined_data)
# [1] 1402  273
sum(is.na(combined_data))
# [1] 193311
sum(combined_data[,2:273], na.rm = TRUE)
# [1] 345205621

combined_data[is.na(combined_data)] <- 0

sum(combined_data[,2:273], na.rm = TRUE)
# [1] 345205621
sum(is.na(combined_data))
# [1] 0

dim(pirna_mapping_all)
# [1] 655   2

dim(quantification_batch)
# [1] 272   2


write.csv(combined_data, "data/formatted/umi_counts.csv", row.names = FALSE)
write.csv(pirna_mapping_all, "data/formatted/pirna_mapping.csv", row.names = FALSE)
write.csv(quantification_batch, "data/formatted/quantification_batch.csv", row.names = FALSE)


#create combat-seq transformed data file
data <- read.table("data/formatted/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                   nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
phenotype <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t")
output_labels <- phenotype %>%
  select(Sample, country)
data.combat_seq <- ComBat_seq(counts = as.matrix(data),
                              batch = factor(output_labels$country))
write.csv(data.combat_seq, "data/formatted/umi_counts_combat_seq.csv")
