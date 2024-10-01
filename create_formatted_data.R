#author : Abhishek Vijayan

library(tidyverse)
library(readxl)

library(sva)
# library(Seurat)

library(ggvenn)
library(ggplot2)

library(ComplexHeatmap)
library(viridis)
library(RColorBrewer)


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



#create combined filter (AU + DK) filtered together + combatseq datasets

create_combined_filter_with_combat_seq_files <- function(comparison, classes){
  data <- read.table("data/formatted/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                     nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
  phenotype <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t")  
  
  output_labels <- phenotype %>%
    rename("Label" = comparison) %>%
    filter(Label %in% classes, age_group == "adult") %>%
    dplyr::select(Sample, Label, country, age_group)
  
  data <- data[, output_labels$Sample]

  keep <- edgeR::filterByExpr(data, group = output_labels$Label)
  data <- data[keep, ]

  data.combat_seq <- ComBat_seq(counts = as.matrix(data),
                                batch = factor(output_labels$country))
  file_path <- paste0("data/formatted/", comparison, "_umi_counts_combat_seq.csv")
  write.csv(data.combat_seq, file_path)
}

create_combined_filter_with_combat_seq_files(comparison = "CFRDVsIGT", 
                                             classes = c("CFRD", "IGT"))
create_combined_filter_with_combat_seq_files(comparison = "CFRDVsNGT", 
                                             classes = c("CFRD", "NGT"))
create_combined_filter_with_combat_seq_files(comparison = "IGTVsNGT", 
                                             classes = c("IGT", "NGT"))


##########################################################################
#create one combined seurat3 output for CFRD, IGT, NGT adult samples
# since 1 file for each comparison causes issues in IGTVsNGT

#creating 2 files - one with norm and find_var_features, one without

data <- read.table("data/formatted/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                   nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
phenotype <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t") %>%
  filter(age_group == "adult") %>%
  filter(!is.na(CFRDVsIGT) | !is.na(CFRDVsNGT) | !is.na(IGTVsNGT)) %>%
  mutate(updated_sample_name = paste(condition, Sample, sep = "_"))

data <- data[, phenotype$Sample]
colnames(data) <- phenotype$updated_sample_name

phenotype <- phenotype %>%
  column_to_rownames("updated_sample_name")

rnaseq <- CreateSeuratObject(counts = data, meta.data = phenotype)
str(rnaseq)

rnaseq.list <- SplitObject(rnaseq, split.by = "country")

for (i in 1:length(rnaseq.list)) {
  rnaseq.list[[i]] <- NormalizeData(rnaseq.list[[i]], verbose = FALSE)
  rnaseq.list[[i]] <- FindVariableFeatures(rnaseq.list[[i]], selection.method = "vst", 
                                           nfeatures = 2000, verbose = FALSE)
}
rnaseq.anchors <- FindIntegrationAnchors(object.list = rnaseq.list)

#k.weight default value is 100
#choose a value less than smallest dataset size - here DK-60, AU-41
rnaseq.integrated <- IntegrateData(anchorset = rnaseq.anchors, k.weight = 41)

# Run the standard workflow for visualization and clustering
rnaseq.integrated <- ScaleData(rnaseq.integrated, verbose = FALSE)
rnaseq.integrated <- RunPCA(rnaseq.integrated, npcs = 30, verbose = FALSE)
rnaseq.integrated <- RunUMAP(rnaseq.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(rnaseq.integrated, reduction = "umap", group.by = "country")
p2 <- DimPlot(rnaseq.integrated, reduction = "umap", group.by = "condition", label = TRUE, 
              repel = TRUE) + NoLegend()
p1 + p2

data.int <- as.data.frame(rnaseq.integrated$integrated@data)
all.equal(colnames(data.int), rownames(phenotype))
colnames(data.int) <- phenotype$Sample
write.csv(data.int, "data/formatted/umi_counts_seurat3_with_norm_and_find_var_feat.csv")



rnaseq.list <- SplitObject(rnaseq, split.by = "country")
rnaseq.anchors <- FindIntegrationAnchors(object.list = rnaseq.list)

#k.weight default value is 100
#choose a value less than smallest dataset size - here DK-60, AU-41
rnaseq.integrated <- IntegrateData(anchorset = rnaseq.anchors, k.weight = 41)

# Run the standard workflow for visualization and clustering
rnaseq.integrated <- ScaleData(rnaseq.integrated, verbose = FALSE)
rnaseq.integrated <- RunPCA(rnaseq.integrated, npcs = 30, verbose = FALSE)
rnaseq.integrated <- RunUMAP(rnaseq.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(rnaseq.integrated, reduction = "umap", group.by = "country")
p2 <- DimPlot(rnaseq.integrated, reduction = "umap", group.by = "condition", label = TRUE, 
              repel = TRUE) + NoLegend()
p1 + p2

data.int_not_full <- as.data.frame(rnaseq.integrated$integrated@data)
all.equal(colnames(data.int_not_full), rownames(phenotype))
colnames(data.int_not_full) <- phenotype$Sample
write.csv(data.int_not_full, 
          "data/formatted/umi_counts_seurat3_without_norm_and_find_var_feat.csv")
##########################################################################

#batch effect removal in DK using AU

# comparison = "CFRDVsIGT"
# classes = c("IGT", "CFRD")

#filtering only 2 conditions of interest and performing batch correction with seurat seems to not work
#throws some error
#therefore, use samples from all 3 conditions
#create 1st file with filter

data <- read.table("data/formatted/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                   nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
phenotype <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t")  

output_labels <- phenotype %>%
  filter(!is.na(CFRDVsIGT) | !is.na(CFRDVsNGT) | !is.na(IGTVsNGT)) %>%
  mutate(updated_sample_name = paste(condition, Sample, sep = "_")) %>%
  column_to_rownames("updated_sample_name")

data <- data[, output_labels$Sample]
colnames(data) <- rownames(output_labels)

output_labels.au <- output_labels %>%
  filter(country == "AU")
output_labels.dk <- output_labels %>%
  filter(country == "DK")

data.au <- data[, rownames(output_labels.au)]
data.dk <- data[, rownames(output_labels.dk)]

keep <- edgeR::filterByExpr(data.au, group = output_labels.au$condition)
data.au <- data.au[keep, ]
data.dk <- data.dk[keep, ]

data.au.seurat <- CreateSeuratObject(counts = data.au, meta.data = output_labels.au)
str(data.au.seurat)
data.dk.seurat <- CreateSeuratObject(counts = data.dk, meta.data = output_labels.dk)
str(data.dk.seurat)

#ref : https://satijalab.org/seurat/reference/transferdata

# perform standard preprocessing on each object
data.au.seurat <- NormalizeData(data.au.seurat)
data.au.seurat <- FindVariableFeatures(data.au.seurat)
data.au.seurat <- ScaleData(data.au.seurat)
data.au.seurat.assay_orig <- as.data.frame(GetAssayData(data.au.seurat[['RNA']]))

data.dk.seurat <- NormalizeData(data.dk.seurat)
data.dk.seurat <- FindVariableFeatures(data.dk.seurat)
data.dk.seurat <- ScaleData(data.dk.seurat)
data.dk.seurat.assay_orig <- as.data.frame(GetAssayData(data.dk.seurat[['RNA']]))

# find anchors
anchors <- FindTransferAnchors(reference = data.au.seurat, query = data.dk.seurat, 
                               k.filter = 20)

# transfer labels
assay_modified <- TransferData(
  anchorset = anchors,
  refdata = GetAssayData(data.au.seurat[['RNA']]),
  k.weight = 20
)

data.dk <- as.data.frame(GetAssayData(assay_modified))
data.au <- as.data.frame(GetAssayData(data.au.seurat[['RNA']]))

all.equal(data.au, data.au.seurat.assay_orig)
all.equal(data.dk, data.dk.seurat.assay_orig)

data.combined <- cbind(data.au, data.dk)
output_labels.combined <- rbind(output_labels.au, output_labels.dk)
all.equal(rownames(output_labels.combined), colnames(data.combined))
#TRUE
colnames(data.combined) <- output_labels.combined$Sample
write.csv(data.combined, "data/formatted/umi_counts_filtered_seurat3_au_ref.csv")


#############################
#2nd file same process without filter
data <- read.table("data/formatted/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                   nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
phenotype <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t")  

output_labels <- phenotype %>%
  filter(!is.na(CFRDVsIGT) | !is.na(CFRDVsNGT) | !is.na(IGTVsNGT)) %>%
  mutate(updated_sample_name = paste(condition, Sample, sep = "_")) %>%
  column_to_rownames("updated_sample_name")

data <- data[, output_labels$Sample]
colnames(data) <- rownames(output_labels)

output_labels.au <- output_labels %>%
  filter(country == "AU")
output_labels.dk <- output_labels %>%
  filter(country == "DK")

data.au <- data[, rownames(output_labels.au)]
data.dk <- data[, rownames(output_labels.dk)]

data.au.seurat <- CreateSeuratObject(counts = data.au, meta.data = output_labels.au)
str(data.au.seurat)
data.dk.seurat <- CreateSeuratObject(counts = data.dk, meta.data = output_labels.dk)
str(data.dk.seurat)

#ref : https://satijalab.org/seurat/reference/transferdata

# perform standard preprocessing on each object
data.au.seurat <- NormalizeData(data.au.seurat)
data.au.seurat <- FindVariableFeatures(data.au.seurat)
data.au.seurat <- ScaleData(data.au.seurat)
data.au.seurat.assay_orig <- as.data.frame(GetAssayData(data.au.seurat[['RNA']]))

data.dk.seurat <- NormalizeData(data.dk.seurat)
data.dk.seurat <- FindVariableFeatures(data.dk.seurat)
data.dk.seurat <- ScaleData(data.dk.seurat)
data.dk.seurat.assay_orig <- as.data.frame(GetAssayData(data.dk.seurat[['RNA']]))

# find anchors
anchors <- FindTransferAnchors(reference = data.au.seurat, query = data.dk.seurat, 
                               k.filter = 20)

# transfer labels
assay_modified <- TransferData(
  anchorset = anchors,
  refdata = GetAssayData(data.au.seurat[['RNA']]),
  k.weight = 20
)

data.dk <- as.data.frame(GetAssayData(assay_modified))
data.au <- as.data.frame(GetAssayData(data.au.seurat[['RNA']]))

all.equal(data.au, data.au.seurat.assay_orig)
all.equal(data.dk, data.dk.seurat.assay_orig)

data.combined <- cbind(data.au, data.dk)
output_labels.combined <- rbind(output_labels.au, output_labels.dk)
all.equal(rownames(output_labels.combined), colnames(data.combined))
#TRUE
colnames(data.combined) <- output_labels.combined$Sample
write.csv(data.combined, "data/formatted/umi_counts_no_filter_seurat3_au_ref.csv")


##############################################

#seurat after filterByExpr

data <- read.table("data/formatted/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                   nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
phenotype <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t") %>%
  filter(age_group == "adult") %>%
  filter(!is.na(CFRDVsIGT) | !is.na(CFRDVsNGT) | !is.na(IGTVsNGT)) %>%
  mutate(updated_sample_name = paste(condition, Sample, sep = "_"))

data <- data[, phenotype$Sample]
colnames(data) <- phenotype$updated_sample_name

phenotype <- phenotype %>%
  column_to_rownames("updated_sample_name")

keep <- edgeR::filterByExpr(data, group = phenotype$condition)
data <- data[keep, ] 

rnaseq <- CreateSeuratObject(counts = data, meta.data = phenotype)
str(rnaseq)

rnaseq.list <- SplitObject(rnaseq, split.by = "country")

for (i in 1:length(rnaseq.list)) {
  rnaseq.list[[i]] <- NormalizeData(rnaseq.list[[i]], verbose = FALSE)
  rnaseq.list[[i]] <- FindVariableFeatures(rnaseq.list[[i]], selection.method = "vst", 
                                           nfeatures = 2000, verbose = FALSE)
}
rnaseq.anchors <- FindIntegrationAnchors(object.list = rnaseq.list)

#k.weight default value is 100
#choose a value less than smallest dataset size - here DK-60, AU-41
rnaseq.integrated <- IntegrateData(anchorset = rnaseq.anchors, k.weight = 41)

# Run the standard workflow for visualization and clustering
rnaseq.integrated <- ScaleData(rnaseq.integrated, verbose = FALSE)
rnaseq.integrated <- RunPCA(rnaseq.integrated, npcs = 30, verbose = FALSE)
rnaseq.integrated <- RunUMAP(rnaseq.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(rnaseq.integrated, reduction = "umap", group.by = "country")
p2 <- DimPlot(rnaseq.integrated, reduction = "umap", group.by = "condition", label = TRUE, 
              repel = TRUE) + NoLegend()
p1 + p2

data.int <- as.data.frame(rnaseq.integrated$integrated@data)
all.equal(colnames(data.int), rownames(phenotype))
colnames(data.int) <- phenotype$Sample
write.csv(data.int, "data/formatted/umi_counts_filtered_seurat3_with_norm_and_find_var_feat.csv")

####################################################################

#seurat after filterByExpr on pre-postmodulator child samples
#objective is to obtain prediction from model with adult samples on CFRDVsIGT, CFRDVsNGT, IGTVsNGT
#so same transcripts should be used

data <- read.table("data/formatted/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                   nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
phenotype <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t") %>%
  filter(age_group == "adult") %>%
  filter(!is.na(CFRDVsIGT) | !is.na(CFRDVsNGT) | !is.na(IGTVsNGT)) %>%
  mutate(updated_sample_name = paste(condition, Sample, sep = "_"))

data.train <- data[, phenotype$Sample]
colnames(data.train) <- phenotype$updated_sample_name

phenotype <- phenotype %>%
  column_to_rownames("updated_sample_name")

keep.train <- edgeR::filterByExpr(data.train, group = phenotype$condition)
data.train <- data.train[keep.train, ] 


phenotype.prepost <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t") %>%
  filter(age_group == "child" & condition == "CF_pre_post_modulator") %>%
  mutate(updated_sample_name = paste("prepost", Sample, sep = "_"))
summary(factor(phenotype.prepost$country))

data.prepost <- data[, phenotype.prepost$Sample]
colnames(data.prepost) <- phenotype.prepost$updated_sample_name

phenotype.prepost <- phenotype.prepost %>%
  column_to_rownames("updated_sample_name")

keep.prepost <- edgeR::filterByExpr(data.prepost, group = phenotype.prepost$condition)
data.prepost <- data.prepost[keep.prepost, ]   


ggvenn(list("PrePost filtered" = rownames(data.prepost),
            "Adult CFRD/IGT/NGT filtered" = rownames(data.train)),
       stroke_size = 0.1, fill_color = c("skyblue3", "yellow3"),
       set_name_size = 4,
       text_size = 3, stroke_linetype = "blank")
ggsave("prediction_pipeline/plots/venn/filtered_overlap.png")



# prepost data with transcripts from filtered data.train
data.prepost <- data[, phenotype.prepost$Sample]
colnames(data.prepost) <- rownames(phenotype.prepost)
data.prepost <- data.prepost[keep.train, ]



rnaseq <- CreateSeuratObject(counts = data.prepost, meta.data = phenotype.prepost)
str(rnaseq)

rnaseq <- NormalizeData(rnaseq, verbose = FALSE)
rnaseq <- FindVariableFeatures(rnaseq, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
rnaseq <- ScaleData(rnaseq, verbose = FALSE)

rnaseq <- RunPCA(rnaseq, npcs = 30, verbose = FALSE)
rnaseq <- RunUMAP(rnaseq, reduction = "pca", dims = 1:30)

p1 <- DimPlot(rnaseq, reduction = "umap", group.by = "pre_post_modulator")
p1

# data.prepost.pp <- as.data.frame(rnaseq$RNA@scale.data)

# all.equal(colnames(data.prepost.pp), rownames(phenotype.prepost))
# colnames(data.prepost.pp) <- phenotype.prepost$Sample
# write.csv(data.prepost.pp, 
#           "data/formatted/umi_counts_prepost_filtered_seurat3.csv")


data.train <- read.csv("data/formatted/umi_counts_filtered_seurat3_with_norm_and_find_var_feat.csv",
                       row.names = 1)
rowSums(data.train)
#not almost 0
#so save rnaseq$RNA@data for prepost

data.prepost.pp <- as.data.frame(rnaseq$RNA@data)
all.equal(colnames(data.prepost.pp), rownames(phenotype.prepost))
colnames(data.prepost.pp) <- phenotype.prepost$Sample

rowSums(data.prepost.pp)

write.csv(data.prepost.pp,
          "data/formatted/umi_counts_prepost_filtered_seurat3.csv")



##############################

#check if dk cohort samples used are postmodulators

dk <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t") %>%
  filter(country == "DK")
summary(factor(dk$pre_post_modulator))
# 0    1 NA's 
# 33   66   34 
summary(factor(dk$age_group))  
# adult 
# 133 

dk_used <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t") %>%
  filter(age_group == "adult") %>%
  filter(!is.na(CFRDVsIGT) | !is.na(CFRDVsNGT) | !is.na(IGTVsNGT)) %>%
  mutate(updated_sample_name = paste(condition, Sample, sep = "_"))
summary(factor(dk_used$pre_post_modulator))
# 0 NA's 
# 42   59 


############################################

#seurat after filterByExpr based on new phenotype

data <- read.table("data/formatted/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                   nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
phenotype <- read.table("data/formatted/phenotype_new.txt", header=TRUE, sep="\t") %>%
  filter(!is.na(CFRDVsIGT) | !is.na(CFRDVsNGT) | !is.na(IGTVsNGT)) %>%
  mutate(updated_sample_name = paste(condition, Sample, sep = "_"))

data <- data[, phenotype$Sample]
colnames(data) <- phenotype$updated_sample_name

phenotype <- phenotype %>%
  column_to_rownames("updated_sample_name")

keep <- edgeR::filterByExpr(data, group = phenotype$condition)
data <- data[keep, ] 

rnaseq <- CreateSeuratObject(counts = data, meta.data = phenotype)
str(rnaseq)

rnaseq.list <- SplitObject(rnaseq, split.by = "country")

for (i in 1:length(rnaseq.list)) {
  rnaseq.list[[i]] <- NormalizeData(rnaseq.list[[i]], verbose = FALSE)
  rnaseq.list[[i]] <- FindVariableFeatures(rnaseq.list[[i]], selection.method = "vst", 
                                           nfeatures = 2000, verbose = FALSE)
}
rnaseq.anchors <- FindIntegrationAnchors(object.list = rnaseq.list)

#k.weight default value is 100
#choose a value less than smallest dataset size - here DK-60, AU-41
rnaseq.integrated <- IntegrateData(anchorset = rnaseq.anchors, k.weight = 41)

# Run the standard workflow for visualization and clustering
rnaseq.integrated <- ScaleData(rnaseq.integrated, verbose = FALSE)
rnaseq.integrated <- RunPCA(rnaseq.integrated, npcs = 30, verbose = FALSE)
rnaseq.integrated <- RunUMAP(rnaseq.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(rnaseq.integrated, reduction = "umap", group.by = "country")
p2 <- DimPlot(rnaseq.integrated, reduction = "umap", group.by = "condition", label = TRUE, 
              repel = TRUE) + NoLegend()
p1 + p2

data.int <- as.data.frame(rnaseq.integrated$integrated@data)
all.equal(colnames(data.int), rownames(phenotype))
colnames(data.int) <- phenotype$Sample
write.csv(data.int, "data/formatted/umi_counts_filtered_seurat3_with_more_samples.csv")



############ create subsets of this

best_features <- read.csv("data/selected_features/best_features_with_is_best.csv")   %>%
  filter(is_best == 1, grepl("CF_EV_adult_filtered_seurat3_norm_find_var_none_", dataset_id)) %>%
  separate(dataset_id, into = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, "comparison"), remove = FALSE) 

#write_subset_file from pipeline_results_analyze_features.R



####################
#check number of mirnas and piRNAs
data <- read.table("data/formatted/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                   nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
transcripts <- rownames(data)
pirnas <- transcripts[grepl("piR", transcripts)]
mirnas <- transcripts[!grepl("piR", transcripts)]


####################

#check and obtain umi counts from the new set of 62 samples

data <- read_excel("../../Qiagen/qiagen_expression_output/Test experiment piRNA matrix.xlsx")

sum(is.na(data$Name))
#1

which(is.na(data$Name))
#3101

data <- read_excel("../../Qiagen/qiagen_expression_output/Test experiment piRNA matrix.xlsx") %>%
  filter(!is.na(Name)) %>%
  mutate(across(!contains("Name"), as.numeric))

###############
#the above was with spikein - but rnaseq in Ramaciotti did not use spikein
#below are the processed files without spikein

data1 <- read_excel("data/formatted/rna_new/qiagen_output/disease_NGT_control/disease_status piRNA matrix.xlsx")
data2 <- read_excel("data/formatted/rna_new/qiagen_output/disease_all_pairs/disease_status_all_comparison piRNA matrix.xlsx")
data3 <- read_excel("data/formatted/rna_new/qiagen_output/disease_random/disease_status_random piRNA matrix.xlsx")

all.equal(data1, data2)
all.equal(data1, data3)
all.equal(data2, data3)

#all true
# therefore expression count doesn't consider condition i.e. normalize it in some way
rm(data1, data2, data3)

all.equal(data, data1)
#there is mismatch b/w this and processing assuming spikein

data <- read_excel("data/formatted/rna_new/qiagen_output/disease_NGT_control/disease_status piRNA matrix.xlsx")
sum(is.na(data$Name))

which(is.na(data$Name))

data <- data %>%
  filter(!is.na(Name)) %>%
  mutate(across(!contains("Name"), as.numeric))
dim(data)
# [1] 6574   63
data_grouped <- data %>% 
  group_by(Name) %>%
  summarize(n = n()) %>%
  filter(n > 1)
#0

length(unique(data$Name))
#6574

data <- data %>%
  column_to_rownames("Name")

rsum <- rowSums(data)

data <- data %>%
  filter(rowSums(data) != 0)
dim(data)
# [1] 5359   62


#check if all samples uploaded are unique - i.e. no re-upload of same data
#the below read file is created by copy paste from qiagen output - with all 4 columns in different rows
rna_new_rough_data <- read.table("data/rna_new/all_sample_names_initial.txt", header = FALSE, sep = "\n")
(dim(rna_new_rough_data)[1]-4)/4

#used -4 to remove header lines
rna_new <- data.frame(matrix(ncol = 4, nrow = (dim(rna_new_rough_data)[1]-4)/4))
colnames(rna_new) <- rna_new_rough_data[c(1:4), 1]

i <- 5
row_num <- floor((i-1) / 4)
while(row_num < nrow(rna_new)){
  # print(i)
  row_num <- floor((i-1) / 4)
  # print(row_num)
  rna_new[row_num, ] <- c(rna_new_rough_data[c(i:(i+3)), 1])
  i <- i + 4
}

length(unique(rna_new$`Sample name`))
#334

#hence verified that all uploaded samples are different


############
#formatting quantified results from Qiagen RNASeqPortal for all samples : 334 (as of Oct 27, 2023)
data <- read_excel("data/formatted/rna_all/Quantify_all_334 piRNA matrix.xlsx")
sum(is.na(data$Name))

which(is.na(data$Name))

data <- data %>%
  filter(!is.na(Name)) %>%
  mutate(across(!contains("Name"), as.numeric))
dim(data)
# [1] 11449   335

data_grouped <- data %>% 
  group_by(Name) %>%
  summarize(n = n()) %>%
  filter(n > 1)
#0

length(unique(data$Name))
# [1] 11449
#so all unique

data <- data %>%
  column_to_rownames("Name")

rsum <- rowSums(data)

data <- data %>%
  filter(rowSums(data) != 0)
dim(data)
# [1] 10654   334

write.csv(data, "data/formatted/rna_all/umi_counts.csv")

#perform some qc to check :
#for each of the transcripts, perc of zeroes of that transcript across all samples
zero_entries <- as.data.frame(data == 0)
perc_zero <- data.frame(perc = 100 * rowSums(zero_entries) / ncol(zero_entries))

summary(perc_zero$perc)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   97.31   99.40   90.57   99.70   99.70

ggplot(perc_zero, aes(x = "", y = perc)) +
  geom_boxplot() +
  labs(x = "",
       y = "% of zeroes across all samples",
       title = "Boxplot of % of zeroes in each transcript across all samples")
ggsave("data/formatted/rna_all/qc/perc_zero_boxplot.png")

ggplot(perc_zero, aes(x = perc)) +
  geom_histogram(binwidth = 1) +
  labs(x = "% of zeroes across all samples",
       y = "number of transcripts",
       title = "Frequency of % of zeroes in each transcript across all samples",
       caption = paste("Total number of transcripts =", nrow(perc_zero))) +
  scale_x_continuous(breaks = seq(0, 100, by = 5)) +
  scale_y_continuous(n.breaks = 10)
ggsave("data/formatted/rna_all/qc/perc_zero_barplot.png")

ggplot(perc_zero, aes(x = perc)) +
  stat_ecdf(geom = "step") +
  labs(x = "% of zeroes across all samples",
       y = "Cumulative proportion of transcripts",
       title = "Cumulative distribution of % of zeroes in each transcript across all samples",
       caption = paste("Total number of transcripts =", nrow(perc_zero))) +
  scale_x_continuous(breaks = seq(0, 100, by = 5)) +
  scale_y_continuous(n.breaks = 10)
ggsave("data/formatted/rna_all/qc/perc_zero_cumulative.png")

#from the cumulative dist plot, about 0.15 transcripts have < 90% zeros across all samples

#verifying
sum(perc_zero$perc < 90) / nrow(perc_zero) 
# [1] 0.1508354

filt_data <- data[perc_zero$perc < 90, ]
nrow(filt_data)
#1607

# about 0.2 transcripts have < 95% zeros across all samples
filt_data <- data[perc_zero$perc < 95, ]
nrow(filt_data)
#2071

#using median
filt_data <- data[perc_zero$perc < 99.4, ]
nrow(filt_data)
#4989

#using 3rd quartile
filt_data <- data[perc_zero$perc < 99.7, ]
nrow(filt_data)
#6615


#use which of the above ?
# maybe compare using expression heatmap

meta_data <- read.csv("data/formatted/tra_metadata_all_2023Oct.csv") %>%
  dplyr::select(c(sample_long_name, condition, cohort, pre_post_modulator, batch_name)) %>%
  mutate(modulator_status = ifelse(is.na(pre_post_modulator) | pre_post_modulator == 0, "pre", "post")) %>%
  dplyr::select(-c(pre_post_modulator)) %>%
  mutate(condition = factor(condition),
         cohort = factor(cohort),
         modulator_status = factor(modulator_status),
         batch_name = factor(batch_name))

create_expression_heatmap <- function(expr_data, meta_data, file_name){
  expr_data <- expr_data[, meta_data$sample_long_name]
  all.equal(meta_data$sample_long_name, colnames(expr_data))
  
  data_to_plot <- as.matrix(log2(expr_data + 2^-10))
  
  ht <- Heatmap(data_to_plot, name = "Log2 Transcriptomics expression",
          col = viridis(n = 10, option = "magma"),
          show_column_names = FALSE,
          show_column_dend = FALSE,
          show_row_names = FALSE,
          show_row_dend = FALSE,
          row_title = paste0("Transcripts (", nrow(data_to_plot), ")"),
          column_title = "Samples",
          bottom_annotation = HeatmapAnnotation(
            "Condition" = meta_data$condition,
            "Modulator Status" = meta_data$modulator_status,
            "Cohort" = meta_data$cohort,
            "Batch name" = meta_data$batch_name,
            col = list("Condition" = c("CFRD" = "red",
                                       "IGT" = "orange",
                                       "NGT" = "yellow",
                                       "HC" = "green",
                                       "UNKNOWN TO PREDICT" = "gray60"),
                       "Modulator Status" = c("pre" = "aquamarine",
                                              "post" = "seagreen"),
                       "Cohort" = c("CPH" = "green",
                                    "RPA_NSW" = "magenta",
                                    "SCH_NSW" = "plum",
                                    "UNSW" = "red"),
                       "Batch name" = c("initial" = "coral",
                                        "new" = "cyan"))
          ))
  plot_dir_path <- "data/formatted/rna_all/qc/heatmap/"
  if(!dir.exists(plot_dir_path)){
    dir.create(plot_dir_path, recursive = TRUE)
  }
  file_path <- paste0(plot_dir_path, file_name)
  png(file_path, units = "cm", width = 20, height = 25, res = 1200)  
  draw(ht)    
  dev.off()
}

create_expression_heatmap(data, meta_data, "1_zero_filtered.png")

filt_data <- data[perc_zero$perc < 99.4, ]
nrow(filt_data)
create_expression_heatmap(filt_data, meta_data, "2_median_filtered.png")

filt_data <- data[perc_zero$perc < 95, ]
create_expression_heatmap(filt_data, meta_data, "3_95_filtered.png")

filt_data <- data[perc_zero$perc < 90, ]
create_expression_heatmap(filt_data, meta_data, "4_90_filtered.png")

nrow(filt_data)
#1607

write.csv(filt_data, "data/formatted/rna_all/umi_counts_filter90.csv")




############
#checking quantified results from Qiagen RNASeqPortal for redemultiplexed new (62) samples (provided by Ramaciotti on Dec 2023)
data.redemux <- read_excel("data/formatted/rna_new/redemultiplexed/QC experiment piRNA matrix.xlsx")

sum(is.na(data.redemux$Name))
#0

data.redemux <- data.redemux %>%
  mutate(across(!contains("Name"), as.numeric))
dim(data.redemux)
# [1] 6574   63

data.redemux_grouped <- data.redemux %>% 
  group_by(Name) %>%
  summarize(n = n()) %>%
  filter(n > 1)
#0

length(unique(data.redemux$Name))
# [1] 6574
#so all unique


data.redemux <- data.redemux %>%
  column_to_rownames("Name")

rsum <- rowSums(data.redemux)

data.redemux <- data.redemux %>%
  filter(rowSums(data.redemux) != 0)
dim(data.redemux)
# [1] 5359   62



data <- read_excel("data/formatted/rna_new/qiagen_output/disease_NGT_control/disease_status piRNA matrix.xlsx")
sum(is.na(data$Name))

which(is.na(data$Name))

data <- data %>%
  filter(!is.na(Name)) %>%
  mutate(across(!contains("Name"), as.numeric))
dim(data)
# [1] 6574   63
data_grouped <- data %>% 
  group_by(Name) %>%
  summarize(n = n()) %>%
  filter(n > 1)
#0

length(unique(data$Name))
#6574

data <- data %>%
  column_to_rownames("Name")
rsum <- rowSums(data)
data <- data %>%
  filter(rowSums(data) != 0)
dim(data)
# [1] 5359   62

data <- data %>%
  rownames_to_column("Name") %>%
  arrange(Name) %>%
  column_to_rownames("Name")

data.redemux <- data.redemux %>%
  rownames_to_column("Name") %>%
  arrange(Name) %>%
  column_to_rownames("Name")

all.equal(rownames(data), rownames(data.redemux))

meta_data <- read.csv("data/formatted/rna_new/meta_data_minimal.csv")

data <- data %>%
  dplyr::select(meta_data$Sample)
data.redemux <- data.redemux %>%
  dplyr::select(meta_data$Sample)

all.equal(data, data.redemux)
# [1] TRUE

#redemultiplexed data same as original 