library(tidyverse)
library(Seurat)
library(ComplexHeatmap)
library(viridis)
library(ggplot2)


base_dir <- "~/UNSW/VafaeeLab/CysticFibrosisGroup/ExoCF/CFRD_EV_biomarker/"
setwd(base_dir)


# heatmap of samples vs normalized gene expression
# Fig 1 J
num_variable_features = 100
variable_features_only = FALSE
perform_filter = TRUE
plot_file_dir = "plots/expression_heatmap"
plot_file_name = "var100_fil_seurat.png"
norm = "seurat_norm"

generate_expression_heatmap <- function(plot_file_name, norm, perform_filter = TRUE,
                                        variable_features_only = FALSE,
                                        num_variable_features = 100,
                                        plot_file_dir = "plots/expression_heatmap"){
  data <- read.table("data/formatted/umi_counts.csv", header=TRUE, sep=",", row.names=1, skip=0,
                     nrows=-1, comment.char="", fill=TRUE, na.strings = "NA")
  phenotype <- read.table("data/formatted/phenotype.txt", header=TRUE, sep="\t")
  phenotype <- phenotype %>%
    mutate(condition2 = case_when((!is.na(pre_post_modulator) & pre_post_modulator == 1) ~ "CFPostModulator",
                                  condition %in% c("CFRD", "IGT", "NGT", "HC") ~ condition,
                                  TRUE ~ "otherCF"),
           .after = "condition"
    )
  data <- data[, phenotype$Sample]
  rownames(data) <- gsub(pattern = '_', replacement = '-', rownames(data))
                                  
  if(perform_filter){
    keep <- edgeR::filterByExpr(data, group = phenotype$condition2)
    data <- data[keep, ]    
  }

  if(variable_features_only){
    group_info <- phenotype %>% select(Sample, condition2)
    data_joined <- group_info %>%
      inner_join(as.data.frame(t(data)) %>% 
                   rownames_to_column('Sample'))
    data_joined.mean <- data_joined %>%
      group_by(condition2) %>%
      summarise(across(where(is.numeric), mean))
    data_joined.var <- data_joined.mean %>%
      summarise(across(where(is.numeric), var))
    data_joined.var <- as.data.frame(t(data_joined.var))
    colnames(data_joined.var) <- 'val'
    data_joined.var <- data_joined.var %>%
      arrange(desc(val))
    transcripts_of_interest <- rownames(data_joined.var)[1:num_variable_features]    
  }

  if(norm == "logCPM"){
    data <- edgeR::cpm(data, log = TRUE)  
  } else if(norm == "seurat_norm"){
    phenotype <- phenotype %>%
      mutate(Sample = paste(condition2, Sample, sep = "_")) %>%
      column_to_rownames("Sample")
    colnames(data) <- rownames(phenotype)
    
    data.seurat <- CreateSeuratObject(counts = data, meta.data = phenotype)
    str(data.seurat)
    data.seurat <- NormalizeData(data.seurat)
    data.seurat <- FindVariableFeatures(data.seurat, selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
    data.seurat <- ScaleData(data.seurat, verbose = FALSE)
    data.seurat.norm <- as.data.frame(GetAssayData(data.seurat[['RNA']]))
    data <- data.seurat.norm
  }
  
  if(variable_features_only){
    data <- data[transcripts_of_interest, ]
  }
  data <- data.matrix(data)
  
  if(!dir.exists(plot_file_dir)){
    dir.create(plot_file_dir, recursive = TRUE)
  }
  plot_file_path = paste(plot_file_dir, plot_file_name, sep = "/")
  png(plot_file_path, units = "cm", width = 20, height = 15, res = 1200)
  
  anno_col = list("Condition" = c("CFRD" = "red",
                                  "IGT" = "orange",
                                  "NGT" = "yellow",
                                  "otherCF" = "magenta",
                                  "CFPostModulator" = "blue",
                                  "HC" = "skyblue"
                                  ),
                  "Country" = c("AU" = "green",
                                "DK" = "navyblue"))
  
  phenotype <- phenotype %>%
    mutate(country = factor(country, levels = c("AU", "DK")),
           condition2 = factor(condition2, levels = c("CFRD", "IGT", "NGT", 
                                                      "otherCF", "CFPostModulator", "HC")))
  transcript_count <- dim(data)[1]
  ht <- Heatmap(data, show_row_names = FALSE, show_column_names = FALSE,
          row_title = paste0("Transcripts (", transcript_count, ")"), column_title = "Samples",
          col = viridis(10), name = norm,
          show_row_dend = FALSE,
          show_column_dend = FALSE, cluster_columns = FALSE,
          column_split = phenotype$condition2,
          bottom_annotation = HeatmapAnnotation(
            "Condition" = phenotype$condition2, 
            "Country" = phenotype$country,
            col = anno_col
          ))
  draw(ht)
  dev.off()
}

generate_expression_heatmap(plot_file_name = "fil_seurat.png", norm = "seurat_norm")
generate_expression_heatmap(plot_file_name = "nofil_seurat.png", norm = "seurat_norm", perform_filter = FALSE)
generate_expression_heatmap(plot_file_name = "fil_logcpm.png", norm = "logCPM")
generate_expression_heatmap(plot_file_name = "nofil_logcpm.png", norm = "logCPM", perform_filter = FALSE)
generate_expression_heatmap(plot_file_name = "fil_seurat_var100.png", norm = "seurat_norm", 
                            variable_features_only = TRUE)
generate_expression_heatmap(plot_file_name = "fil_logcpm_var100.png", norm = "logCPM",
                            variable_features_only = TRUE)


# pheno_sub <- phenotype %>%
#   filter(condition %in% c("CFRD", "IGT", "NGT")) %>%
#   filter(!is.na(pre_post_modulator) & pre_post_modulator == 1) %>%
#   select(Sample, country, age_group, condition, pre_post_modulator)
# 
# phenotype %>%
#   filter(condition == "HC") %>%
#   select(Sample, country, age_group, condition, pre_post_modulator)
