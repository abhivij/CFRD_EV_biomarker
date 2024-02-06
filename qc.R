#to check the expression levels of proteins and transcripts and samplewise expression 

library(tidyverse)
library(ggrepel)
library(ggvenn)
library(sva)
library(ggplot2)



# comparison <- "CFRDVsIGT"
# classes <- c("IGT", "CFRD")
# omics_type <- "transcriptomics"
# perform_filter = TRUE
# norm <- "log_cpm"
# batch_effect_correction <- "combat"
# data_file_path <- "data/formatted/rna_all/umi_counts_filter90.csv"
# phenotype_file_path <- "data/formatted/tra_phenotype_2024Jan.txt" 
# plot_type <- "per_sample"



create_per_feature_or_sample_plot <- function(comparison, classes,
                                              omics_type, norm,
                                              perform_filter = FALSE,
                                              batch_effect_correction = "none",
                                              plot_type = "per_feature",
                                              plot_dir_path,
                                              data_file_path = NA,
                                              phenotype_file_path = NA,
                                              plot_file_name_prefix = NA){
  if(omics_type == "proteomics"){
    if(is.na(data_file_path)){
      data_file_path <- "data/proteomics/data_333samples_imputed_mf.csv"
    }
    if(is.na(phenotype_file_path)){
      phenotype_file_path <- "data/formatted/prot_phenotype_333_2024Jan.txt"  
    }
  } else if(omics_type == "transcriptomics"){
    
    if(is.na(data_file_path)){
      data_file_path <- "data/formatted/rna_all/umi_counts_filter90.csv"
    }
    if(is.na(phenotype_file_path)){
      phenotype_file_path <- "data/formatted/tra_phenotype_2024Jan.txt"  
    }
  }
  
  data <- read.csv(data_file_path, row.names = 1)
  phenotype <- read.table(phenotype_file_path, header=TRUE, sep="\t")
  
  output_labels <- phenotype %>%
    dplyr::rename("Label" = comparison) %>%
    dplyr::filter(Label %in% classes) %>%
    dplyr::select(Sample, Label, country, cohort, age_group)

  #currently data format : (transcripts x samples)
  
  data <- data[, output_labels$Sample]
  title <- paste0(norm, " BEC-", batch_effect_correction)
  
  if(perform_filter){
    keep <- edgeR::filterByExpr(data, group = output_labels$Label)
    data <- data[keep, ]
  }
  
  if(norm == "quantile"){
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
    
  } else if(norm == "log_cpm"){
    #calculating norm log cpm
    data <- edgeR::cpm(data, log=TRUE)
  } else if(norm == "none"){
    print("just log")
    data <- log2(data + 2^-10)
  }
  
  data <- as.data.frame(t(as.matrix(data)))
  
  output_labels <- output_labels %>%
    mutate(Label = factor(Label), country = factor(country))
  
  group_counts <- output_labels %>%
    dplyr::mutate(Label = paste(Label, country, sep = "_")) %>%
    group_by(Label) %>%
    summarise(n = n())
  
  group_counts_text <- paste(apply(group_counts, MARGIN = 1, FUN = function(x){paste(x[1], x[2], sep = ":")}),
                             collapse = "  ")
  
  all.equal(rownames(data), output_labels$Sample)
  
  if(batch_effect_correction == "combat"){
    data <- as.data.frame(t(as.matrix(data)))
    data.combat = ComBat(dat=data, batch=output_labels$country)
    data.combat <- as.data.frame(t(as.matrix(data.combat)))
    data <- data.combat
  }
  
  if(plot_type == "per_feature"){
    
    data_to_plot <- data %>%
      rownames_to_column(var = "Sample") %>%
      pivot_longer(cols = !Sample, names_to = "features", values_to = "expr") %>%
      inner_join(output_labels %>%
                   mutate(Label = paste(Label, country, sep = "_")) %>%
                   dplyr::select(c(Sample, Label))
      )
    
    data_to_plot <- data_to_plot %>%
      group_by(features, Label) %>%
      summarize(med_expr = median(expr), 
                q25up = quantile(expr, probs = 0.75), 
                q25down = quantile(expr, probs = 0.25)) %>%
      ungroup()
    
    median_expr_data <- data_to_plot %>%
      group_by(features) %>%
      summarize(med_med_expr = median(med_expr)) %>%
      arrange(desc(med_med_expr))
    
    data_to_plot <- data_to_plot %>%
      pivot_longer(cols = c(med_expr, q25up, q25down), names_to = "expr_type", values_to = "expr") %>%
      mutate(features = factor(features, levels = median_expr_data$features))
    
    #for simplicity decided to just show median expression
    data_to_plot <- data_to_plot %>%
      dplyr::filter(expr_type == "med_expr")
    
    data_to_plot <- data_to_plot %>%
      inner_join(group_counts) %>%
      mutate(Label = paste0(Label, ' (', n, ')')) %>%
      dplyr::select(-c(n)) %>%
      mutate(Label = factor(Label))
    
    y_lab <- "Expression"
    if(omics_type == "transcriptomics"){
      x_lab <- "transcripts"
    } else{
      x_lab <- "proteins"
      
      # data_to_plot <- data_to_plot %>%
      #   left_join(all_protein_names %>% dplyr::select(c(from_id, primary_gene_id)), 
      #             by = c("biomarker" = "from_id")) %>%
      #   dplyr::select(-c(biomarker)) %>%
      #   dplyr::rename(c("biomarker" = "primary_gene_id"))
    }
    x_lab <- paste0(x_lab, "(", length(median_expr_data$features), ")")
    
    # biomarker_agg <- data_to_plot %>%
    #   group_by(biomarker) %>%
    #   summarize(med_expr = median(norm_expr)) %>%
    #   arrange(desc(med_expr))
    # 
    # # data_to_plot <- data_to_plot %>%
    # #   mutate(Label = factor(Label, levels = rev(classes)),
    # #          biomarker = factor(biomarker, biomarker_agg$biomarker)) 
    # data_to_plot <- data_to_plot %>%
    #   mutate(biomarker = factor(biomarker, biomarker_agg$biomarker))
    # # print("here")
    
    # ggplot(data_to_plot) +
    #   geom_point(aes(x = features, y = expr, fill = Label, shape = expr_type), stroke = 0.1) +
    #   theme(axis.text.x = element_blank()) +
    #   scale_shape_manual(name = "Expression type", values = c(21, 24, 25)) +
    #   guides(fill = guide_legend(override.aes = list(shape = 21, colour = NA))) +
    #   ggtitle(paste(sub("Vs", " Vs ", comparison), x_lab)) +
    #   xlab(x_lab) +
    #   ylab(y_lab)
    
    ggplot(data_to_plot) +
      geom_point(aes(x = features, y = expr, color = Label)) +
      theme(axis.text.x = element_blank()) +
      ggtitle(paste(sub("Vs", " Vs ", comparison), x_lab, title)) +
      xlab(x_lab) +
      ylab(y_lab)    
  }  else if(plot_type == "per_sample"){
    
    data_to_plot <- data %>%
      rownames_to_column(var = "Sample") %>%
      pivot_longer(cols = !Sample, names_to = "features", values_to = "expr") %>%
      inner_join(output_labels %>%
                   mutate(Label = paste(Label, country, sep = "_")) %>%
                   dplyr::select(c(Sample, Label))
      )
    data_to_plot <- data_to_plot %>%
      mutate(Label = factor(Label))
    
    median_expr_data <- data_to_plot %>%
      group_by(Sample, Label) %>%
      summarize(med_expr = median(expr)) %>%
      arrange(Label, desc(med_expr))
    
    data_to_plot <- data_to_plot %>%
      mutate(Sample = factor(Sample, levels = median_expr_data$Sample))
    
    data_to_plot <- data_to_plot %>%
      inner_join(group_counts) %>%
      mutate(Label = paste0(Label, ' (', n, ')')) %>%
      dplyr::select(-c(n)) %>%
      mutate(Label = factor(Label))
    
    y_lab <- "Expression"
    x_lab <- "Sample"
    # x_lab <- paste0(x_lab, "(", length(data_to_plot$Sample), ")")
    
    ggplot(data_to_plot) +
      geom_boxplot(aes(x = Sample, y = expr, color = Label)) +
      theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.99, size = rel(0.75))) +
      ggtitle(paste(sub("Vs", " Vs ", comparison), x_lab, title)) +
      xlab("Sample") +
      ylab(y_lab)    
  }
  
  
  if(!dir.exists(plot_dir_path)){
    dir.create(plot_dir_path, recursive = TRUE)
  }
  file_name <- paste(comparison, x_lab, title, ".png", sep = "_")
  if(!is.na(plot_file_name_prefix)){
    file_name <- paste0(plot_file_name_prefix, file_name)
  }
  file_path <- paste(plot_dir_path, file_name, sep = "/")
  ggplot2::ggsave(file_path, units = "cm", width = 30)
  
}



