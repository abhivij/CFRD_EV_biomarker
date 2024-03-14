library(tidyverse)
library(viridis)
library(umap)
library(ggrepel)
library(VennDiagram)
library(SummarizedExperiment)

library(ggplot2)
library(tidyverse)
library(gridExtra)
library(RColorBrewer)
library(ggfortify)
library(DEP)
library(factoextra)
library(UpSetR)
library(EnhancedVolcano)
library(magrittr)
library(data.table)
library(SummarizedExperiment)
library(reshape2)
library(circlize)
library(ComplexHeatmap)
library(rstatix)
library(edgeR)
library(ggrepel)
library(gplots)
library(ggvenn)
library(ggpubr)
library(xlsx)
options(scipen=999) #removes scientific notation


create_volcano_plot <- function(result, title, file_name, dir_path = "",
                                p_val_cutoff = 0.05, logFC_cutoff = 5,
                                p_val_column = "pVal",
                                k = 5) {
  # if(use_p_val){
  #   p_val_column <- "pVal"
  # } else{
  #   p_val_column <- "adjPVal"
  # }
  
  result$diffexpr <- "NO"
  result$diffexpr[result$logFC > logFC_cutoff & result[, p_val_column] < p_val_cutoff] <- "UP"
  result$diffexpr[result$logFC < -logFC_cutoff & result[, p_val_column] < p_val_cutoff] <- "DOWN"
  
  
  rownames(result) <- result$Molecule
  up_reg_result <- result %>%
    filter(logFC > 0 & .[[p_val_column]] < p_val_cutoff)
  if(dim(up_reg_result)[1] == 0){
    up_reg_result <- result %>%
      filter(logFC > 0)    
  }
  up_reg_indices <- order(up_reg_result$logFC,decreasing = TRUE)[1:k]
  up_reg_indices <- rownames(up_reg_result[up_reg_indices,])
  
  down_reg_result <- result %>%
    filter(logFC < 0 & .[[p_val_column]] < p_val_cutoff)
  if(dim(down_reg_result)[1] == 0){
    down_reg_result <- result %>%
      filter(logFC < 0) 
  }  
  down_reg_indices <- order(-down_reg_result$logFC, decreasing = TRUE)[1:k] 
  down_reg_indices <- rownames(down_reg_result[down_reg_indices,])
  
  result$de_mol <- NA
  result[up_reg_indices, "de_mol"] <- result[up_reg_indices, "Molecule"]
  result[down_reg_indices, "de_mol"] <- result[down_reg_indices, "Molecule"]
  
  
  colours <- c("red", "blue", "grey")
  names(colours) <- c("UP", "DOWN", "NO")
  
  # https://stackoverflow.com/questions/20005233/how-to-use-simultaneously-superscript-and-variable-in-a-axis-label-with-ggplot2
  
  #just superscript or subscript can be added in axis label using expression
  #using value from variable along with it requires bquote
  
  vol_plot <- ggplot(data=result, aes(x=logFC, y=-log10(!!as.symbol(p_val_column)), col=diffexpr, label=de_mol)) +
    geom_point() + 
    geom_text_repel(max.overlaps = 20, colour = "black", size = 2) +
    geom_vline(xintercept=c(-logFC_cutoff, logFC_cutoff), col="green", linetype = "dashed") +
    geom_hline(yintercept=-log10(p_val_cutoff), col="green", linetype = "dashed") +
    scale_colour_manual(values = colours) +
    labs(title = title, colour = "Diff Expr",
         x = expression(paste(log[2], "FC")),
         y = bquote(-log[10] ~ .(p_val_column)))
  vol_plot
  ggsave(paste(dir_path, file_name, sep = "/"), vol_plot) 
}



#renames condition_column_name to condition and adds a column named replicate having unique values 
#within each of the condition 
#eg: cond1 : replicate from 1 to 5, cond2 : replicate from 1 to 10, cond3 : replicate from 1 to 6
add_replicate_column <- function(sample_info, condition_column_name){
  sample_info_temp <- sample_info %>%
    dplyr::rename("condition" = condition_column_name) %>%
    group_by(condition) %>%
    mutate(replicate = row_number()) %>%
    ungroup()
  return (sample_info_temp)
}

#create histogram of protein count
plot_protein_count <- function(se, title, plot_file_name, 
                               is_condition_factor = FALSE,
                               condition_levels = c()){
  
  data <- get_df_long(se) %>%
    select(label, condition, replicate, name, intensity) %>%
    mutate(label = sub("LFQ.intensity.", "", label, fixed = TRUE)) %>%
    tidyr::separate(label, into = c(NA, NA, NA, "label"), sep = "-") %>%
    filter(!is.na(intensity))
  
  grouped_info <- data %>%
    group_by(label, condition, replicate) %>%
    summarize(n = n())
  
  if(is_condition_factor){
    grouped_info <- grouped_info %>%
      mutate(condition = factor(condition, levels = condition_levels))
  }
  
  grouped_info <- grouped_info %>%
    arrange(condition, replicate)
  
  grouped_info <- grouped_info %>%
    mutate(label = factor(label, levels = grouped_info$label))
  
  ggplot(grouped_info) +
    geom_bar(aes(x = label, y = n, fill = condition), stat = "identity") +
    xlab("Sample") +
    ylab("Protein Count") +
    ggtitle(title) +
    scale_y_continuous(breaks = seq(0, max(grouped_info$n), by = 100)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  ggsave(plot_file_name)
}

make_se <- function(proteins_unique, columns, expdesign) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(is.data.frame(proteins_unique),
                          is.integer(columns),
                          is.data.frame(expdesign))
  
  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(proteins_unique))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(proteins_unique)),
         "'.\nRun make_unique() to obtain the required columns",
         call. = FALSE)
  }
  if(any(!c("label", "condition", "replicate") %in% colnames(expdesign))) {
    stop("'label', 'condition' and/or 'replicate' columns",
         "are not present in the experimental design",
         call. = FALSE)
  }
  if(any(!apply(proteins_unique[, columns], 2, is.numeric))) {
    stop("specified 'columns' should be numeric",
         "\nRun make_se_parse() with the appropriate columns as argument",
         call. = FALSE)
  }
  
  # If input is a tibble, convert to data.frame
  if(tibble::is_tibble(proteins_unique))
    proteins_unique <- as.data.frame(proteins_unique)
  if(tibble::is_tibble(expdesign))
    expdesign <- as.data.frame(expdesign)
  
  # Select the assay data
  rownames(proteins_unique) <- proteins_unique$name
  raw <- proteins_unique[, columns]
  raw[raw == 0] <- NA
  raw <- log2(raw)
  
  # Generate the colData from the experimental design
  # and match these with the assay data
  expdesign <- mutate(expdesign, condition = make.names(condition)) %>%
    unite(ID, condition, replicate, remove = FALSE)
  rownames(expdesign) <- expdesign$ID
  
  matched <- match(make.names(expdesign$label),
                   make.names(colnames(raw)))
  if(any(is.na(matched))) {
    stop("None of the labels in the experimental design match ",
         "with column names in 'proteins_unique'",
         "\nRun make_se() with the correct labels in the experimental design",
         "and/or correct columns specification")
  }
  
  colnames(raw)[matched] <- expdesign$ID
  raw <- raw[, !is.na(colnames(raw))][rownames(expdesign)]
  
  # Select the rowData
  row_data <- proteins_unique[, -columns]
  rownames(row_data) <- row_data$name
  
  # Generate the SummarizedExperiment object
  se <- SummarizedExperiment(assays = as.matrix(raw),
                             colData = expdesign,
                             rowData = row_data)
  return(se)
} #Slightly edited from DEP as it wasn't reading names correctly


# sinfo <- read.table(file= "data/sample_info_cf.txt", sep="\t", header=T)
# # sinfo <- sinfo %>%
# #   filter(oxygen_level == "Normoxic")
# cond1 <- "NN"
# cond2 <- "NH"
# file_path <- "data/de_proteins/ALI/CF/NN_NH/"
# 
# 
# diffex(data_unique, sinfo, cond1 = "NN", cond2 = "NH", file_path = "data/de_proteins/ALI/HCplusCF/NN_NH/")

# 
# imputed_se <- imputed
# cond1 = "PreModulator_CFRD"
# cond2 = "PreModulator_IGT"
# file_path = "data/de_results_2024/proteomics/DEP_regular/p/"
# protein_names_file_path = "data/proteomics/protein_names.csv"
# fc = 1.5
# k = 10
# x_lim = NA
# y_lim = NA
# fc = 1.2
# pval_cutoff = 0.05
# use_adj_pval = FALSE

diffex <- function(imputed_se, cond1, cond2, file_path,
                   protein_names_file_path = "data/protein_names.csv",
                   k = 10,
                   x_lim = NA,
                   y_lim = NA,
                   fc = 1.2,
                   pval_cutoff = 0.05,
                   use_adj_pval = FALSE) {
  
  #extract the directory path from file_path and create the directory if required
  comp = strsplit(file_path, split = "/")[[1]]
  dir_path = paste0(comp[1:(length(comp)-1)], collapse = "/")
  if(!dir.exists(dir_path)){
    dir.create(dir_path, recursive = TRUE)
  }
  
  lfc_threshold <- log2(fc)
  diff <- test_diff(imputed_se, type = "manual", test = paste(cond1, cond2, sep = "_vs_"))
  dep <- add_rejections(diff, lfc = lfc_threshold)
  data_results <- get_results(dep)
  data_results <- data_results[, c(1,7,3,4)]
  names(data_results) <- c("Protein", "logFC", "Pvalue", "P.adj")
  
  if(use_adj_pval){
    column_name <- "P.adj"
  } else{
    column_name <- "Pvalue"
  }
  
  #arrange by significance needs to be done so as to show 
  # the labels in the order - downregulated, not significant, upregulated in the volcano plot 
  data_results <- data_results %>%
    dplyr::rename(c("pvalcolumn" = column_name)) %>%
    mutate(significance = case_when(logFC <= -lfc_threshold & pvalcolumn <= pval_cutoff ~ 'Downregulated',
                                    logFC >= lfc_threshold & pvalcolumn <= pval_cutoff ~ 'Upregulated',
                                    TRUE ~ 'Not significant')) %>% 
    mutate(colour = case_when(significance == 'Downregulated' ~ '#E8495C',
                              significance == 'Upregulated' ~ '#38ACE2',
                              TRUE ~ 'grey')) %>%
    arrange(significance)
  
  if(use_adj_pval){
    data_results <- data_results %>%
      dplyr::rename(c("P.adj" = "pvalcolumn")) 
  } else{
    data_results <- data_results %>%
      dplyr::rename(c("Pvalue" = "pvalcolumn")) 
  }
    
  
  #include protein names
  protein_names <- read.csv(protein_names_file_path) %>%
    select(name, protein_name)
  data_results <- data_results %>%
    left_join(protein_names, by = c("Protein" = "name"))
  
  upreg <- data_results %>%
    filter(significance == 'Upregulated') %>%
    arrange(desc(logFC))
  downreg <- data_results %>%
    filter(significance == "Downregulated") %>%
    arrange(logFC)
  
  if(nrow(upreg) >= k){
    top_proteins <- upreg[1:k, "Protein"]
  } else{
    top_proteins <- upreg[1:nrow(upreg), "Protein"]
  }
  if(nrow(downreg) >= k){
    top_proteins <- c(top_proteins, downreg[1:k, "Protein"])
  } else{
    top_proteins <- c(top_proteins, downreg[1:nrow(downreg), "Protein"])
  }
  
  data_results <- data_results %>%
    mutate(label = case_when(Protein %in% top_proteins ~ Protein,
                             TRUE ~ NA_character_))
  
  # data_results <- data_results %>%
  #   filter(!is.na(label)) %>%
  #   arrange(desc(logFC))
  
  sig <- rbind(upreg, downreg) %>%
    select(-c(significance, colour)) %>%
    arrange(desc(logFC))
  # print(sig)
  # print(paste("Upregulated:", nrow(sig[sig$logFC > 0,]), sep=" "))
  # print(paste("Downregulated:", nrow(sig[sig$logFC < 0,]), sep=" "))
  
  keyvals <- data_results$colour
  names(keyvals) <- data_results$significance
  
  upreg_count <- sum(names(keyvals) == 'Upregulated')
  downreg_count <- sum(names(keyvals) == 'Downregulated')
  ns_count <- sum(names(keyvals) == 'Not significant')
  
  names(keyvals) <- gsub('Upregulated', paste0('Upregulated (', upreg_count, ')'), names(keyvals))
  names(keyvals) <- gsub('Downregulated', paste0('Downregulated (', downreg_count, ')'), names(keyvals))
  names(keyvals) <- gsub('Not significant', paste0('Not significant (', ns_count, ')'), names(keyvals))
  
  if(is.na(x_lim)){
    x_lim <- c(min(data_results$logFC)-0.5,
               max(data_results$logFC)+0.5)
  }
  if(is.na(y_lim)){
    y_lim <- c(0, max(-log10(data_results$Pvalue))+0.5)
  }
  volcanoplot <- EnhancedVolcano(data_results,
                                 lab = data_results$label,
                                 x = 'logFC',
                                 y = column_name,
                                 xlim = x_lim,
                                 ylim = y_lim,
                                 pCutoff = 0.05,
                                 FCcutoff = lfc_threshold,
                                 colCustom = keyvals,
                                 colAlpha = 1,
                                 title = paste(cond1,"vs",cond2, sep=" "),
                                 legendPosition = 'bottom',
                                 subtitle = '',
                                 labSize = 3,
                                 drawConnectors = T,
                                 max.overlaps = Inf,
                                 arrowheads = F)
  
  plot(volcanoplot)
  
  ggsave(paste0(file_path, cond1,".vs.",cond2,".jpg"), units = "cm", width = 25)
  ggsave(paste0(file_path, cond1,".vs.",cond2,".pdf"), units = "cm", width = 25)
  
  write.table(sig, file = paste0(file_path, "with_name_", cond1,".vs.",cond2,".txt"), sep = "\t", quote = F, row.names = F)
  write.table(sig %>% select(-c(protein_name)), 
              file = paste0(file_path,cond1,".vs.",cond2,".txt"), sep = "\t", quote = F, row.names = F)
  
  write.table(data_results, 
              file = paste0(file_path, "all_prots_with_name_", cond1,".vs.",cond2,".txt"), 
              sep = "\t", quote = F, row.names = F)
} #Cond1 = numerator, cond2 = denominator


data = combined_data
title = paste("Pathways in", conditions[1])
file_name = paste0("Pathways", conditions[1], ".png")
dir_path = "de_results_2024/proteomics/IPA_premod/"
data_dir_path = "de_results_2024/proteomics/IPA_premod/formatted_files/"
entries_of_interest_file_name = NA
ylab_substr = "Pathways"
max_count_filtered_pathways = 30

create_bar_plot <- function(data, title, file_name,
                            dir_path = "plots/IPA/",
                            data_dir_path = "IPA/files/",
                            ylab_substr = "Pathways",
                            entries_of_interest_file_name = NA,
                            max_count_filtered_pathways = 30){
  if(!is.na(entries_of_interest_file_name)){
    entries_of_interest <- read.table(entries_of_interest_file_name, header = FALSE, sep = "\t", col.names = "entries")
    data <- data %>%
      filter(pathways %in% entries_of_interest$entries)
  }
  if(!dir.exists(dir_path)){
    dir.create(dir_path, recursive = TRUE)
  }
  if(!dir.exists(data_dir_path)){
    dir.create(data_dir_path, recursive = TRUE)
  }
  
  data_to_plot <- data %>%
    mutate(zscore = as.numeric(zscore)) %>%
    arrange(desc(zscore)) 
  
  data_to_plot <- data_to_plot %>%
    mutate(pathways = factor(pathways, levels = rev(c(data_to_plot[, "pathways"]))))
  
  num_pathways <- length(data_to_plot$pathways)
  print(dim(data_to_plot))
  
  write.table(data_to_plot, paste0(data_dir_path, 
                                   sub(pattern = ".png", ".tsv", file_name)), 
              row.names = FALSE, sep = "\t")
  
  ggplot(data_to_plot, aes(x = adjpval, y = pathways, fill = zscore)) +
    geom_bar(stat = "identity") +
    scale_fill_viridis_c(option = "inferno") +
    xlab("-log10 adj.pval") +
    ylab(paste0(ylab_substr," (count = ", num_pathways, ")")) +
    geom_vline(xintercept = 1.301, linetype = 2) +
    labs(caption = "Vertical dashed line indicates -log10(0.05)") +
    ggtitle(title) +
    theme(axis.title.x = element_text(size=rel(1.2)),
          axis.title.y = element_text(size=rel(1.2)),
          plot.title  = element_text(size=rel(1.3)),
          plot.title.position = "plot")
  
  if(!dir.exists(dir_path)){
    dir.create(dir_path, recursive = TRUE)
  }
  ggsave(paste(dir_path, file_name), units = "cm", width = 25, height = 30)  
  
  if(nrow(data_to_plot) > max_count_filtered_pathways){

    count.activation <- ceiling((sum(data_to_plot$zscore > 0) / num_pathways) * max_count_filtered_pathways)
    count.inhibition <- max_count_filtered_pathways - count.activation
        
    temp <- data_to_plot %>%
      dplyr::filter(zscore > 0) %>%
      arrange(desc(adjpval))
    pathways.activated <- temp$pathways[c(1:count.activation)]
    
    temp <- data_to_plot %>%
      dplyr::filter(zscore < 0) %>%
      arrange(desc(adjpval))
    pathways.inhibited <- temp$pathways[c(1:count.inhibition)]
    
    data_to_plot.filt <- data_to_plot %>%
      dplyr::filter(pathways %in% c(pathways.activated, pathways.inhibited)) %>%
      arrange(desc(zscore)) 
    data_to_plot.filt <- data_to_plot.filt %>%
      mutate(pathways = factor(pathways, levels = rev(c(data_to_plot.filt[, "pathways"]))))
    print(dim(data_to_plot.filt))
    
    ggplot(data_to_plot.filt, aes(x = adjpval, y = pathways, fill = zscore)) +
      geom_bar(stat = "identity") +
      scale_fill_viridis_c(option = "inferno") +
      xlab("-log10 adj.pval") +
      ylab(paste0(ylab_substr," (count = ", max_count_filtered_pathways, ")")) +
      geom_vline(xintercept = 1.301, linetype = 2) +
      labs(caption = "Vertical dashed line indicates -log10(0.05)") +
      ggtitle(title) +
      theme(axis.title.x = element_text(size=rel(1.2)),
            axis.title.y = element_text(size=rel(1.2)),
            plot.title  = element_text(size=rel(1.3)),
            plot.title.position = "plot")
    
    fil_dir_path <- paste0(dir_path, "filtered/")
    
    if(!dir.exists(fil_dir_path)){
      dir.create(fil_dir_path, recursive = TRUE)
    }
    ggsave(paste(fil_dir_path, file_name), units = "cm", width = 25, height = 30)  
  }
}


# #create boxplots of DE prot/transcripts to show shift from CFRD to NGT
#de_file_path : file path of sig_<de file> with columns : Molecule, logFC, PVal, adjPVal
create_DE_boxplot <- function(data, phenotype, conditions_of_interest,
                              x_lab, output_dir_path,
                              de_file_path, 
                              k = 10){
  
  de_data <- read.table(de_file_path, sep = "\t", header = TRUE)  
  
  de.upreg <- de_data %>%
    arrange(desc(logFC))
  de.upreg <- de.upreg[c(1:k), "Molecule"]
  
  de.downreg <- de_data %>%
    arrange(logFC)
  de.downreg <- de.downreg[c(1:k), "Molecule"]
  
  phenotype.sub <- phenotype %>%
    dplyr::filter(modstatus_condition %in% conditions_of_interest)
  data.sub <- as.data.frame(t(data[, phenotype.sub$Sample])) %>%
    rownames_to_column("Sample") %>%
    inner_join(phenotype.sub %>% 
                 dplyr::select(c(Sample, modstatus_condition)) %>%
                 dplyr::rename(c("Condition" = "modstatus_condition"))) %>%
    dplyr::relocate(Condition, .after = Sample)
  
  if(!dir.exists(output_dir_path)){
    dir.create(path = output_dir_path, recursive = TRUE)
  }

  i <- 1
  for(de in list(de.upreg, de.downreg)){
    data_to_plot <- data.sub[, c("Sample", "Condition", de)]
    data_to_plot <- data_to_plot %>%
      pivot_longer(!c(Sample, Condition), names_to = "Molecule") %>%
      mutate(Molecule = factor(Molecule, levels = de)) %>%
      mutate(Condition = factor(Condition, levels = conditions_of_interest))
    
    if(i == 1){
      text <- "Upregulated"
    } else{
      text <- "Downregulated"
    }
    title <- paste("Top", k, text, x_lab)
    
    ggplot(data_to_plot, aes(x = Molecule, y = value)) +
      geom_boxplot(aes(fill = Condition)) +
      xlab(x_lab) +
      ylab("Expression") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = rel(1.2)),
            axis.title.x = element_text(size = rel(1.2)),
            axis.text.y = element_text(size = rel(1.2)),
            axis.title.y = element_text(size = rel(1.2)),
            legend.title = element_text(size = rel(1.2), face = "bold"),
            legend.text = element_text(size = rel(1.2)),
            plot.title = element_text(size = rel(1.4), face = "bold", hjust = 0.5)
            ) +
      ggtitle(title)
    
    output_file_path <- paste0(output_dir_path, title, ".png")
    ggsave(output_file_path)
    i <- i + 1
  }
}











