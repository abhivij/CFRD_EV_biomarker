#helper functions related to DE analysis - to create volcano plots and write result files

library(EnhancedVolcano)
library(tidyverse)

# results = result
# k = 10
# plot_title = "PREOPE Vs MET"
# output_dir_path = "DE_results_2024/proteomics/1_only_condition/p/"
# plot_file_name = "PREOPEVsMET.png"
# fc_cutoff = 1.5
# pval_cutoff = 0.05
# use_adj_pval = FALSE
# x_lim = NA
# y_lim = NA
# molecule_names_file_path = "Data/Protein/formatted_data/all_protein_names.csv"
# molecule_names_file_columns = c(1, 3, 4)
# plot_width_cm = 25

#molecule_names_file_columns : c(<column_id of molecule_id>, <column_id of molecule_name>, 
#                                                            <optional : column_id of alternate molecule_id to be used>)
plot_volcano_and_save_DE <- function(
    results, plot_title, plot_file_name, output_dir_path,
    k = 10, fc_cutoff = 1.5, pval_cutoff = 0.05, use_adj_pval = FALSE,
    x_lim = NA, y_lim = NA,
    molecule_names_file_path = NA,
    molecule_names_file_columns = NA,
    plot_width_cm = 25
) {
  
  if(!dir.exists(output_dir_path)){
    dir.create(output_dir_path, recursive = TRUE)
  }
  
  lfc_cutoff <- log2(fc_cutoff)
  
  if(use_adj_pval){
    column_name <- "adjPVal"
  } else{
    column_name <- "PVal"
  }
  
  #arrange by significance needs to be done so as to show 
  # the labels in the order - downregulated, not significant, upregulated in the volcano plot 
  results <- results %>%
    dplyr::rename(c("pvalcolumn" = column_name)) %>%
    mutate(significance = case_when(logFC <= -lfc_cutoff & pvalcolumn <= pval_cutoff ~ 'Downregulated',
                                    logFC >= lfc_cutoff & pvalcolumn <= pval_cutoff ~ 'Upregulated',
                                    TRUE ~ 'Not significant')) %>% 
    mutate(colour = case_when(significance == 'Downregulated' ~ '#38ACE2',
                              significance == 'Upregulated' ~ '#E8495C',
                              TRUE ~ 'grey')) %>%
    arrange(significance)
  
  if(use_adj_pval){
    results <- results %>%
      dplyr::rename(c("adjPVal" = "pvalcolumn")) 
    y_lab <- bquote(~-Log[10] ~ italic(adjP))
    p_val_caption <- 'adj p value cutoff :'
  } else{
    results <- results %>%
      dplyr::rename(c("PVal" = "pvalcolumn")) 
    y_lab <- bquote(~-Log[10] ~ italic(P))
    p_val_caption <- 'p value cutoff :'
  }
  
  #include protein names
  if(!is.na(molecule_names_file_path)){
    molecule_names <- read.csv(molecule_names_file_path) %>%
      dplyr::select(all_of(molecule_names_file_columns))
    
    if(length(molecule_names_file_columns) == 3){
      colnames(molecule_names) <- c("Molecule", "MoleculeName", "name")  
    } else if(length(molecule_names_file_columns) == 2){
      colnames(molecule_names) <- c("Molecule", "MoleculeName")      
    } else{
      print("invalid molecule_names_file_columns")
      return
    }
    
    results <- results %>%
      left_join(molecule_names, by = "Molecule")
    if("name" %in% colnames(results)){
      results <- results %>%
        mutate(Molecule = case_when(!is.na(name) ~ name,
                                    TRUE ~ Molecule)) %>%
        dplyr::select(-c(name))
    }
  }
  
  upreg <- results %>%
    filter(significance == 'Upregulated') %>%
    arrange(desc(logFC))
  downreg <- results %>%
    filter(significance == "Downregulated") %>%
    arrange(logFC)
  
  if(nrow(upreg) >= k){
    top_proteins <- upreg[1:k, "Molecule"]
  } else{
    top_proteins <- upreg[1:nrow(upreg), "Molecule"]
  }
  if(nrow(downreg) >= k){
    top_proteins <- c(top_proteins, downreg[1:k, "Molecule"])
  } else{
    top_proteins <- c(top_proteins, downreg[1:nrow(downreg), "Molecule"])
  }
  
  results <- results %>%
    mutate(label = case_when(Molecule %in% top_proteins ~ Molecule,
                             TRUE ~ NA_character_))
  
  sig <- rbind(upreg, downreg) %>%
    select(-c(significance, colour)) %>%
    arrange(desc(logFC))
  # print(sig)
  # print(paste("Upregulated:", nrow(sig[sig$logFC > 0,]), sep=" "))
  # print(paste("Downregulated:", nrow(sig[sig$logFC < 0,]), sep=" "))
  
  keyvals <- results$colour
  names(keyvals) <- results$significance
  
  upreg_count <- sum(names(keyvals) == 'Upregulated')
  downreg_count <- sum(names(keyvals) == 'Downregulated')
  ns_count <- sum(names(keyvals) == 'Not significant')
  
  names(keyvals) <- gsub('Upregulated', paste0('Upregulated (', upreg_count, ')'), names(keyvals))
  names(keyvals) <- gsub('Downregulated', paste0('Downregulated (', downreg_count, ')'), names(keyvals))
  names(keyvals) <- gsub('Not significant', paste0('Not significant (', ns_count, ')'), names(keyvals))
  
  if(is.na(x_lim)){
    x_lim <- c(min(results$logFC)-0.5,
               max(results$logFC)+0.5)
  }
  if(is.na(y_lim)){
    y_lim <- c(0, max(-log10(results$PVal))+0.5)
  }
  
  volcanoplot <- EnhancedVolcano(results,
                                 lab = results$label,
                                 x = 'logFC',
                                 y = column_name,
                                 xlim = x_lim,
                                 ylim = y_lim,
                                 ylab = y_lab,
                                 pCutoff = pval_cutoff,
                                 FCcutoff = lfc_cutoff,
                                 colCustom = keyvals,
                                 colAlpha = 1,
                                 title = plot_title,
                                 legendPosition = 'bottom',
                                 subtitle = '',
                                 caption = paste(p_val_caption, pval_cutoff,
                                                 '   ',
                                                 'fold change cutoff :', fc_cutoff),
                                 labSize = 3,
                                 drawConnectors = T,
                                 max.overlaps = Inf,
                                 arrowheads = F)
  
  plot(volcanoplot)
  ggsave(paste0(output_dir_path, plot_file_name), units = "cm", width = plot_width_cm)
  # ggsave(paste0(file_path, cond1,".vs.",cond2,".pdf"), units = "cm", width = 25)
  
  de_results_file_name <- sub(".png", ".csv", plot_file_name)
  de_results_file_name <- sub(".jpg", ".csv", de_results_file_name)
  de_results_file_name <- sub(".pdf", ".csv", de_results_file_name)
  
  write.table(sig, file = paste0(output_dir_path, 
                                 "sig_", de_results_file_name), sep = "\t", quote = F, row.names = F)
  if("MoleculeName" %in% colnames(sig)){
    write.table(sig %>% select(-c(MoleculeName)), 
                file = paste0(output_dir_path, "sig_no_name_", de_results_file_name), 
                sep = "\t", quote = F, row.names = F)    
  }
  write.table(results, file = paste0(output_dir_path, 
                                     "all_", de_results_file_name), sep = "\t", quote = F, row.names = F)
}


