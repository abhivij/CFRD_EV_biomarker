library(tidyverse)
library(viridis)
library(umap)
library(ggrepel)
library(VennDiagram)

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

