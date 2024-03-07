library(GOplot)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rrvgo)

go_BP_plot <- function(go.bp, file_path){
  
  # go.bp <- de.down.ego.bp.result
  # file_path <- paste0(plot_dir_path, "down_", file_name, ".pdf")
  
  simMatrix <- calculateSimMatrix(go.bp$ID,
                                  orgdb="org.Hs.eg.db",
                                  ont="BP",
                                  method="Rel")
  scores <- setNames(-log10(go.bp$qvalue), go.bp$ID)
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.7,
                                  orgdb="org.Hs.eg.db")
  length(factor(reducedTerms$parentTerm))
  
  pdf(file = file_path, height = 10, width = 10)
  scatterPlot(simMatrix, reducedTerms)
  dev.off()
  
  pdf(file = sub("goBPscatter", "goBPtree", file_path), height = 10, width = 10)
  treemapPlot(reducedTerms)
  dev.off()
}


dir_path <- "de_results_2024/proteomics/premod/"

for(file_name in c("PreModulator_CFRDVsPreModulator_IGT", "PreModulator_CFRDVsPreModulator_NGT",
                     "PreModulator_IGTVsPreModulator_NGT")){
  print(file_name)
  # file_name <- "PreModulator_CFRDVsPreModulator_IGT"
  full_file_name <- paste0("sig_no_name_", file_name, ".csv")
  
  de <- read.table(paste0(dir_path, 
                          "p/", 
                          full_file_name), 
                   sep = "\t", header = TRUE)
  de.up <- de %>%
    filter(logFC > 0)
  de.down <- de %>%
    filter(logFC < 0)
  
  de.up.geneids <- mapIds(org.Hs.eg.db, de.up$Molecule, 'ENTREZID', 'SYMBOL')
  de.down.geneids <- mapIds(org.Hs.eg.db, de.down$Molecule, 'ENTREZID', 'SYMBOL')
  
  de.up.ego.bp <- enrichGO(
    gene          = de.up.geneids,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pvalueCutoff  = 0.05
  )
  de.down.ego.bp <- enrichGO(
    gene          = de.down.geneids,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pvalueCutoff  = 0.05
  )
  
  de.up.ego.bp.result <- de.up.ego.bp@result %>%
    filter(qvalue < 0.05)
  rownames(de.up.ego.bp.result) <- c()
  
  de.down.ego.bp.result <- de.down.ego.bp@result %>%
    filter(qvalue < 0.05)
  rownames(de.down.ego.bp.result) <- c()
  
  plot_dir_path <- paste0(dir_path, "goBPtree/")
  if(!dir.exists(plot_dir_path)){
    dir.create(plot_dir_path, recursive = TRUE)
  }
  
  plot_dir_path <- paste0(dir_path, "goBPscatter/")
  if(!dir.exists(plot_dir_path)){
    dir.create(plot_dir_path, recursive = TRUE)
  }
  
  print(dim(de.up.ego.bp.result))
  if(nrow(de.up.ego.bp.result) > 2){
    go.bp <- de.up.ego.bp.result
    simMatrix <- calculateSimMatrix(go.bp$ID,
                                    orgdb="org.Hs.eg.db",
                                    ont="BP",
                                    method="Rel")
    scores <- setNames(-log10(go.bp$qvalue), go.bp$ID)
    reducedTerms <- reduceSimMatrix(simMatrix,
                                    scores,
                                    threshold=0.7,
                                    orgdb="org.Hs.eg.db")
    length(factor(reducedTerms$parentTerm))
    
    file_path <- paste0(plot_dir_path, "up_", file_name, ".pdf")
    
    pdf(file = file_path, height = 10, width = 10)
    scatterPlot(simMatrix, reducedTerms)
    dev.off()
    
    pdf(file = sub("goBPscatter", "goBPtree", file_path), height = 10, width = 10)
    treemapPlot(reducedTerms)
    dev.off()
  }
  
  print(dim(de.down.ego.bp.result))
  if(nrow(de.down.ego.bp.result) > 2){
    go.bp <- de.down.ego.bp.result
    simMatrix <- calculateSimMatrix(go.bp$ID,
                                    orgdb="org.Hs.eg.db",
                                    ont="BP",
                                    method="Rel")
    scores <- setNames(-log10(go.bp$qvalue), go.bp$ID)
    reducedTerms <- reduceSimMatrix(simMatrix,
                                    scores,
                                    threshold=0.7,
                                    orgdb="org.Hs.eg.db")
    length(factor(reducedTerms$parentTerm))
    
    file_path <- paste0(plot_dir_path, "down_", file_name, ".pdf")
    
    pdf(file = file_path, height = 10, width = 10)
    scatterPlot(simMatrix, reducedTerms)
    dev.off()
    
    pdf(file = sub("goBPscatter", "goBPtree", file_path), height = 10, width = 10)
    treemapPlot(reducedTerms)
    dev.off()
  }
  
}



dir_path <- "de_results_2024/proteomics/postmod_premod/"

for(file_name in c("PostModulator_CFRDVsPreModulator_CFRD", "PostModulator_IGTVsPreModulator_IGT",
                   "PostModulator_NGTVsPreModulator_NGT", "PostModulatorVsPreModulator")){
  print(file_name)
  # file_name <- "PreModulator_CFRDVsPreModulator_IGT"
  full_file_name <- paste0("sig_no_name_", file_name, ".csv")
  
  de <- read.table(paste0(dir_path, 
                          "p/", 
                          full_file_name), 
                   sep = "\t", header = TRUE)
  de.up <- de %>%
    filter(logFC > 0)
  de.down <- de %>%
    filter(logFC < 0)
  
  de.up.geneids <- mapIds(org.Hs.eg.db, de.up$Molecule, 'ENTREZID', 'SYMBOL')
  de.down.geneids <- mapIds(org.Hs.eg.db, de.down$Molecule, 'ENTREZID', 'SYMBOL')
  
  de.up.ego.bp <- enrichGO(
    gene          = de.up.geneids,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pvalueCutoff  = 0.05
  )
  de.down.ego.bp <- enrichGO(
    gene          = de.down.geneids,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pvalueCutoff  = 0.05
  )
  
  de.up.ego.bp.result <- de.up.ego.bp@result %>%
    filter(qvalue < 0.05)
  rownames(de.up.ego.bp.result) <- c()
  
  de.down.ego.bp.result <- de.down.ego.bp@result %>%
    filter(qvalue < 0.05)
  rownames(de.down.ego.bp.result) <- c()
  
  plot_dir_path <- paste0(dir_path, "goBPtree/")
  if(!dir.exists(plot_dir_path)){
    dir.create(plot_dir_path, recursive = TRUE)
  }
  
  plot_dir_path <- paste0(dir_path, "goBPscatter/")
  if(!dir.exists(plot_dir_path)){
    dir.create(plot_dir_path, recursive = TRUE)
  }
  
  print(dim(de.up.ego.bp.result))
  if(nrow(de.up.ego.bp.result) > 2){
    go.bp <- de.up.ego.bp.result
    simMatrix <- calculateSimMatrix(go.bp$ID,
                                    orgdb="org.Hs.eg.db",
                                    ont="BP",
                                    method="Rel")
    scores <- setNames(-log10(go.bp$qvalue), go.bp$ID)
    reducedTerms <- reduceSimMatrix(simMatrix,
                                    scores,
                                    threshold=0.7,
                                    orgdb="org.Hs.eg.db")
    length(factor(reducedTerms$parentTerm))
    
    file_path <- paste0(plot_dir_path, "up_", file_name, ".pdf")
    
    pdf(file = file_path, height = 10, width = 10)
    scatterPlot(simMatrix, reducedTerms)
    dev.off()
    
    pdf(file = sub("goBPscatter", "goBPtree", file_path), height = 10, width = 10)
    treemapPlot(reducedTerms)
    dev.off()
  }
  
  print(dim(de.down.ego.bp.result))
  if(nrow(de.down.ego.bp.result) > 2){
    go.bp <- de.down.ego.bp.result
    simMatrix <- calculateSimMatrix(go.bp$ID,
                                    orgdb="org.Hs.eg.db",
                                    ont="BP",
                                    method="Rel")
    scores <- setNames(-log10(go.bp$qvalue), go.bp$ID)
    reducedTerms <- reduceSimMatrix(simMatrix,
                                    scores,
                                    threshold=0.7,
                                    orgdb="org.Hs.eg.db")
    length(factor(reducedTerms$parentTerm))
    
    file_path <- paste0(plot_dir_path, "down_", file_name, ".pdf")
    
    pdf(file = file_path, height = 10, width = 10)
    scatterPlot(simMatrix, reducedTerms)
    dev.off()
    
    pdf(file = sub("goBPscatter", "goBPtree", file_path), height = 10, width = 10)
    treemapPlot(reducedTerms)
    dev.off()
  }
  
}



#############################
#IPA

#canonical pathways
comparison_adjpval <- read.table(file = "de_results_2024/proteomics/IPA_premod/comparison_analysis/can_path_bhpval.txt",
                                 sep="\t", header=T, row.names=1, skip = 2, quote = "")
comparison_zscore <- read.table(file = "de_results_2024/proteomics/IPA_premod/comparison_analysis/can_path_zscore.txt",
                                sep="\t", header=T, row.names=1, skip = 2, quote = "",
                                colClasses = c("character", "numeric", "numeric", "numeric"),
                                na.strings = c("N/A"))
conditions <- c("Premodulator CFRD Vs Premodulator IGT", 
                "Premodulator CFRD Vs Premodulator NGT", 
                "Premodulator IGT Vs Premodulator NGT")
colnames(comparison_adjpval) <- conditions
colnames(comparison_zscore) <- conditions

#Premodulator CFRD Vs Premodulator IGT
combined_data <- comparison_zscore %>%
  dplyr::rename(zscore = conditions[1]) %>%
  select(zscore) %>%
  rownames_to_column("pathways") %>%
  filter(!(is.na(zscore))) %>%
  inner_join(
    comparison_adjpval %>%
      dplyr::rename(adjpval = conditions[1]) %>%
      select(adjpval) %>%
      rownames_to_column("pathways") %>%
      filter(adjpval > -log10(0.05))
  )
create_bar_plot(data = combined_data,
                title = paste("Pathways in", conditions[1]),
                file_name = paste0("Pathways", conditions[1], ".png"),
                dir_path = "de_results_2024/proteomics/IPA_premod/", 
                data_dir_path = "de_results_2024/proteomics/IPA_premod/formatted_files/")

#Premodulator CFRD Vs Premodulator NGT
combined_data <- comparison_zscore %>%
  dplyr::rename(zscore = conditions[2]) %>%
  select(zscore) %>%
  rownames_to_column("pathways") %>%
  filter(!(is.na(zscore))) %>%
  inner_join(
    comparison_adjpval %>%
      dplyr::rename(adjpval = conditions[2]) %>%
      select(adjpval) %>%
      rownames_to_column("pathways") %>%
      filter(adjpval > -log10(0.05))
  )
create_bar_plot(data = combined_data,
                title = paste("Pathways in", conditions[2]),
                file_name = paste0("Pathways", conditions[2], ".png"),
                dir_path = "de_results_2024/proteomics/IPA_premod/", 
                data_dir_path = "de_results_2024/proteomics/IPA_premod/formatted_files/")

#Premodulator IGT Vs Premodulator NGT
combined_data <- comparison_zscore %>%
  dplyr::rename(zscore = conditions[3]) %>%
  select(zscore) %>%
  rownames_to_column("pathways") %>%
  filter(!(is.na(zscore))) %>%
  inner_join(
    comparison_adjpval %>%
      dplyr::rename(adjpval = conditions[3]) %>%
      select(adjpval) %>%
      rownames_to_column("pathways") %>%
      filter(adjpval > -log10(1))
  )
create_bar_plot(data = combined_data,
                title = paste("Pathways in", conditions[3]),
                file_name = paste0("Pathways", conditions[3], ".png"),
                dir_path = "de_results_2024/proteomics/IPA_premod/", 
                data_dir_path = "de_results_2024/proteomics/IPA_premod/formatted_files/")



comparison_adjpval <- read.table(file = "de_results_2024/proteomics/IPA_postmod_premod/comparison_analysis/can_path_bhpval.txt",
                                 sep="\t", header=T, row.names=1, skip = 2, quote = "")
comparison_zscore <- read.table(file = "de_results_2024/proteomics/IPA_postmod_premod/comparison_analysis/can_path_zscore.txt",
                                sep="\t", header=T, row.names=1, skip = 2, quote = "",
                                colClasses = c("character", "numeric", "numeric", "numeric"),
                                na.strings = c("N/A"))
conditions <- c("Postmodulator Vs Premodulator", 
                "Postmodulator CFRD Vs Premodulator CFRD", 
                "Postmodulator IGT Vs Premodulator IGT",
                "Postmodulator NGT Vs Premodulator NGT")
colnames(comparison_adjpval) <- conditions
colnames(comparison_zscore) <- conditions


combined_data <- comparison_zscore %>%
  dplyr::rename(zscore = conditions[1]) %>%
  select(zscore) %>%
  rownames_to_column("pathways") %>%
  filter(!(is.na(zscore))) %>%
  inner_join(
    comparison_adjpval %>%
      dplyr::rename(adjpval = conditions[1]) %>%
      select(adjpval) %>%
      rownames_to_column("pathways") %>%
      filter(adjpval > -log10(0.05))
  )
create_bar_plot(data = combined_data,
                title = paste("Pathways in", conditions[1]),
                file_name = paste0("Pathways", conditions[1], ".png"),
                dir_path = "de_results_2024/proteomics/IPA_postmod_premod/", 
                data_dir_path = "de_results_2024/proteomics/IPA_postmod_premod/formatted_files/")


combined_data <- comparison_zscore %>%
  dplyr::rename(zscore = conditions[2]) %>%
  select(zscore) %>%
  rownames_to_column("pathways") %>%
  filter(!(is.na(zscore))) %>%
  inner_join(
    comparison_adjpval %>%
      dplyr::rename(adjpval = conditions[2]) %>%
      select(adjpval) %>%
      rownames_to_column("pathways") %>%
      filter(adjpval > -log10(0.05))
  )
create_bar_plot(data = combined_data,
                title = paste("Pathways in", conditions[2]),
                file_name = paste0("Pathways", conditions[2], ".png"),
                dir_path = "de_results_2024/proteomics/IPA_postmod_premod/", 
                data_dir_path = "de_results_2024/proteomics/IPA_postmod_premod/formatted_files/")


combined_data <- comparison_zscore %>%
  dplyr::rename(zscore = conditions[3]) %>%
  select(zscore) %>%
  rownames_to_column("pathways") %>%
  filter(!(is.na(zscore))) %>%
  inner_join(
    comparison_adjpval %>%
      dplyr::rename(adjpval = conditions[3]) %>%
      select(adjpval) %>%
      rownames_to_column("pathways") %>%
      filter(adjpval > -log10(0.05))
  )
create_bar_plot(data = combined_data,
                title = paste("Pathways in", conditions[3]),
                file_name = paste0("Pathways", conditions[3], ".png"),
                dir_path = "de_results_2024/proteomics/IPA_postmod_premod/", 
                data_dir_path = "de_results_2024/proteomics/IPA_postmod_premod/formatted_files/")


combined_data <- comparison_zscore %>%
  dplyr::rename(zscore = conditions[4]) %>%
  select(zscore) %>%
  rownames_to_column("pathways") %>%
  filter(!(is.na(zscore))) %>%
  inner_join(
    comparison_adjpval %>%
      dplyr::rename(adjpval = conditions[4]) %>%
      select(adjpval) %>%
      rownames_to_column("pathways") %>%
      filter(adjpval > -log10(0.05))
  )
create_bar_plot(data = combined_data,
                title = paste("Pathways in", conditions[4]),
                file_name = paste0("Pathways", conditions[4], ".png"),
                dir_path = "de_results_2024/proteomics/IPA_postmod_premod/", 
                data_dir_path = "de_results_2024/proteomics/IPA_postmod_premod/formatted_files/")


#diseases and biofunctions
comparison_adjpval <- read.table(file = "de_results_2024/proteomics/IPA_premod/comparison_analysis/diseases_bhpval.txt",
                                 sep="\t", header=T, row.names=1, skip = 2, quote = "")
comparison_zscore <- read.table(file = "de_results_2024/proteomics/IPA_premod/comparison_analysis/diseases_zscore.txt",
                                sep="\t", header=T, row.names=1, skip = 2, quote = "",
                                colClasses = c("character", "numeric", "numeric", "numeric"),
                                na.strings = c("N/A"))
conditions <- c("Premodulator CFRD Vs Premodulator IGT", 
                "Premodulator CFRD Vs Premodulator NGT", 
                "Premodulator IGT Vs Premodulator NGT")
colnames(comparison_adjpval) <- conditions
colnames(comparison_zscore) <- conditions

#Premodulator CFRD Vs Premodulator IGT
combined_data <- comparison_zscore %>%
  dplyr::rename(zscore = conditions[1]) %>%
  select(zscore) %>%
  rownames_to_column("pathways") %>%
  filter(!(is.na(zscore))) %>%
  inner_join(
    comparison_adjpval %>%
      dplyr::rename(adjpval = conditions[1]) %>%
      select(adjpval) %>%
      rownames_to_column("pathways") %>%
      filter(adjpval > -log10(0.05))
  )
create_bar_plot(data = combined_data, ylab_substr = "Diseases & Biofunctions",
                title = paste("Diseases & Biofunctions in", conditions[1]),
                file_name = paste0("Diseases & Biofunctions", conditions[1], ".png"),
                dir_path = "de_results_2024/proteomics/IPA_premod/", 
                data_dir_path = "de_results_2024/proteomics/IPA_premod/formatted_files/")

#Premodulator CFRD Vs Premodulator NGT
combined_data <- comparison_zscore %>%
  dplyr::rename(zscore = conditions[2]) %>%
  select(zscore) %>%
  rownames_to_column("pathways") %>%
  filter(!(is.na(zscore))) %>%
  inner_join(
    comparison_adjpval %>%
      dplyr::rename(adjpval = conditions[2]) %>%
      select(adjpval) %>%
      rownames_to_column("pathways") %>%
      filter(adjpval > -log10(0.05))
  )
create_bar_plot(data = combined_data, ylab_substr = "Diseases & Biofunctions",
                title = paste("Diseases & Biofunctions in", conditions[2]),
                file_name = paste0("Diseases & Biofunctions", conditions[2], ".png"),
                dir_path = "de_results_2024/proteomics/IPA_premod/", 
                data_dir_path = "de_results_2024/proteomics/IPA_premod/formatted_files/")

#Premodulator IGT Vs Premodulator NGT
combined_data <- comparison_zscore %>%
  dplyr::rename(zscore = conditions[3]) %>%
  select(zscore) %>%
  rownames_to_column("pathways") %>%
  filter(!(is.na(zscore))) %>%
  inner_join(
    comparison_adjpval %>%
      dplyr::rename(adjpval = conditions[3]) %>%
      select(adjpval) %>%
      rownames_to_column("pathways") %>%
      filter(adjpval > -log10(1))
  )
create_bar_plot(data = combined_data, ylab_substr = "Diseases & Biofunctions",
                title = paste("Diseases & Biofunctions in", conditions[3]),
                file_name = paste0("Diseases & Biofunctions", conditions[3], ".png"),
                dir_path = "de_results_2024/proteomics/IPA_premod/", 
                data_dir_path = "de_results_2024/proteomics/IPA_premod/formatted_files/")



comparison_adjpval <- read.table(file = "de_results_2024/proteomics/IPA_postmod_premod/comparison_analysis/diseases_bhpval.txt",
                                 sep="\t", header=T, row.names=1, skip = 2, quote = "")
comparison_zscore <- read.table(file = "de_results_2024/proteomics/IPA_postmod_premod/comparison_analysis/diseases_zscore.txt",
                                sep="\t", header=T, row.names=1, skip = 2, quote = "",
                                colClasses = c("character", "numeric", "numeric", "numeric"),
                                na.strings = c("N/A"))
conditions <- c("Postmodulator Vs Premodulator", 
                "Postmodulator CFRD Vs Premodulator CFRD", 
                "Postmodulator IGT Vs Premodulator IGT",
                "Postmodulator NGT Vs Premodulator NGT")
colnames(comparison_adjpval) <- conditions
colnames(comparison_zscore) <- conditions


combined_data <- comparison_zscore %>%
  dplyr::rename(zscore = conditions[1]) %>%
  select(zscore) %>%
  rownames_to_column("pathways") %>%
  filter(!(is.na(zscore))) %>%
  inner_join(
    comparison_adjpval %>%
      dplyr::rename(adjpval = conditions[1]) %>%
      select(adjpval) %>%
      rownames_to_column("pathways") %>%
      filter(adjpval > -log10(0.05))
  )
create_bar_plot(data = combined_data, ylab_substr = "Diseases & Biofunctions",
                title = paste("Diseases & Biofunctions in", conditions[1]),
                file_name = paste0("Diseases & Biofunctions", conditions[1], ".png"),
                dir_path = "de_results_2024/proteomics/IPA_postmod_premod/", 
                data_dir_path = "de_results_2024/proteomics/IPA_postmod_premod/formatted_files/")


combined_data <- comparison_zscore %>%
  dplyr::rename(zscore = conditions[2]) %>%
  select(zscore) %>%
  rownames_to_column("pathways") %>%
  filter(!(is.na(zscore))) %>%
  inner_join(
    comparison_adjpval %>%
      dplyr::rename(adjpval = conditions[2]) %>%
      select(adjpval) %>%
      rownames_to_column("pathways") %>%
      filter(adjpval > -log10(0.05))
  )
create_bar_plot(data = combined_data, ylab_substr = "Diseases & Biofunctions",
                title = paste("Diseases & Biofunctions in", conditions[2]),
                file_name = paste0("Diseases & Biofunctions", conditions[2], ".png"),
                dir_path = "de_results_2024/proteomics/IPA_postmod_premod/", 
                data_dir_path = "de_results_2024/proteomics/IPA_postmod_premod/formatted_files/")


combined_data <- comparison_zscore %>%
  dplyr::rename(zscore = conditions[3]) %>%
  select(zscore) %>%
  rownames_to_column("pathways") %>%
  filter(!(is.na(zscore))) %>%
  inner_join(
    comparison_adjpval %>%
      dplyr::rename(adjpval = conditions[3]) %>%
      select(adjpval) %>%
      rownames_to_column("pathways") %>%
      filter(adjpval > -log10(0.05))
  )
create_bar_plot(data = combined_data, ylab_substr = "Diseases & Biofunctions",
                title = paste("Diseases & Biofunctions in", conditions[3]),
                file_name = paste0("Diseases & Biofunctions", conditions[3], ".png"),
                dir_path = "de_results_2024/proteomics/IPA_postmod_premod/", 
                data_dir_path = "de_results_2024/proteomics/IPA_postmod_premod/formatted_files/")


combined_data <- comparison_zscore %>%
  dplyr::rename(zscore = conditions[4]) %>%
  select(zscore) %>%
  rownames_to_column("pathways") %>%
  filter(!(is.na(zscore))) %>%
  inner_join(
    comparison_adjpval %>%
      dplyr::rename(adjpval = conditions[4]) %>%
      select(adjpval) %>%
      rownames_to_column("pathways") %>%
      filter(adjpval > -log10(0.05))
  )
create_bar_plot(data = combined_data, ylab_substr = "Diseases & Biofunctions",
                title = paste("Diseases & Biofunctions in", conditions[4]),
                file_name = paste0("Diseases & Biofunctions", conditions[4], ".png"),
                dir_path = "de_results_2024/proteomics/IPA_postmod_premod/", 
                data_dir_path = "de_results_2024/proteomics/IPA_postmod_premod/formatted_files/")
