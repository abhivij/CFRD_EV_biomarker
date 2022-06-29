library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(multiMiR)

# geneList.uniprot <- unique(gsub("-[0-9]+", "", names(geneList_)))
# map <- select(org.Hs.eg.db, geneList.uniprot, "ENTREZID", "UNIPROT")
# geneList <- na.omit(map$ENTREZID)
# 
# ego <- enrichGO(gene           = geneList,
#                 OrgDb         = org.Hs.eg.db,
#                 #keyType       = 'UNIPROT',
#                 ont           = "BP",
#                 pAdjustMethod = "fdr",
#                 pvalueCutoff  = 0.05)


####################################CF Vs HC#################################### 

input_file_path = "de_results/DE_CFVsHC.csv"

go_BP(input_file_path = "de_results/DE_CFVsHC.csv", 
      output_dir = "functional_analysis", 
      output_file_name = "cfvshc_go_bp", plot_title = "CF Vs HC")

go_BP(input_file_path = "de_results/DE_CFRDVsIGT.csv", 
      output_dir = "functional_analysis", 
      output_file_name = "cfrdvsigt_go_bp", plot_title = "CFRD Vs IGT")

go_BP(input_file_path = "de_results/DE_CFRDVsNGT.csv", 
      output_dir = "functional_analysis", 
      output_file_name = "cfrdvsngt_go_bp", plot_title = "CFRD Vs NGT")

go_BP(input_file_path = "de_results/DE_IGTVsNGT.csv", 
      output_dir = "functional_analysis", 
      output_file_name = "igtvsngt_go_bp", plot_title = "IGT Vs NGT")



obtain_mirna_go_bp <- function(mirna_vec,
                               plot_title,
                               output_dir,
                               output_plot_name,
                               output_result_file_name){
  
  mir.targets <- get_multimir(
    org     = "hsa",
    mirna  = mirna_vec,
    table   = "validated",
    summary = TRUE
  )  
  
  mir.targets.data <- mir.targets@data[, c(2, 3, 5)] %>%
    unique() %>%
    filter(target_entrez != "")
  
  mir.targets.data.grouped <- mir.targets.data %>%
    group_by(target_entrez) %>%
    summarise(n = n()) %>%
    arrange(desc(n))
  
  gene_list <- unique(mir.targets.data$target_entrez)
  
  ego.bp <- enrichGO(
    gene          = gene_list,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pvalueCutoff  = 0.05
  )
  ego.bp.result <- ego.bp@result
  
  options(scipen = 2, digits = 3)
  
  barplot(ego.bp, x = "Count", showCategory = 15,
          title = plot_title)
  ggsave(paste0(output_dir, "/", output_plot_name), width = 20)
  # heatplot(ego.bp, showCategory = 10)
  
  write.csv(ego.bp.result, paste0(output_dir, "/", output_result_file_name), row.names = FALSE)  
}


go_BP <- function(input_file_path, output_dir, output_file_name, plot_title){
  de <- read.csv(input_file_path) %>%
    filter(FDR < 0.05)
  
  de.mir <- de %>%
    filter(grepl("miR", Molecule) | grepl("hsa-let", Molecule))
  de.pir <- de %>%
    filter(grepl("piR", Molecule))
  
  assertthat::are_equal(
    dim(
      de %>%
        filter(!Molecule %in% de.mir$Molecule) %>%
        filter(!Molecule %in% de.pir$Molecule)
    )[1], 
    0)
  
  de.mir.up <- de.mir %>%
    filter(logFC > 0)
  print("up reg FDR < 0.05")
  print(dim(de.mir.up))
  de.mir.down <- de.mir %>%
    filter(logFC < 0)
  print("down reg FDR < 0.05")
  print(dim(de.mir.down))

  obtain_mirna_go_bp(mirna_vec = de.mir.up$Molecule,
                     plot_title = paste("Biological Process : upregulated miRNAs", plot_title),
                     output_dir = output_dir,
                     output_plot_name = paste0(output_file_name, "_up.png"),
                     output_result_file_name = paste0(output_file_name, "_up.csv"))  
  
  obtain_mirna_go_bp(mirna_vec = de.mir.down$Molecule,
                     plot_title = paste("Biological Process : downregulated miRNAs", plot_title),
                     output_dir = output_dir,
                     output_plot_name = paste0(output_file_name, "_down.png"),
                     output_result_file_name = paste0(output_file_name, "_down.csv"))
  
}





de <- read.csv("de_results/DE_CFVsHC.csv") %>%
  filter(FDR < 0.05)
de.mir <- de %>%
  filter(grepl("miR", Molecule) | grepl("hsa-let", Molecule))
de.pir <- de %>%
  filter(grepl("piR", Molecule))

assertthat::are_equal(
  dim(
    de %>%
      filter(!Molecule %in% de.mir$Molecule) %>%
      filter(!Molecule %in% de.pir$Molecule)
    )[1], 
  0)

de.mir.up <- de.mir %>%
  filter(logFC > 0) %>%
  arrange(desc(logFC))
de.mir.down <- de.mir %>%
  filter(logFC < 0) %>%
  arrange(desc(-logFC))


for(mir in de.mir.up[c(4:5), "Molecule"]){
  print(mir)
  obtain_mirna_go_bp(mirna_vec = mir,
                     plot_title = paste("Biological Process : upregulated :", mir, "CF Vs HC"),
                     output_dir = "functional_analysis",
                     output_plot_name = paste0("CFVsHC_", mir, "_up.png"),
                     output_result_file_name = paste0("CFVsHC_", mir, "_up.csv"))  
}

for(mir in de.mir.down[c(1:5), "Molecule"]){
  print(mir)
  obtain_mirna_go_bp(mirna_vec = mir,
                     plot_title = paste("Biological Process : downregulated :", mir, "CF Vs HC"),
                     output_dir = "functional_analysis",
                     output_plot_name = paste0("CFVsHC_", mir, "_down.png"),
                     output_result_file_name = paste0("CFVsHC_", mir, "_down.csv"))  
}


###########################
de.mir.down.targets <- get_multimir(
  org     = "hsa",
  mirna  = de.mir.down$Molecule,
  table   = "validated",
  summary = TRUE
)
de.mir.down.targets@summary


#### not removing the below commented code to show the reason why only single identifier i.e. entrez gene id was used

# de.mir.down.targets.data <- de.mir.down.targets@data[, c(2:5)] %>%
#   unique() %>%
#   filter(target_entrez != "")

# grouped.de.mir.down.targets.data <- de.mir.down.targets.data %>%
#   group_by(target_symbol, target_entrez) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   arrange(desc(n))
# 
# grouped2.de.mir.down.targets.data <- grouped.de.mir.down.targets.data %>%
#   group_by(target_entrez) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   arrange(desc(n))
# 
# target_entrez_multiple_entries <- grouped2.de.mir.down.targets.data %>%
#   filter(n > 1) 
# 
# grouped.de.mir.down.targets.data %>%
#   filter(target_entrez %in% target_entrez_multiple_entries$target_entrez) %>%
#   arrange(target_entrez)

# gene_list <- unique(de.mir.down.targets.data$target_entrez)

#entrez gene id and gene symbol are not bijective
#using entrez gene id only


de.mir.down.targets.data <- de.mir.down.targets@data[, c(2, 3, 5)] %>%
  unique() %>%
  filter(target_entrez != "")

de.mir.down.targets.data.grouped <- de.mir.down.targets.data %>%
  group_by(target_entrez) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

gene_list <- unique(de.mir.down.targets.data$target_entrez)

ego.bp <- enrichGO(
  gene          = gene_list,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pvalueCutoff  = 0.05
)
ego.bp.result <- ego.bp@result

options(scipen = 2, digits = 3)

barplot(ego.bp, x = "Count", showCategory = 15,
        title = "Biological Process : CF Vs HC downregulated miRNAs")
ggsave("functional_analysis/cfvshc_down_go_bp.png")
# heatplot(ego.bp, showCategory = 10)

write.csv(ego.bp.result, "functional_analysis/cfvshc_down_go_bp.csv", row.names = FALSE)

ek <- enrichKEGG(gene          = gene_list,
                 organism = "hsa",
                 pvalueCutoff  = 0.05)
ek.result <- ek@result

barplot(ek, x = "Count", showCategory = 15)
ggsave("functional_analysis/cfvshc_down_kegg.png")
