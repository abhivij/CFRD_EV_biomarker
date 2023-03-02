library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(multiMiR)
library(xlsx)
library(ggvenn)

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


mirna_vec = cond1higher
plot_title = "GO processes related to targets of miRNAs higher in CFRD among CFRD Vs IGT biomarkers"
output_dir = "prediction_pipeline/functional_analysis"
output_plot_name = "CFRDVsIGT_CFRD.png"
output_result_file_name = "CFRDVsIGT_CFRD.csv"

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
  
  if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
  }
  barplot(ego.bp, x = "Count", showCategory = 15,
          title = plot_title)
  ggsave(paste0(output_dir, "/", output_plot_name), width = 15, height = 10)
  # heatplot(ego.bp, showCategory = 10)
  write.csv(ego.bp.result, paste0(output_dir, "/", output_result_file_name), row.names = FALSE)  
}


go_BP <- function(input_file_path, output_dir, output_file_name, plot_title){
  de <- read.csv(input_file_path) %>%
    filter(FDR < 0.05)
  
  de.mir <- de %>%
    filter(grepl("miR", Molecule) | grepl("hsa-let", Molecule))
  print("de mir")
  print(dim(de.mir))
  
  de.pir <- de %>%
    filter(grepl("piR", Molecule))
  print("de pir")
  print(dim(de.pir))
  
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
  print("up reg FDR < 0.05")
  print(dim(de.mir.up))
  
  de.mir.down <- de.mir %>%
    filter(logFC < 0) %>%
    arrange(desc(-logFC))
  print("down reg FDR < 0.05")
  print(dim(de.mir.down))
  
  
  #write fdr < 0.05 de mirnas to xlsx file
  dir_path <- strsplit(input_file_path, split = "/")[[1]][1]
  de_file_name <- strsplit(input_file_path, split = "/")[[1]][2] 
  
  dir_path <- paste(dir_path, "sig", sep = "_")
  if(!dir.exists(dir_path)){
    dir.create(dir_path, recursive = TRUE)
  }
  de_file_name <- sub(".csv", ".xlsx", de_file_name)
  
  write.xlsx(de.mir.up, paste(dir_path, de_file_name, sep = "/"), 
             sheetName = "mir_upreg",
             col.names = TRUE, row.names = FALSE, append = TRUE)
  write.xlsx(de.mir.down, paste(dir_path, de_file_name, sep = "/"), 
             sheetName = "mir_downreg",
             col.names = TRUE, row.names = FALSE, append = TRUE)

  
  
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




#############
#functional analsysis of identified biomarkers

#CFRDVsIGT
cond1higher <- c('hsa-miR-192-5p')
cond2higher <- c('hsa-miR-181d-5p')

obtain_mirna_go_bp(mirna_vec = cond1higher, 
                   plot_title = "GO processes related to targets of miRNAs higher in CFRD among CFRD Vs IGT biomarkers",
                   output_dir = "prediction_pipeline/functional_analysis",
                   output_plot_name = "CFRDVsIGT_CFRD.png",
                   output_result_file_name = "CFRDVsIGT_CFRD.csv")

obtain_mirna_go_bp(mirna_vec = cond2higher, 
                   plot_title = "GO processes related to targets of miRNAs higher in IGT among CFRD Vs IGT biomarkers",
                   output_dir = "prediction_pipeline/functional_analysis",
                   output_plot_name = "CFRDVsIGT_IGT.png",
                   output_result_file_name = "CFRDVsIGT_IGT.csv")

#CFRDVsNGT
cond1higher <- c("hsa-miR-122-5p", "hsa-miR-342-3p", 
                 "hsa-miR-486-3p", "hsa-miR-182-5p",
                 "hsa-miR-17-5p", "hsa-miR-4443" 
                 )
cond2higher <- c("hsa-miR-181c-3p", "hsa-miR-181d-5p",
                 "hsa-miR-27b-5p", "hsa-miR-3617-5p",
                 "hsa-miR-4286", "hsa-miR-4497",
                 "hsa-miR-501-5p", "hsa-miR-519d-5p",
                 "hsa-miR-92b-5p")

obtain_mirna_go_bp(mirna_vec = cond1higher, 
                   plot_title = "GO processes related to targets of miRNAs higher in CFRD among CFRD Vs NGT biomarkers",
                   output_dir = "prediction_pipeline/functional_analysis",
                   output_plot_name = "CFRDVsNGT_CFRD.png",
                   output_result_file_name = "CFRDVsNGT_CFRD.csv")

obtain_mirna_go_bp(mirna_vec = cond2higher, 
                   plot_title = "GO processes related to targets of miRNAs higher in NGT among CFRD Vs NGT biomarkers",
                   output_dir = "prediction_pipeline/functional_analysis",
                   output_plot_name = "CFRDVsNGT_NGT.png",
                   output_result_file_name = "CFRDVsNGT_NGT.csv")

#IGTVsNGT
cond1higher <- c("hsa-miR-17-5p")
cond2higher <- c("hsa-miR-3617-5p", "hsa-miR-326")

obtain_mirna_go_bp(mirna_vec = cond1higher, 
                   plot_title = "GO processes related to targets of miRNAs higher in IGT among IGT Vs NGT biomarkers",
                   output_dir = "prediction_pipeline/functional_analysis",
                   output_plot_name = "IGTVsNGT_IGT.png",
                   output_result_file_name = "IGTVsNGT_IGT.csv")

obtain_mirna_go_bp(mirna_vec = cond2higher, 
                   plot_title = "GO processes related to targets of miRNAs higher in NGT among IGT Vs NGT biomarkers",
                   output_dir = "prediction_pipeline/functional_analysis",
                   output_plot_name = "IGTVsNGT_NGT.png",
                   output_result_file_name = "IGTVsNGT_NGT.csv")


mir.targets1 <- get_multimir(
  org     = "hsa",
  mirna  = c("hsa-miR-17-5p"),
  table   = "validated",
  summary = TRUE
)  
mir.targets1.data <- mir.targets1@data[, c(2, 3, 5)] %>%
  unique() %>%
  filter(target_entrez != "")

mir.targets2 <- get_multimir(
  org     = "hsa",
  mirna  = c("hsa-miR-17-3p"),
  table   = "validated",
  summary = TRUE
)  
mir.targets2.data <- mir.targets2@data[, c(2, 3, 5)] %>%
  unique() %>%
  filter(target_entrez != "")


mir.targets3 <- get_multimir(
  org     = "hsa",
  mirna  = c("hsa-miR-17"),
  table   = "validated",
  summary = TRUE
)  
mir.targets3.data <- mir.targets3@data[, c(2, 3, 5)] %>%
  unique() %>%
  filter(target_entrez != "")



##############
#compare pancreas-specific
pancreas_all <- read_tsv("data/pancreas_all.tsv")
pancreas_enriched <- read_tsv("data/pancreas_enriched.tsv")
pancreas_specific <- read_tsv("data/pancreas_specific.tsv")

ggvenn(list("Pancreas all" = pancreas_all$Gene,
            "Pancreas enriched" = pancreas_enriched$Gene,
            "Pancreas specific" = pancreas_specific$Gene),
       stroke_size = 0.1, fill_color = c("skyblue3", "orange3", "yellow3"),
       set_name_size = 4,
       text_size = 3, stroke_linetype = "blank")
ggsave("plots/venn/pancreas_ref_data.png")


#################


mirna_vec = cond1higher
plot_title = "GO processes related to targets of miRNAs higher in CFRD among CFRD Vs IGT biomarkers"
output_dir = "prediction_pipeline/functional_analysis"
output_plot_name = "CFRDVsIGT_CFRD.png"
output_result_file_name = "CFRDVsIGT_CFRD.csv"

comparison <- "CFRDVsIGT"
output_dir_path = "prediction_pipeline/mirna_targets"


#mirna_vec : a single element vector
get_mirna_targets <- function(mirna_vec,
                              comparison,
                              output_dir_path = "prediction_pipeline/mirna_targets"){
  
  mir.targets <- get_multimir(
    org     = "hsa",
    mirna  = mirna_vec,
    table   = "validated",
    summary = TRUE
  )  
  
  mir.targets.data <- mir.targets@data %>%
    arrange(target_entrez, target_ensembl)
  
  ############
  #using only ones with ensemblid 
  #merging based on target_symbol and gene, merges haplotypic ones
  
  # dim(mir.targets.data %>%
  #   filter(target_ensembl == ""))
  # #164
  # dim(mir.targets.data %>%
  #       filter(target_entrez == ""))
  # #191
  # 
  # therefore use ensembl
  
  # doi <- mir.targets.data %>%
  #   filter(target_entrez != "") %>%
  #   inner_join(pancreas_enriched, by = c("target_ensembl" = "Ensembl"))
  # doi.other <- mir.targets.data %>%
  #   # filter(target_entrez == "") %>%
  #   inner_join(pancreas_enriched, by = c("target_symbol" = "Gene"))  %>%
  #   filter(target_ensembl != Ensembl)
  # 
  # doi %>%
  #   filter(target_ensembl %in% doi.other$Ensembl)
  # 
  # mir.targets.data.group_count <- mir.targets.data %>%
  #   group_by(target_entrez) %>%
  #   summarize(count = n()) %>%
  #   arrange(desc(count))
  
  ############
  
  # #as seen below in commented code, some targets are present in multiple dbs, 
  # #and sometimes in each db multiple times, for different experiments
  # #so just use unique target_ensembl
  # 
  # mir.targets.data.processed <- mir.targets.data %>%
  #   filter(target_ensembl != "")
  # multiple <- mir.targets.data.processed %>%
  #   group_by(target_ensembl) %>%
  #   summarize(count = n()) %>%
  #   arrange(desc(count))
  # doi <- mir.targets.data.processed %>%
  #   filter(target_ensembl == 'ENSG00000139687')
  # #6 entries
  
  #the below is zero for targets from "hsa-miR-192-5p"
  #assuming this for others and using only those with target_ensembl id present.
  
  # doi.other <- mir.targets.data %>%
  #   filter(target_ensembl == "") %>%
  #   inner_join(pancreas_enriched, by = c("target_symbol" = "Gene"))
  
  unique_targets <- mir.targets.data %>%
    filter(target_ensembl != "") %>%
    distinct(target_ensembl)
  
  #read targets of interest from reference
  pancreas_all <- read_tsv("data/pancreas_all.tsv")
  pancreas_enriched <- read_tsv("data/pancreas_enriched.tsv")
  pancreas_specific <- read_tsv("data/pancreas_specific.tsv")
  ###################
  targets_of_interest.all <- pancreas_all %>%
    dplyr::filter(Ensembl %in% unique_targets$target_ensembl) %>%
    dplyr::select(c(1:11, 16:19))
    
  targets_of_interest.enriched <- pancreas_enriched %>%
    dplyr::filter(Ensembl %in% unique_targets$target_ensembl) %>%
    dplyr::select(c(1:11, 16:19))
  
  targets_of_interest.specific <- pancreas_specific %>%
    dplyr::filter(Ensembl %in% unique_targets$target_ensembl) %>%
    dplyr::select(c(1:11, 16:19))
  
  if(!exists(output_dir_path)){
    dir.create(output_dir_path, recursive = TRUE)
  }
  
  #assuming mirna_vec has only 1 element and extracting first
  write.xlsx(as.data.frame(targets_of_interest.all), 
             paste0(output_dir_path, "/", comparison, "_", mirna_vec[1], ".xlsx"),
             sheetName = "pancreas all", showNA = FALSE,
             col.names = TRUE, row.names = FALSE)
  write.xlsx(as.data.frame(targets_of_interest.enriched), 
             paste0(output_dir_path, "/", comparison, "_", mirna_vec[1], ".xlsx"),
             sheetName = "pancreas enriched", showNA = FALSE,
             col.names = TRUE, row.names = FALSE, append = TRUE)
  write.xlsx(as.data.frame(targets_of_interest.specific), 
             paste0(output_dir_path, "/", comparison, "_", mirna_vec[1], ".xlsx"),
             sheetName = "pancreas specifc", showNA = FALSE,
             col.names = TRUE, row.names = FALSE, append = TRUE)
  
  summary_row <- data.frame("mirna" = mirna_vec[1],
                            "comparison" = comparison,
                            "pancreas specific targets" = nrow(targets_of_interest.specific),
                            "pancreas enriched targets" = nrow(targets_of_interest.enriched),
                            "pancreas expressed targets" = nrow(targets_of_interest.all))
  write.table(summary_row,
              paste0(output_dir_path, "/pancreatic_targets_summary.csv"),
              row.names = FALSE, append = TRUE, sep = ",",
              col.names = !file.exists(paste0(output_dir_path, "/pancreatic_targets_summary.csv")))
}


get_mirna_targets(mirna_vec = c("hsa-miR-181d-5p"),
                  comparison = "CFRDVsIGT")
get_mirna_targets(mirna_vec = c("hsa-miR-192-5p"),
                  comparison = "CFRDVsIGT")

get_mirna_targets(mirna_vec = c("hsa-miR-122-5p"),
                  comparison = "CFRDVsNGT")
get_mirna_targets(mirna_vec = c("hsa-miR-342-3p"),
                  comparison = "CFRDVsNGT")
get_mirna_targets(mirna_vec = c("hsa-miR-486-3p"),
                  comparison = "CFRDVsNGT")
get_mirna_targets(mirna_vec = c("hsa-miR-182-5p"),
                  comparison = "CFRDVsNGT")
get_mirna_targets(mirna_vec = c("hsa-miR-17-5p"),
                  comparison = "CFRDVsNGT")
get_mirna_targets(mirna_vec = c("hsa-miR-181d-5p"),
                  comparison = "CFRDVsNGT")
get_mirna_targets(mirna_vec = c("hsa-miR-27b-5p"),
                  comparison = "CFRDVsNGT")
get_mirna_targets(mirna_vec = c("hsa-miR-181c-3p"),
                  comparison = "CFRDVsNGT")
get_mirna_targets(mirna_vec = c("hsa-miR-4443"),
                  comparison = "CFRDVsNGT")
get_mirna_targets(mirna_vec = c("hsa-miR-92b-5p"),
                  comparison = "CFRDVsNGT")
get_mirna_targets(mirna_vec = c("hsa-miR-3617-5p"),
                  comparison = "CFRDVsNGT")
get_mirna_targets(mirna_vec = c("hsa-miR-501-5p"),
                  comparison = "CFRDVsNGT")
get_mirna_targets(mirna_vec = c("hsa-miR-4286"),
                  comparison = "CFRDVsNGT")
get_mirna_targets(mirna_vec = c("hsa-miR-4497"),
                  comparison = "CFRDVsNGT")
get_mirna_targets(mirna_vec = c("hsa-miR-519d-5p"),
                  comparison = "CFRDVsNGT")


get_mirna_targets(mirna_vec = c("hsa-miR-17-5p"),
                  comparison = "IGTVsNGT")
get_mirna_targets(mirna_vec = c("hsa-miR-3617-5p"),
                  comparison = "IGTVsNGT")
get_mirna_targets(mirna_vec = c("hsa-miR-326"),
                  comparison = "IGTVsNGT")


#sort summary

summary <- read.csv("prediction_pipeline/mirna_targets/pancreatic_targets_summary.csv")
summary_sorted <- summary %>%
  arrange(comparison, 
          desc(pancreas.specific.targets),
          desc(pancreas.enriched.targets),
          desc(pancreas.expressed.targets))
write.csv(summary_sorted, "prediction_pipeline/mirna_targets/pancreatic_targets_summary_sorted.csv",
          row.names = FALSE)
