library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(tidyverse)
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

de <- read.csv("de_results/DE_CFVsHC.csv")

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
de.mir.down <- de.mir %>%
  filter(logFC < 0)

de.mir.down.targets <- get_multimir(
  org     = "hsa",
  mirna  = de.mir.down$Molecule,
  table   = "validated",
  summary = TRUE
)
de.mir.down.targets@summary

de.mir.down.targets.data <- de.mir.down.targets@data[, c(2:6)] %>%
  unique() %>%
  filter(target_entrez != "")
gene_list <- unique(de.mir.down.targets.data$target_entrez)


ego <- enrichGO(gene          = gene_list,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "fdr",
                pvalueCutoff  = 0.05)
tmp <- data.frame(ego)
tmp2 <- ego@result

all.equal(tmp, tmp2[1:2168, ])
