library(tidyverse)
library(edgeR)
library(EnhancedVolcano)

base_dir <- "/home/abhivij/UNSW/VafaeeLab/CysticFibrosisGroup/ExoCF/CFRD_EV_biomarker/"
setwd(base_dir)


#read transcriptomics data (miRNA)
umi_counts <- read.csv("data/formatted/umi_counts.csv", row.names = 1)
colnames(umi_counts) <- gsub("^X", "", colnames(umi_counts))


#read metadata
meta_data <- read.csv("data/formatted/meta_data.csv")


sum(is.na(umi_counts))



#CF Vs HC all samples
meta_data_subset <- meta_data %>%
  mutate(de_condition = case_when(!is.na(pre_post_modulator) ~ "modulator_effect",
                                  condition == "CF_pre_post_modulator" ~ "modulator_effect",
                                  condition == "HC" ~ "HC",
                                  TRUE ~ "CF")) %>%
  arrange(de_condition)

counts <- umi_counts[, meta_data_subset$sample_long_name]
group <- factor(meta_data_subset$de_condition, levels = c("HC", "CF", "modulator_effect"))

y <- DGEList(counts = counts, 
             group = group)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y
# head(y$counts)

y <- calcNormFactors(y)
y

design <- model.matrix(~group)
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design)

qlf.2vs1 <- glmQLFTest(fit, coef=2)
de_rnas <- topTags(qlf.2vs1, n = Inf)$table
EnhancedVolcano(de_rnas,
                lab = row.names(de_rnas),
                x = 'logFC',
                y = 'FDR', pCutoff = 0.05,
                title = "CF Vs HC",
                labSize = 3,
                col = c("grey", "grey", "blue", "red"),
                colAlpha = 1)

de_rnas <- de_rnas %>%
 rownames_to_column("Molecule")

result = de_rnas
title = "CF Vs HC"
file_name = "volcano_CFVsHC.png"
dir_path = "plots"
p_val_cutoff = 0.05 
logFC_cutoff = 5

p_val_column = "FDR"
k = 5

###############

#sample design matrix
group <- factor(c(1,1,2,2,3,3))
design <- model.matrix(~group)
