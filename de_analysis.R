library(tidyverse)
library(edgeR)
library(EnhancedVolcano)
library("ggvenn")

base_dir <- "/home/abhivij/UNSW/VafaeeLab/CysticFibrosisGroup/ExoCF/CFRD_EV_biomarker/"
setwd(base_dir)

source("utils.R")


#read transcriptomics data (miRNA)
umi_counts <- read.csv("data/formatted/umi_counts.csv", row.names = 1)
colnames(umi_counts) <- gsub("^X", "", colnames(umi_counts))


#read metadata
meta_data <- read.csv("data/formatted/meta_data.csv")


sum(is.na(umi_counts))



#CF Vs HC all samples
meta_data_subset <- meta_data %>%
  filter(is.na(pre_post_modulator) | pre_post_modulator != 1) %>%
  mutate(de_condition = case_when(condition == "HC" ~ "HC",
                                  TRUE ~ "CF")) %>%
  arrange(condition)

counts <- umi_counts[, meta_data_subset$sample_long_name]
group <- factor(meta_data_subset$de_condition, levels = c("HC", "CF"))


de_analysis(
  counts = counts,
  group = group,
  results_dir_path = "de_results",
  plot_title = "CF Vs HC",
  p_val_cutoff = 0.05,
  logFC_cutoff = 0.5,
  k = 5
)

# counts
# group
# results_dir_path = "de_results"
# plot_title = "CF Vs HC"
# p_val_cutoff = 0.05
# logFC_cutoff = 0.5
# k = 5
de_analysis <- function(counts, group, results_dir_path,
                        plot_title, p_val_cutoff, logFC_cutoff, contrast, k = 5){
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
  
  # contrast = c(0,-1,0,0)
  qlf <- glmQLFTest(fit, contrast=contrast)
  de_rnas <- topTags(qlf, n = Inf)$table
  # 
  # contrast = c(0,-1,0,0)
  # qlf <- glmQLFTest(fit, contrast=contrast)
  # de_rnas2 <- topTags(qlf, n = Inf)$table
  # 
  # qlf <- glmQLFTest(fit, coef=2)
  # de_rnas3 <- topTags(qlf, n = Inf)$table
  
  de_rnas <- de_rnas %>%
    rownames_to_column("Molecule")
  
  p_val_column = "FDR"
  create_volcano_plot(
    result = de_rnas,
    title = plot_title,
    file_name = paste0("volcano_", gsub(" ", "", plot_title, fixed = TRUE),
                       "_", p_val_column, ".png"),
    dir_path = results_dir_path,
    p_val_cutoff = p_val_cutoff, 
    logFC_cutoff = logFC_cutoff,
    p_val_column = p_val_column,
    k = k
  )
  p_val_column = "PValue"
  create_volcano_plot(
    result = de_rnas,
    title = plot_title,
    file_name = paste0("volcano_", gsub(" ", "", plot_title, fixed = TRUE),
                       "_", p_val_column, ".png"),
    dir_path = results_dir_path,
    p_val_cutoff = p_val_cutoff, 
    logFC_cutoff = logFC_cutoff,
    p_val_column = p_val_column,
    k = k
  )
  file_name <- paste0("DE_", gsub(" ", "", plot_title, fixed = TRUE), ".csv")
  significant_de <- de_rnas %>%
    filter(PValue < 0.05) %>%
    arrange(logFC)
  write.csv(significant_de, paste(results_dir_path, file_name, sep = "/"), row.names = FALSE)  
}


###############

#sample design matrix
group <- factor(c(1,1,2,2,3,3))
design <- model.matrix(~group)



##############


#CF Vs HC child samples

meta_data_subset <- meta_data %>%
  filter(is.na(pre_post_modulator) | pre_post_modulator != 1) %>%
  filter(age_group == "child") %>%
  mutate(de_condition = case_when(condition == "HC" ~ "HC",
                                  TRUE ~ "CF")) %>%
  arrange(de_condition)


counts <- umi_counts[, meta_data_subset$sample_long_name]
group <- factor(meta_data_subset$de_condition, levels = c("HC", "CF"))


de_analysis(
  counts = counts,
  group = group,
  results_dir_path = "de_results",
  plot_title = "CF Vs HC child samples",
  p_val_cutoff = 0.05,
  logFC_cutoff = 0.5,
  k = 5
)



#CF Vs HC adult samples
meta_data_subset <- meta_data %>%
  filter(is.na(pre_post_modulator) | pre_post_modulator != 1) %>%
  filter(age_group == "adult") %>%
  mutate(de_condition = case_when(condition == "HC" ~ "HC",
                                  TRUE ~ "CF")) %>%
  arrange(de_condition)

counts <- umi_counts[, meta_data_subset$sample_long_name]
group <- factor(meta_data_subset$de_condition, levels = c("HC", "CF"))


de_analysis(
  counts = counts,
  group = group,
  results_dir_path = "de_results",
  plot_title = "CF Vs HC adult samples",
  p_val_cutoff = 0.05,
  logFC_cutoff = 0.5,
  k = 5
)


de <- read.csv("de_results/DE_CFVsHC.csv")
de_adult <- read.csv("de_results/DE_CFVsHCadultsamples.csv")
de_child <- read.csv("de_results/DE_CFVsHCchildsamples.csv")

de_up <- de %>%
  filter(logFC > 0)
de_down <- de %>%
  filter(logFC < 0)

de_adult_up <- de_adult %>%
  filter(logFC > 0)
de_adult_down <- de_adult %>%
  filter(logFC < 0)

de_child_up <- de_child %>%
  filter(logFC > 0)
de_child_down <- de_child %>%
  filter(logFC < 0)


ggvenn(list("CF Vs HC Up Reg" = de_up$Molecule, 
            "CF Vs HC Adult Up Reg" = de_adult_up$Molecule,
            "CF Vs HC Child Up Reg" = de_child_up$Molecule),
       stroke_size = 0.1,
       set_name_size = 5,
       text_size = 3, )
ggsave("de_results/CFVsHC_up_reg_venn.png")

Reduce(intersect, list(de_up$Molecule, de_adult_up$Molecule, de_child_up$Molecule))


ggvenn(list("CF Vs HC Down Reg" = de_down$Molecule, 
            "CF Vs HC Adult Down Reg" = de_adult_down$Molecule,
            "CF Vs HC Child Down Reg" = de_child_down$Molecule),
       stroke_size = 0.1,
       set_name_size = 5,
       text_size = 3, )
ggsave("de_results/CFVsHC_down_reg_venn.png")

Reduce(intersect, list(de_down$Molecule, de_adult_down$Molecule, de_child_down$Molecule))



############

moi <- read.csv("data/formatted/validated_mirnas_in_data.csv")[,1]
moi <- gsub(".", "-", moi, fixed = TRUE)


ggvenn(list("CF Vs HC Up Reg" = de_up$Molecule, 
            "Validated IGFBP7 miRNAs" = moi),
       stroke_size = 0.1,
       set_name_size = 5,
       text_size = 3, )
ggsave("de_results/CFVsHC_up_validated.png")
intersect(de_up$Molecule, moi)

ggvenn(list("CF Vs HC Down Reg" = de_down$Molecule, 
            "Validated IGFBP7 miRNAs" = moi),
       stroke_size = 0.1,
       set_name_size = 5,
       text_size = 3, )
ggsave("de_results/CFVsHC_down_validated.png")
intersect(de_down$Molecule, moi)

###############

#CFRD Vs IGT

meta_data_subset <- meta_data %>%
  filter(is.na(pre_post_modulator) | pre_post_modulator != 1) %>%
  filter(condition %in% c("CFRD", "IGT", "NGT", "HC")) %>%
  mutate(de_condition = condition) %>%
  arrange(de_condition)

summary(factor(meta_data_subset$de_condition))


counts <- umi_counts[, meta_data_subset$sample_long_name]
group <- factor(meta_data_subset$de_condition, levels = c("CFRD", "IGT", "NGT", "HC"))


# counts = counts
# group = group
# results_dir_path = "de_results"
# plot_title = "CFRD Vs IGT"
# p_val_cutoff = 0.05
# logFC_cutoff = 0.5
# k = 5
# coef = 1:2

de_analysis(
  counts = counts,
  group = group,
  results_dir_path = "de_results",
  plot_title = "CFRD Vs NGT",
  p_val_cutoff = 0.05,
  logFC_cutoff = 0.5,
  k = 5,
  contrast = c(0, -1, 0, 0)
)

# contrast = c(1, -1, 0, 0)


#CFRD Vs NGT
de_analysis(
  counts = counts,
  group = group,
  results_dir_path = "de_results",
  plot_title = "CFRD Vs NGT",
  p_val_cutoff = 0.05,
  logFC_cutoff = 0.5,
  k = 5,
  contrast = c(0, 0, -1, 0)
)
de_analysis(
  counts = counts,
  group = group,
  results_dir_path = "de_results",
  plot_title = "CFRD Vs IGT adult samples",
  p_val_cutoff = 0.05,
  logFC_cutoff = 0.5,
  k = 5,
  contrast = c(0, -1, 0, 0)
)


#IGT Vs NGT
de_analysis(
  counts = counts,
  group = group,
  results_dir_path = "de_results",
  plot_title = "IGT Vs NGT",
  p_val_cutoff = 0.05,
  logFC_cutoff = 0.5,
  k = 5,
  contrast = c(0, 1, -1, 0)
)



de_c_i <- read.csv("de_results/DE_CFRDVsIGT.csv")
de_c_n <- read.csv("de_results/DE_CFRDVsNGT.csv")
de_i_n <- read.csv("de_results/DE_IGTVsNGT.csv")

de_c_i_up <- de_c_i %>%
  filter(logFC > 0)
de_c_i_down <- de_c_i %>%
  filter(logFC < 0)

de_c_n_up <- de_c_n %>%
  filter(logFC > 0)
de_c_n_down <- de_c_n %>%
  filter(logFC < 0)

de_i_n_up <- de_i_n %>%
  filter(logFC > 0)
de_i_n_down <- de_i_n %>%
  filter(logFC < 0)


ggvenn(list("CFRD Vs IGT Up Reg" = de_c_i_up$Molecule, 
            "CFRD Vs NGT Up Reg" = de_c_n_up$Molecule,
            "IGT Vs NGT Up Reg" = de_i_n_up$Molecule),
       stroke_size = 0.1,
       set_name_size = 5,
       text_size = 3, )
ggsave("de_results/CFRDVsIGTVsNGT_up_reg_venn.png")


ggvenn(list("CFRD Vs IGT Down Reg" = de_c_i_down$Molecule, 
            "CFRD Vs NGT Down Reg" = de_c_n_down$Molecule,
            "IGT Vs NGT Down Reg" = de_i_n_down$Molecule),
       stroke_size = 0.1,
       set_name_size = 5,
       text_size = 3, )
ggsave("de_results/CFRDVsIGTVsNGT_down_reg_venn.png")


########################
#modulator vs premodulator

meta_data_subset <- meta_data %>%
  filter(!is.na(pre_post_modulator)) %>%
  mutate(de_condition = factor(pre_post_modulator)) %>%
  arrange(de_condition)
levels(meta_data_subset$de_condition) <- c("pre", "post")

summary(factor(meta_data_subset$de_condition))


counts <- umi_counts[, meta_data_subset$sample_long_name]
group <- factor(meta_data_subset$de_condition, levels = c("pre", "post"))


de_analysis(
  counts = counts,
  group = group,
  results_dir_path = "de_results",
  plot_title = "pre Vs post",
  p_val_cutoff = 0.05,
  logFC_cutoff = 0.5,
  k = 5,
  contrast = c(0, -1)
)
