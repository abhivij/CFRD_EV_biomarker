library(tidyverse)
library(xlsx)

#write best biomarkers - combined + combat into excel file

best_features <- read.csv("data/selected_features/best_features_with_is_best.csv") %>%
  filter(is_best == 1) %>%
  filter(grepl(pattern = "CF_EV_tra_combat_", x = dataset_id))

# all_protein_names <- read.csv("Data/Protein/formatted_data/all_protein_names.csv")

file_name <- "best_features_tra_prot.xlsx"


for(i in c(1:nrow(best_features))){
  dataset_id <- best_features[i, "dataset_id"]
  biomarkers <- strsplit(best_features[i, "biomarkers"], split = "|", fixed = TRUE)[[1]]
  
  biomarkers_df <- data.frame(biomarkers = biomarkers)
  if(grepl("tra", dataset_id)){
    #replace _ in mirna names to - : that's how it is represnted in mirbase
    #but pirnabank has _ in names
    biomarkers_df <- biomarkers_df %>%
      mutate(biomarkers = ifelse(!grepl("piR", biomarkers), 
                                 gsub(".", "-", biomarkers, fixed = TRUE),
                                 biomarkers))
    sheet_name <- paste0("tra_", sub("CF_EV_tra_combat_", "", dataset_id))
  } else{
    #this case is proteomics
    #add protein name
    biomarkers_df <- biomarkers_df %>%
      left_join(all_protein_names, by = c("biomarkers" = "from_id")) %>%
      relocate(primary_gene_id, .before = protein_name)
    sheet_name <- paste0("tra_", sub("CF_EV_prot_combat_", "", dataset_id))
  }
  
  #print(biomarkers_df)
  # 
  print(dim(biomarkers_df))
  write.xlsx(biomarkers_df, file = paste0("data/selected_features/", file_name), 
             append = TRUE, sheetName = sheet_name, row.names = FALSE)
}
