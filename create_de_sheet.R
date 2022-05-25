library("xlsx")

dir_path <- "de_results"
# files <- list.files(path=dir_path, pattern="*.csv")
files <- c("DE_CFVsHC.csv",
           "DE_CFVsHCadultsamples.csv",
           "DE_CFVsHCchildsamples.csv",
           "DE_CFRDVsIGT.csv",
           "DE_CFRDVsNGT.csv",
           "DE_IGTVsNGT.csv",
           "DE_preVspost.csv")

for(f in files){
  sheet_name <- gsub("DE_|.csv", "", f)
  print(sheet_name)
  data <- read.csv(paste(dir_path, f, sep = "/"))
  # print(head(data))
  write.xlsx(data, paste(dir_path, "all_de_results.xlsx", sep = "/"), sheetName = sheet_name,
             col.names = TRUE, row.names = FALSE, append = TRUE)
}
