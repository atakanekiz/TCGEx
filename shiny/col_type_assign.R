library(dplyr)

tcga_projects <- c("ACC","BLCA","CESC","CHOL","COAD",
                   "DLBC","ESCA","GBM","KICH","KIRP","LAML",
                   "LIHC","MESO","PAAD","PCPG","READ","TGCT",
                   "THYM","UCS","UVM","UCEC","THCA","STAD","SKCM","PRAD",
                   "OV","LUSC","LUAD","LGG","KIRC","HNSC","BRCA","SARC")

for (i in tcga_projects) {
  input_file <- paste0("projects/", i, ".rds")
  output_file <- paste0("na_free_and_correct_type_of_column_data/", i, ".rds")
  
  if (!file.exists(input_file)) {
    print(paste("Giriş dosyası bulunamadı:", input_file))
    next
  }
  
  data <- readRDS(input_file)
  
  data <- data %>%
    mutate(across(starts_with("meta."), ~replace(., tolower(.) %in% c("na", "-", "not available", "notavailable","NA","Not Available","not Available"), NA)))
  
  if (!dir.exists("na_free_and_correct_type_of_column_data")) {
    dir.create("na_free_and_correct_type_of_column_data")
  }
  
  saveRDS(data, output_file)
}
