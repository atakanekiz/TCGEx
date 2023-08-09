library(dplyr)
library(data.table)

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
    mutate(across(starts_with("meta."), ~replace(., tolower(.) %in% c("not-reported","na", "-","not reported","Not Reported","Not Reported", "not available", "notavailable","NA","Not Available","not Available"), NA)))
  
  if (!dir.exists("na_free_and_correct_type_of_column_data")) {
    dir.create("na_free_and_correct_type_of_column_data")
  }
  
  
  # "meta." ile başlayan sütunları bulalım
  meta_columns <- grep("^meta\\.", names(data))
  
  # Her bir "meta." ile başlayan sütunu işleyelim
  for (col in meta_columns) {
    # Sütundaki unique değerleri ve unique değer sayısını alalım
    unique_values <- unique(data[[col]])
    num_unique_values <- length(unique_values)
    numeric_unique_values<- sum(is.numeric(data[[meta.ploidy]]))
    
    # Sütundaki unique değerler sayısal ise ve sayısı 15'ten büyükse "factor" yerine "numeric" yapalım
    if (numeric_unique_values >= 15 && num_unique_values >= 15) {
      data[[col]] <- as.numeric(data[[col]])
    }
  }
  
  saveRDS(data, output_file)
}

