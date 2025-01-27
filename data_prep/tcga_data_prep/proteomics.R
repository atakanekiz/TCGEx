#------------------------------------------------------------------------------- DOWNLOAD THE PROTEIN DATA

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)

tcga_projects <- c(
  "TCGA-ACC",  # Adrenocortical carcinoma
  "TCGA-BLCA", # Bladder Urothelial Carcinoma
  "TCGA-BRCA", # Breast invasive carcinoma
  "TCGA-CESC", # Cervical squamous cell carcinoma and endocervical adenocarcinoma
  "TCGA-CHOL", # Cholangiocarcinoma
  "TCGA-COAD", # Colon adenocarcinoma
  "TCGA-DLBC", # Lymphoid Neoplasm Diffuse Large B-cell Lymphoma
  "TCGA-ESCA", # Esophageal carcinoma
  "TCGA-GBM",  # Glioblastoma multiforme
  "TCGA-HNSC", # Head and Neck squamous cell carcinoma
  "TCGA-KICH", # Kidney Chromophobe
  "TCGA-KIRC", # Kidney renal clear cell carcinoma
  "TCGA-KIRP", # Kidney renal papillary cell carcinoma
  # "TCGA-LAML", # Acute Myeloid Leukemia (NA)
  "TCGA-LGG",  # Brain Lower Grade Glioma
  "TCGA-LIHC", # Liver hepatocellular carcinoma
  "TCGA-LUAD", # Lung adenocarcinoma
  "TCGA-LUSC", # Lung squamous cell carcinoma
  "TCGA-MESO", # Mesothelioma
  "TCGA-OV",   # Ovarian serous cystadenocarcinoma
  "TCGA-PAAD", # Pancreatic adenocarcinoma
  "TCGA-PCPG", # Pheochromocytoma and Paraganglioma
  "TCGA-PRAD", # Prostate adenocarcinoma
  "TCGA-READ", # Rectum adenocarcinoma
  "TCGA-SARC", # Sarcoma
  "TCGA-SKCM", # Skin Cutaneous Melanoma
  "TCGA-STAD", # Stomach adenocarcinoma
  "TCGA-TGCT", # Testicular Germ Cell Tumors
  "TCGA-THCA", # Thyroid carcinoma
  "TCGA-THYM", # Thymoma
  "TCGA-UCEC", # Uterine Corpus Endometrial Carcinoma
  "TCGA-UCS",  # Uterine Carcinosarcoma
  "TCGA-UVM"   # Uveal Melanoma
)

for (project in tcga_projects) {
  cat("Processing project:", project, "\n")
  
  # query
  query_protein <- GDCquery(
    project = project,
    data.category = "Proteome Profiling",
    data.type = "Protein Expression Quantification",
    experimental.strategy = "Reverse Phase Protein Array"
  )
  
  # download
  tryCatch({
    GDCdownload(query_protein)
    cat("Download completed for project:", project, "\n")
  }, error = function(e) {
    cat("Error in downloading project:", project, "\n", e$message, "\n")
  })
}

#------------------------------------------------------------------------------- PREPARE AND MERGE THE PROTEIN DATA

library(ggpubr)
library(dplyr)
library(tidyr)

#create direction

output_dir <- "tcga_proteomics_data"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# query

for (project in tcga_projects) {
  cat("Processing project:", project, "\n")
  
  
  query_protein <- GDCquery(
    project = project,
    data.category = "Proteome Profiling",
    data.type = "Protein Expression Quantification",
    experimental.strategy = "Reverse Phase Protein Array"
  )
  
  tryCatch({
    # data prepare
    proteomics_data <- GDCprepare(query_protein)
    
    # data name
    file_name <- file.path(output_dir, paste0(project, "_proteomics_data.rds"))
    
    # save data
    saveRDS(proteomics_data, file = file_name)
    cat("Data saved for project:", project, "at", file_name, "\n")
  }, error = function(e) {
    cat("Error processing project:", project, "\n", e$message, "\n")
  })
}


# Prepare the proteomics data into a usable format

  #Find the duplicates (set_id differences)

    data_dir <- "/tcga_proteomics_data"

    rds_files <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE)

        # list the result
        duplicate_summary <- list()

      # loop

for (file in rds_files) {
  cat("Processing file:", file, "\n")
  tryCatch({
    data <- readRDS(file)
    
    # find the duplicates
    if ("peptide_target" %in% colnames(data)) {
      duplicates <- data %>%
        group_by(peptide_target) %>%
        filter(n() > 1) %>%
        ungroup()
      
      # summary
      if (nrow(duplicates) > 0) {
        duplicate_summary[[file]] <- duplicates
        cat("Duplicates found in:", file, " - Rows:", nrow(duplicates), "\n")
      } else {
        cat("No duplicates found in:", file, "\n")
      }
    } else {
      cat("Column `peptide_target` not found in:", file, "\n")
    }
  }, error = function(e) {
    cat("Error processing file:", file, "\n", e$message, "\n")
  })
}

#set_id appendix
        
       data <- readRDS("tcga_proteomics_data/TCGA-####_proteomics_data.rds")        
      
         data <- data %>%
          mutate(peptide_target = ifelse(
            duplicated(peptide_target) | duplicated(peptide_target, fromLast = TRUE),
            paste0(peptide_target, ".", set_id),
            peptide_target
          ))        

        saveRDS(data, "TCGA-UCEC_proteomics_data.rds")




# proteomics_data <- proteomics_data[rowMeans(is.na(proteomics_data)) <= 0.75, ]


proteomics_t <- data %>% 
  select(6:ncol(data)) %>% 
  t() %>% 
  as.data.frame() %>% 
  setNames(make.names(paste0("prt.", data$peptide_target))) %>% 
  tibble::rownames_to_column(var="meta.sample")


----------------------------------------------------------------------------------cancer <- readRDS(".rds")

# grep("sample",colnames(cancer), ignore.case = T, value = T)

# cancer$meta.sample[1:10]

merged <- left_join(cancer, proteomics_t, by=c("meta.sample"))

grep("akt",colnames(merged), ignore.case = T, value = T)

ggscatter(merged, "AKT2", "prt.Akt2")

saveRDS(merged, "UCEC-cpm_TCGA.rds")




#get rid of the duplicated samples - THCA

duplicated_samples <- query_protein[[1]][[1]][duplicated(query_protein[[1]][[1]]$cases), ]

print(duplicated_samples)


query_protein[[1]][[1]] <- query_protein[[1]][[1]][!duplicated(query_protein[[1]][[1]]$cases), ]

proteomics_data <- GDCprepare(query_protein)






