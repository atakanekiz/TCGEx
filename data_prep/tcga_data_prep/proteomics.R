if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)

query_protein<- GDCquery(project="TCGA-###",
                         data.category="Proteome Profiling",
                         data.type="Protein Expression Quantification",
                         experimental.strategy='Reverse Phase Protein Array'
) 

GDCdownload(query_protein)

#-------------------------------------------------------------------------------
library(ggpubr)
library(dplyr)
library(tidyr)

# Prepare the proteomics data into a usable format
proteomics_data <- GDCprepare(query_protein)

anyDuplicated(proteomics_data$peptide_target)

proteomics_data <- proteomics_data[rowMeans(is.na(proteomics_data)) <= 0.75, ]

# peptide_target make unique

proteomics_data$peptide_target <- make.unique(as.character(proteomics_data$peptide_target), sep = ".")


proteomics_t <- proteomics_data %>% 
  select(6:ncol(proteomics_data)) %>% 
  t() %>% 
  as.data.frame() %>% 
  setNames(make.names(paste0("prt.", proteomics_data$peptide_target))) %>% 
  tibble::rownames_to_column(var="meta.sample")


----------------------------------------------------------------------------------cancer <- readRDS(".rds")

# grep("sample",colnames(cancer), ignore.case = T, value = T)

# cancer$meta.sample[1:10]

merged <- left_join(cancer, proteomics_t, by=c("meta.sample"))

grep("akt",colnames(merged), ignore.case = T, value = T)

ggscatter(merged, "AKT2", "prt.Akt2")

saveRDS(merged, "###_prt.rds")

















# 
# 
# duplicated_samples <- query_protein[[1]][[1]][duplicated(query_protein[[1]][[1]]$cases), ]
# 
# print(duplicated_samples)
# 
# 
# query_protein[[1]][[1]] <- query_protein[[1]][[1]][!duplicated(query_protein[[1]][[1]]$cases), ]
# 
# proteomics_data <- GDCprepare(query_protein)
# 
# 




