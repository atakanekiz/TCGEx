library(dplyr)

tcga_projects<-c("ACC","BLCA","CESC","CHOL","COAD",
                 "DLBC","ESCA","GBM","KICH","KIRP","LAML",
                 "LIHC","MESO","PAAD","PCPG","READ","TGCT",
                 "THYM","UCS","UVM","UCEC","THCA","STAD","SKCM","PRAD",
                 "OV","LUSC","LUAD","LGG","KIRC","HNSC","BRCA","SARC")

for (i in tcga_projects) {
  data <- readRDS(paste0("projects/", i, ".rds"))
  
  # Check if the 'meta.TOTAL.MUTATIONS' column exists in the data
  if ("meta.TOTAL.MUTATIONS" %in% colnames(data)) {
    meta_total_mutations <- as.numeric(data$meta.TOTAL.MUTATIONS)
    data$meta.TOTAL.MUTATIONS <- meta_total_mutations
    saveRDS(data, paste0("projects2/", i, ".rds"))
  } else {
    # The 'meta.TOTAL.MUTATIONS' column does not exist in the data,
    # so move on to the next tcga_project
    next
  }
}
