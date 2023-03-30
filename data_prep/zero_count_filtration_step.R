library(dplyr)

tcga_projects<-c("ACC","BLCA","CESC","CHOL","COAD",
                 "DLBC","ESCA","GBM","KICH","KIRP","LAML",
                 "LIHC","MESO","PAAD","PCPG","READ","TGCT",
                 "THYM","UCS","UVM","UCEC","THCA","STAD","SKCM","PRAD",
                 "OV","LUSC","LUAD","LGG","KIRC","HNSC","BRCA","SARC")

for(i in tcga_projects){
  
  data<-readRDS(paste0("datas/", i, ".rds"))
  
  without_metadata<-data %>% select(!starts_with("meta."))
  
  metadata<-data %>% select(starts_with("meta."))
  
  gene_cols <- unlist(lapply(without_metadata, is.numeric))
  
  # Subset numeric columns of data (genes only)
  
  genes_data <- without_metadata[ , gene_cols]                        
  
  barcodes <- without_metadata[ , !gene_cols]
  
  names(barcodes) <- paste0("meta.", names(barcodes))
  
  # Calculate zero percent of columns
  
  zero_percentage <- apply(genes_data, 2, function(x) sum(x == 0, na.rm=T) / length(x) * 100)
  
  # Detection of columns with zero percent greater than 25%
  
  deleting_cols <- names(genes_data)[zero_percentage > 25]
  
  # Removing the columns to be deleted from the data frame
  
  filtred_gene_data <- genes_data[, !names(genes_data) %in% deleting_cols]
  
  filtered_data<-cbind.data.frame(filtred_gene_data, metadata, barcodes)
  
  saveRDS(filtered_data,paste0("zero_count_filtration/", i, ".rds"))

  
}

 # Turn the data to data.table

 library(data.table)
 library(zstdlite)

for(i in tcga_projects){

 data_as_data_table<-readRDS(paste0("zero_count_filtration/", i, ".rds"))

 data_as_data_table<-data.table(data_as_data_table)

 saveRDS(data_as_data_table,paste0("zero_count_filtration/", i, ".rds"))
 
 # Compression Step
 
 compressed_filtered_data <- zstd_serialize(data_as_data_table, level = 22)
 
 saveRDS(compressed_filtered_data, paste0("zero_count_filtrated_compressed/", i, ".rds"))

 }


