library(data.table)
library(openxlsx)
library(readxl)

rna<-fread("data_mrna_seq_fpkm.txt")

# rna<-fread("data_RNA_Seq_expression_tpm_all_sample_Zscores.txt")

clin_pat<-fread("data_clinical_patient.txt")

clin_samp<-fread("data_clinical_sample.txt")

clin <- merge(clin_pat, clin_samp, by = "#Patient Identifier", all = TRUE)

clin<-clin[4:41 , 1:34]

rm(clin_pat,clin_samp)

colnames(rna)[1] <- "gene"


suppressPackageStartupMessages(library(org.Hs.eg.db))

hs <- org.Hs.eg.db

conv_names <- AnnotationDbi::select(hs,
                                    keys=as.character(rna$gene),
                                    columns = c("ENTREZID", "SYMBOL"),
                                    keytype = "ENTREZID")


sum(duplicated(conv_names$SYMBOL))

rna$gene <- conv_names$SYMBOL

rna <- rna[!is.na(rna$gene), ]


anyDuplicated(rna$gene)

duplicate_genes <- rna[duplicated(gene)]

duplicate_indices <- rna[duplicated(gene), which = TRUE]

rna[, gene := make.unique(gene, sep = ".")]

anyDuplicated(rna$gene)

rm(duplicate_genes,duplicate_indices)


rna_res<-as.data.table(t(rna))

colnames(rna_res) <- as.character(rna_res[1, ])

rna_res <- rna_res[-1, ]

rna$gene<-NULL

rna_res$`Sample Identifier` <-colnames(rna)

rm(rna)

clin_names <- colnames(clin)
clin_names[clin_names != "Sample Identifier"] <- paste0("meta.", clin_names[clin_names != "Sample Identifier"])
colnames(clin) <- clin_names


total <- cbind(rna_res, clin[match(rna_res$`Sample Identifier`, clin$`Sample Identifier`), -1])

rm(clin_names)

colnames(total)[colnames(total) == "Sample Identifier"] <- "meta.Sample Identifier"

for (col in colnames(total)) {
  if (!startsWith(col, "meta.")) {
    total[[col]] <- as.numeric(total[[col]])
  } else {
    total[[col]] <- as.factor(total[[col]])
  }
}

rm(rna_res,col)

colnames(total)[colnames(total) == "meta.Overall Survival (Months)"] <- "meta.days_to_event"

colnames(total)[colnames(total) == "meta.Overall Survival Status"] <- "meta.vital_status"

colnames(total)[colnames(total) == "meta.Sex"] <- "meta.gender"

#column class is changed to numeric if its originally numeric

meta_cols <- grep("^meta\\.", names(total))

for (col in meta_cols) {
  if (is.factor(total[[col]])) {
    numeric_values <- as.numeric(levels(total[[col]])[total[[col]]])
    if (all(!is.na(numeric_values))) {
      total[[col]] <- numeric_values
    }
  } else if (is.character(total[[col]])) {
    numeric_values <- as.numeric(total[[col]])
    if (all(!is.na(numeric_values))) {
      total[[col]] <- numeric_values
    }
  }
}

rm(col,meta_cols,numeric_values)

total$meta.vital_status <- ifelse(total$meta.vital_status == "0:LIVING", 0, 
                                  ifelse(total$meta.vital_status == "1:DECEASED", 1, total$meta.vital_status))


# total$meta.days_to_event <- as.numeric(levels(total$meta.days_to_event)[total$meta.days_to_event])

total$meta.days_to_event <- as.numeric(total$meta.days_to_event) * 30


# duplicated_columns <- total[, duplicated(names(total))]
# 
# anyDuplicated(colnames(total))
# 
# for (col in names(total)[duplicated_columns]) {
#   avg1 <- mean(total[[col]], na.rm = TRUE)
#   avg2 <- mean(total[[paste0(col, ".1")]], na.rm = TRUE)  # dup cols
#   if (!is.na(avg1) && !is.na(avg2) && avg1 > avg2) {
#     total[[paste0(col, ".1")]] <- NULL  
#   } else {
#     total[[col]] <- NULL  
#   }
# }

for (col in colnames(total)) {
  if (!startsWith(col, "meta.")) {
    total[[col]] <- as.numeric(total[[col]])
  }
}

# total$`meta.LDH (treatment Start)` <- as.numeric(levels(total$`meta.LDH (treatment Start)`)[total$`meta.LDH (treatment Start)`])
# total$`meta.Cycles On Therapy` <- as.numeric(levels(total$`meta.Cycles On Therapy`)[total$`meta.Cycles On Therapy`])

total$meta.patient<-total$`meta.Sample Identifier`
    
total[, meta.project_id := as.factor(sapply(1:.N, function(x) paste(sample(c(letters, 0:9), 7, replace = TRUE), collapse = "")))]

# total$meta.definition<-total$`meta.Best Radiographic Response (RECIST 1.1)`

colnames(total)[23439] <- "meta.Sample Identifier.1"

# total$meta.sample_type<-total$meta.definition

# total$meta.therapy <- total$meta.Immunotherapy

# total$meta.therapy <- gsub("MK 3475|MK3475", "MK3475 Therapy", total$meta.therapy)
# total$meta.therapy <- gsub("Nivo|nivolumab|Nivolumab", "Nivolumab Therapy", total$meta.therapy)
# total$meta.therapy <- gsub("Pembro|pembrolizumab|Pembrolizumab", "Pembrolizumab Therapy", total$meta.therapy)
# 
# total$meta.therapy <- as.factor(total$meta.therapy)
# 
# total$meta.definition<- NULL

setnames(total, gsub(" ", ".", names(total)))

setnames(total, gsub("\\(", ".", names(total)))

setnames(total, gsub(")", ".", names(total)))

# total$meta.ICB <- total$meta.therapy
# 
# total$meta.ICB <- gsub("MK3475 Therapy", "Pembrolizumab Therapy", total$meta.ICB)
# 
# total$meta.ICB<-as.factor(total$meta.ICB)


load("C:/Users/memre/OneDrive - Izmir Yuksek Teknoloji Enstitusu/IYTE/Masters Thesis/R/projects/TCGEx/shiny/genesets/msigdb_hallmark.rda")

ifng_genes<- msigdb_hallmark[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]

ifng_genes <- intersect(ifng_genes, colnames(total))

total[, meta_msigdb_ifng_score := rowMeans(.SD, na.rm = TRUE), .SDcols = ifng_genes]


ayers_ifng_genes <- c("CD3D", "IDO1", "CIITA", "CD3E","CCL5","GZMK","CD2","HLA-DRA","CXCL13","IL2RG","NKG7","HLA-E","CXCR6","LAG3","TAGAP","CXCL10","STAT1","GZMB")

ayers_ifng_genes <- intersect(ayers_ifng_genes, colnames(total))

total[, meta_ayers_ifng_score := rowMeans(.SD, na.rm = TRUE), .SDcols = ayers_ifng_genes]

# total[, meta_categorized_response := ifelse(total$meta.sample_type %in% c("Complete Response", "Partial Response"), 1, 
#                                             ifelse(total$meta.sample_type == "Progressive Disease", 0, NA))]

saveRDS(total, "mel_tsam_liang_2017_FPKM_all.rds")

benefit<-fread("data_timeline_treatment.txt")

clin$PATIENT_ID <- clin$`meta.#Patient Identifier`

result <- merge(benefit, clin[, c("PATIENT_ID", "Sample Identifier")], by = "PATIENT_ID", all.x = TRUE)

result$meta.patient<-result$`Sample Identifier`

result$meta.therapy<-result$AGENT

result$meta.therapy_response<-result$TREATMENT_BEST_RESPONSE




# saveRDS(total, "mel_dfci_2019_zTPM_all.rds")


#################################################################################################################################################################
uploaded_data <- readRDS("C:/Users/memre/OneDrive - Izmir Yuksek Teknoloji Enstitusu/IYTE/PhD/R/CbioPortal/mel_dfci_2019/msk_2014_RPKM_all.rds")

df_num = uploaded_data %>% select(where(is.numeric))

df_nonnum = uploaded_data %>% select(-where(is.numeric))

df_nongene = df_num %>% 
  select(starts_with("meta."),starts_with("meta_"))

df_gene <- df_num %>%
  select(-starts_with("meta."), -starts_with("meta_"))

df_gene <- df_gene[, apply(df_gene,2,var)!=0, with=F]

df_gene = log(df_gene +1, base=10)

uploaded_data = cbind(df_gene, df_nonnum,df_nongene)

saveRDS(uploaded_data, "mel_dfci_2019_TPM_LOG_all.rds")

nivo_data <- uploaded_data[meta.ICB == "Nivolumab Therapy", ]

saveRDS(nivo_data, "mel_dfci_2019_TPM_LOG_nivo.rds")

pembro_data <- uploaded_data[meta.ICB == "Pembrolizumab Therapy", ]

saveRDS(pembro_data, "mel_dfci_2019_TPM_LOG_pembro.rds")


################################################################################################################################################################
##"meta_" log


load("C:/Users/memre/OneDrive - Izmir Yuksek Teknoloji Enstitusu/IYTE/Masters Thesis/R/projects/TCGEx/shiny/genesets/msigdb_hallmark.rda")

ifng_genes<- msigdb_hallmark[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]

ifng_genes <- intersect(ifng_genes, colnames(total))

total[, meta_msigdb_ifng_score := rowMeans(.SD, na.rm = TRUE), .SDcols = ifng_genes]


ayers_ifng_genes <- c("CD3D", "IDO1", "CIITA", "CD3E","CCL5","GZMK","CD2","HLA-DRA","CXCL13","IL2RG","NKG7","HLA-E","CXCR6","LAG3","TAGAP","CXCL10","STAT1","GZMB")

ayers_ifng_genes <- intersect(ayers_ifng_genes, colnames(total))

total[, meta_ayers_ifng_score := rowMeans(.SD, na.rm = TRUE), .SDcols = ayers_ifng_genes]

saveRDS(total, "mel_dfci_2019_TPM_LOG_meta_LOG_all.rds")

nivo_data <- total[meta.ICB == "Nivolumab Therapy", ]

saveRDS(total, "mel_dfci_2019_TPM_LOG_meta_LOG_nivo.rds")

pembro_data <- total[meta.ICB == "Pembrolizumab Therapy", ]

saveRDS(pembro_data, "mel_dfci_2019_TPM_LOG_meta_LOGpembro.rds")
#################################################################################################################################################################

mir_columns <- names(total)[grepl("^MIR", names(total))]

mir_columns <- gsub("-", ".", mir_columns)

wb <- createWorkbook()

addWorksheet(wb, "MirColumns")

writeData(wb, sheet = "MirColumns", x = data.frame(MIR_Columns = mir_columns), startCol = 1, startRow = 1)

saveWorkbook(wb, file = "Mir_Columns_List.xlsx")

#################################################################################################################################################################

#immune response

total <- readRDS("C:/Users/memre/OneDrive - Izmir Yuksek Teknoloji Enstitusu/IYTE/PhD/R/CbioPortal/mel_tsam_liang_2017/mel_tsam_liang_2017_FPKM_all.rds")

response<-fread("data_timeline_treatment.txt")

clin_samp<-fread("data_clinical_sample.txt")

clin_samp<- clin_samp[,1:2]

setnames(clin_samp, "#Patient Identifier", "PATIENT_ID")

merged_response <- merge(clin_samp, response, by = "PATIENT_ID")

rm(clin_samp,response)


setnames(merged_response, gsub(" ", ".", names(merged_response)))

setnames(merged_response, gsub("\\(", ".", names(merged_response)))

setnames(merged_response, gsub(")", ".", names(merged_response)))


setnames(merged_response, colnames(merged_response), paste0("meta.", colnames(merged_response)))


# totally<- merge(total, merged_response, by = "meta.Sample.Identifier")

merged_data <- merge(merged_response, total, by = "meta.Sample.Identifier", all.y = TRUE)


merged_data[, meta.response := 
              ifelse(grepl("^PD", meta.TREATMENT_BEST_RESPONSE), "PD",
              ifelse(grepl("^PR", meta.TREATMENT_BEST_RESPONSE), "PR",
              ifelse(grepl("^CR", meta.TREATMENT_BEST_RESPONSE), "CR",
                     ifelse(grepl("?SD/PD", meta.TREATMENT_BEST_RESPONSE), "PD",
                            ifelse(grepl("^Brief", meta.TREATMENT_BEST_RESPONSE), "SD",
              ifelse(grepl("^SD", meta.TREATMENT_BEST_RESPONSE), "SD", NA_character_)))))), 
           ]

merged_data$meta.response<-as.factor(merged_data$meta.response)

merged_data[, meta_categorized_response := ifelse(merged_data$meta.response %in% c("CR","PR"), 1, 
                                            ifelse(merged_data$meta.response == "PD", 0, NA))]

# merged_data$meta_categorized_response<-NULL

rm(merged_response,total)


merged_data$meta.AGENT <- gsub("Nivo|nivolumab|Nivolumab", "Nivolumab Therapy", merged_data$meta.AGENT)
merged_data$meta.AGENT <- gsub("Pembro|pembrolizumab|Pembrolizumab|MK 3475|MK3475|MK-3475", "Pembrolizumab Therapy", merged_data$meta.AGENT)
merged_data$meta.AGENT <- gsub("Ipilimumab|Ipi+pd1|Vem + ipi", "Ipilimumab Therapy", merged_data$meta.AGENT)


saveRDS(merged_data, "mel_tsam_liang_2017_FPKM_all_updated.rds")


#prepare for tcgex

mel_tsam_liang_2017_FPKM_all_updated$meta.project_id<-mel_tsam_liang_2017_FPKM_all_updated$meta.Cancer.Type.Detailed

mel_tsam_liang_2017_FPKM_all_updated$meta.definition<-mel_tsam_liang_2017_FPKM_all_updated$meta.Sample.Type

mel_tsam_liang_2017_FPKM_all_updated$meta.age_at_diagnosis<-mel_tsam_liang_2017_FPKM_all_updated$meta.Age.at.Diagnosis

mel_tsam_liang_2017_FPKM_all_updated$meta.Age.at.Diagnosis<-NULL

saveRDS(mel_tsam_liang_2017_FPKM_all_updated,"mel_tsam_liang_2017_tcgex.rds")




