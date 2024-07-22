setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(data.table)

rna<-fread("GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv")

clin_corrected<-fread("SampleTableCorrected.9.19.16.csv")

clin<-fread("bms038_clinical_data.csv")

cytolytic_score<-fread("GSE91061_BMS038109Sample_Cytolytic_Score_20161026.txt")

#as.factor

columns_to_factor <- setdiff(names(clin), c("PFS_SOR", "OS_SOR", "OS", "OSWK", "PFS", "PFSWK"))

clin[, (columns_to_factor) := lapply(.SD, as.factor), .SDcols = columns_to_factor]

columns_to_factor <- setdiff(names(clin_corrected), "IDX")

clin_corrected[, (columns_to_factor) := lapply(.SD, as.factor), .SDcols = columns_to_factor]

rm(columns_to_factor)

#merge

clin<-clin[, -(2:8)]

clin_corrected$PatientID.replicate<-NULL

setnames(clin_corrected, old = names(clin_corrected)[11], new = "benefit")

clin_corrected$Response<-as.factor(clin_corrected$Response)

clinical <- merge(clin_corrected, clin, by = "PatientID", all.x = TRUE)

rm(clin,clin_corrected)

#cytolytic

clinical$BamName <- sub("\\.bam$", "", clinical$BamName)

colnames(cytolytic_score)[colnames(cytolytic_score) == "V1"] <- "BamName"

#merge

clinical <- merge(clinical, cytolytic_score, by = "BamName", all.x = TRUE)

rm(cytolytic_score)


#gene annotation

colnames(rna)[1] <- "gene"

col_rna<- colnames(rna[,2:110])

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

rm(conv_names,hs)

rna<-as.data.table(t(rna))

colnames(rna) <- as.character(rna[1, ])

rna <- rna[-1, ]

rna$BamName <-col_rna

rm(col_rna)

#meta.cols

colnames(rna)[colnames(rna) == "BamName"] <- "meta.BamName"

colnames(clinical) <- paste0("meta.", colnames(clinical))

#merge

riaz <- merge(rna, clinical, by = "meta.BamName", all.x = TRUE)

rm(clinical,new_colnames,rna)

colnames(riaz)[colnames(riaz) == "meta.V2"] <- "meta.cytolytic_score"

#col arrangements

colnames(riaz)[colnames(riaz) == "meta.OS"] <- "meta.days_to_event"

colnames(riaz)[colnames(riaz) == "meta.OS_SOR"] <- "meta.vital_status"

colnames(riaz)[colnames(riaz) == "meta.PreOn"] <- "meta.definition"


#as numeric

for (col in colnames(riaz)) {
  if (!startsWith(col, "meta.")) {
    riaz[[col]] <- as.numeric(riaz[[col]])
  }
}

#meta.ifng


load("C:/Users/memre/OneDrive - Izmir Yuksek Teknoloji Enstitusu/IYTE/Masters Thesis/R/projects/TCGEx/shiny/genesets/msigdb_hallmark.rda")

ifng_genes<- msigdb_hallmark[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]

ifng_genes <- intersect(ifng_genes, colnames(riaz))

riaz[, meta_msigdb_ifng_score := rowMeans(.SD, na.rm = TRUE), .SDcols = ifng_genes]


ayers_ifng_genes <- c("CD3D", "IDO1", "CIITA", "CD3E","CCL5","GZMK","CD2","HLA-DRA","CXCL13","IL2RG","NKG7","HLA-E","CXCR6","LAG3","TAGAP","CXCL10","STAT1","GZMB")

ayers_ifng_genes <- intersect(ayers_ifng_genes, colnames(riaz))

riaz[, meta_ayers_ifng_score := rowMeans(.SD, na.rm = TRUE), .SDcols = ayers_ifng_genes]

riaz[, meta_categorized_response := ifelse(riaz$meta.BOR %in% c("CR", "PR"), 1, 
                                            ifelse(riaz$meta.BOR == "PD", 0, NA))]



riaz[, meta.project_id := as.factor(sapply(1:.N, function(x) paste(sample(c(letters, 0:9), 7, replace = TRUE), collapse = "")))]

riaz[, meta.gender := as.factor(sapply(1:.N, function(x) paste(sample(c(letters, 0:9), 7, replace = TRUE), collapse = "")))]

riaz$meta.patient<-riaz$meta.PatientID



saveRDS(riaz, "riaz.rds")




#################################################################################################################################################################

library(readxl)

mir_columns <- names(riaz)[grepl("^MIR", names(riaz))]

mir_columns <- gsub("-", ".", mir_columns)

wb <- createWorkbook()

addWorksheet(wb, "MirColumns")

writeData(wb, sheet = "MirColumns", x = data.frame(MIR_Columns = mir_columns), startCol = 1, startRow = 1)

saveWorkbook(wb, file = "Mir_Columns_List.xlsx")

#prepare for tcgex

riaz $meta.project_id<-riaz $meta.PopCateg

# riaz $meta.definition<-riaz $meta.Biopsy.Context

# riaz $meta.age_at_diagnosis<-riaz $meta.Age

# riaz $meta.Age<-NULL

saveRDS(riaz ,"riaz_2017_tcgex.rds")
