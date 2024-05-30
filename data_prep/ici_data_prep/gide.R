setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(data.table)
library(readxl)
library(openxlsx)
library(dplyr)

rna<-fread("cancercell_normalized_counts_genenames.txt")

#preparing pd1 data

pd1_data <- rna%>%select(starts_with("PD1")) 

cols_pd1<-colnames(pd1_data)

pd1_data<-as.data.table(t(pd1_data))

colnames(pd1_data)<-rna$Gene

pd1_data$meta.patient<-cols_pd1

pd1_data$meta.project_id <- gsub(".*PD1_(\\d+)(?:_EDT|_PRE).*", "\\1", pd1_data$meta.patient)

#preparing dual ICB data

dual_data <- rna%>%select(starts_with("ipi")) 

cols_dual<-colnames(dual_data)

dual_data<-as.data.table(t(dual_data))

colnames(dual_data)<-rna$Gene

dual_data$meta.patient<-cols_dual

rm(cols_dual,cols_pd1)

dual_data$meta.project_id <- gsub(".*PD1_(\\d+)(?:_EDT|_PRE).*", "\\1", dual_data$meta.patient)


#preparing dual ICB clinical data

clin_dual<-read_excel("1-s2.0-S1535610819300376-mmc3.xlsx")

setnames(clin_dual, gsub(" ", ".", names(clin_dual)))

setnames(clin_dual, gsub("\\(", ".", names(clin_dual)))

setnames(clin_dual, gsub(")", ".", names(clin_dual)))

colnames(clin_dual) <- paste0("meta.", colnames(clin_dual))

clin_dual$meta.project_id<-clin_dual$meta.Patient.no.

clin_dual$meta.gender<-clin_dual$meta.Sex

clin_dual <- as.data.frame(lapply(clin_dual, function(x) {
  if(is.character(x)) {
    as.factor(x)
  } else {
    x
  }
}))

colnames(clin_dual)[colnames(clin_dual) == "meta.Overall.Survival..Days."] <- "meta.days_to_event"

colnames(clin_dual)[colnames(clin_dual) == "meta.Last.Followup.Status"] <- "meta.vital_status"

clin_dual$meta.vital_status <- ifelse(clin_dual$meta.vital_status == "Alive", 0, 
                                  ifelse(clin_dual$meta.vital_status == "Dead, melanoma", 1, clin_dual$meta.vital_status))

#preparing dual PD1 clinical data

clin_pd1<-read_excel("1-s2.0-S1535610819300376-mmc2.xlsx")

setnames(clin_pd1, gsub(" ", ".", names(clin_pd1)))

setnames(clin_pd1, gsub("\\(", ".", names(clin_pd1)))

setnames(clin_pd1, gsub(")", ".", names(clin_pd1)))

colnames(clin_pd1) <- paste0("meta.", colnames(clin_pd1))

clin_pd1$meta.project_id<-clin_pd1$meta.Patient.no.

clin_pd1$meta.gender<-clin_pd1$meta.Sex

clin_pd1 <- as.data.frame(lapply(clin_pd1, function(x) {
  if(is.character(x)) {
    as.factor(x)
  } else {
    x
  }
}))

colnames(clin_pd1)[colnames(clin_pd1) == "meta.Overall.Survival..Days."] <- "meta.days_to_event"

colnames(clin_pd1)[colnames(clin_pd1) == "meta.Last.Followup.Status"] <- "meta.vital_status"

clin_pd1$meta.vital_status <- ifelse(clin_pd1$meta.vital_status == "Alive", 0, 
                                  ifelse(clin_pd1$meta.vital_status == "Dead, melanoma", 1,
                                  ifelse(clin_pd1$meta.vital_status == "Dead",1,clin_pd1$meta.vital_status)))

#merge pd1

col_index <- 5491
col_name <- names(pd1_data)[col_index] 
new_col_name <- paste0(col_name, ".1")  
names(pd1_data)[col_index] <- new_col_name 

col_index <- 18791
col_name <- names(pd1_data)[col_index] 
new_col_name <- paste0(col_name, ".1")  
names(pd1_data)[col_index] <- new_col_name  

clin_pd1$meta.project_id<-as.character(clin_pd1$meta.project_id)

pd1_gide <- merge(pd1_data, clin_pd1, by = "meta.project_id", all.x = TRUE)

pd1_gide$meta.Biopsy_Time <- substr(pd1_gide$meta.patient, nchar(pd1_gide$meta.patient) - 2, nchar(pd1_gide$meta.patient))

pd1_gide$meta.Biopsy_Time <- as.factor(pd1_gide$meta.Biopsy_Time)


#merge dual

clin_dual$meta.project_id<-as.character(clin_dual$meta.project_id)

anyDuplicated(colnames(dual_data))

col_index <- 5491
col_name <- names(dual_data)[col_index] 
new_col_name <- paste0(col_name, ".1")  
names(dual_data)[col_index] <- new_col_name 

col_index <- 18791
col_name <- names(dual_data)[col_index] 
new_col_name <- paste0(col_name, ".1")  
names(dual_data)[col_index] <- new_col_name 

dual_gide <- merge(dual_data, clin_pd1, by = "meta.project_id", all.x = TRUE)

dual_gide$meta.Biopsy_Time <- substr(dual_gide$meta.patient, nchar(dual_gide$meta.patient) - 2, nchar(dual_gide$meta.patient))

dual_gide$meta.Biopsy_Time <- as.factor(dual_gide$meta.Biopsy_Time)

#meta.ifng pd1

load("C:/Users/memre/OneDrive - Izmir Yuksek Teknoloji Enstitusu/IYTE/Masters Thesis/R/projects/TCGEx/shiny/genesets/msigdb_hallmark.rda")

ifng_genes<- msigdb_hallmark[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]

ifng_genes <- intersect(ifng_genes, colnames(pd1_gide))

pd1_gide[, meta_msigdb_ifng_score := rowMeans(.SD, na.rm = TRUE), .SDcols = ifng_genes]


ayers_ifng_genes <- c("CD3D", "IDO1", "CIITA", "CD3E","CCL5","GZMK","CD2","HLA-DRA","CXCL13","IL2RG","NKG7","HLA-E","CXCR6","LAG3","TAGAP","CXCL10","STAT1","GZMB")

ayers_ifng_genes <- intersect(ayers_ifng_genes, colnames(pd1_gide))

pd1_gide[, meta_ayers_ifng_score := rowMeans(.SD, na.rm = TRUE), .SDcols = ayers_ifng_genes]

pd1_gide[, meta_categorized_response := ifelse(pd1_gide$meta.Best.RECIST.response %in% c("CR", "PR"), 1, 
                                           ifelse(pd1_gide$meta.Best.RECIST.response == "PD", 0, NA))]

pd1_gide$meta.definition<- pd1_gide$meta.Biopsy_Time

pd1_gide[, meta.ICB.Response := ifelse(pd1_gide$meta_categorized_response ==1, "responder", ifelse(pd1_gide$meta_categorized_response == 0, "non-responder", NA))]


#meta.ifng dual

load("C:/Users/memre/OneDrive - Izmir Yuksek Teknoloji Enstitusu/IYTE/Masters Thesis/R/projects/TCGEx/shiny/genesets/msigdb_hallmark.rda")

ifng_genes<- msigdb_hallmark[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]

ifng_genes <- intersect(ifng_genes, colnames(dual_gide))

dual_gide[, meta_msigdb_ifng_score := rowMeans(.SD, na.rm = TRUE), .SDcols = ifng_genes]


ayers_ifng_genes <- c("CD3D", "IDO1", "CIITA", "CD3E","CCL5","GZMK","CD2","HLA-DRA","CXCL13","IL2RG","NKG7","HLA-E","CXCR6","LAG3","TAGAP","CXCL10","STAT1","GZMB")

ayers_ifng_genes <- intersect(ayers_ifng_genes, colnames(dual_gide))

dual_gide[, meta_ayers_ifng_score := rowMeans(.SD, na.rm = TRUE), .SDcols = ayers_ifng_genes]

dual_gide[, meta_categorized_response := ifelse(dual_gide$meta.Best.RECIST.response %in% c("CR", "PR"), 1, 
                                           ifelse(dual_gide$meta.Best.RECIST.response == "PD", 0, NA))]

dual_gide$meta.definition<- dual_gide$meta.Biopsy_Time

dual_gide[, meta.ICB.Response := ifelse(dual_gide$meta_categorized_response ==1, "responder", ifelse(dual_gide$meta_categorized_response == 0, "non-responder", NA))]

#log normalize

log_normalize <- function(x) {
  if (is.numeric(x)) {
    return(log(x+1, base=10))
  } else {
    return(x) 
  }
}

for (col_name in names(pd1_gide)) {
  if (!startsWith(col_name, "meta.")) {
    pd1_gide[[col_name]] <- log_normalize(pd1_gide[[col_name]])
  }
}

#Save the data

saveRDS(pd1_gide, "pd1_gide.rds")

saveRDS(dual_gide, "dual_gide.rds")

#lasso mirnas

library(readxl)

mir_columns <- names(dual_gide)[grepl("^MIR", names(dual_gide))]

mir_columns <- gsub("-", ".", mir_columns)

wb <- createWorkbook()

addWorksheet(wb, "MirColumns")

writeData(wb, sheet = "MirColumns", x = data.frame(MIR_Columns = mir_columns), startCol = 1, startRow = 1)

saveWorkbook(wb, file = "Mir_Columns_List.xlsx")


#prepare for tcgex

dual_gide[, meta.project_id := "melanoma"]

dual_gide$meta.definition<-dual_gide$meta.Biopsy_Time

dual_gide$meta.age_at_diagnosis<-dual_gide$meta.Age..Years.

dual_gide$meta.Age..Years.<-NULL

saveRDS(dual_gide ,"dual_gide_tcgex.rds")
