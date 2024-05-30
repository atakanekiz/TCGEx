setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(data.table)

rna<-fread("data_mrna_seq_rpkm.txt")

# rna<-fread("data_mrna_seq_rpkm_zscores_ref_all_samples.txt")

clin_pat<-fread("data_clinical_patient.txt")

clin_samp<-fread("data_clinical_sample.txt")

clin <- merge(clin_pat, clin_samp, by = "#Patient Identifier", all = TRUE)

clin<-clin[5:42 , 1:55]



colnames(rna)[1] <- "gene"

suppressPackageStartupMessages(library(org.Hs.eg.db))

hs <- org.Hs.eg.db

conv_names <- AnnotationDbi::select(hs, 
                                    keys=as.character(rna$gene),
                                    columns = c("ENTREZID", "SYMBOL", "external_gene_name"),
                                    keytype = "ENTREZID")


sum(duplicated(conv_names$SYMBOL))

rna$gene <- conv_names$SYMBOL

rna <- rna[!is.na(rna$gene), ]

anyDuplicated(rna$gene)


rna<-as.data.table(t(rna))

colnames(rna) <- as.character(rna[1, ])

rna <- rna[-1, ]

rna_res<-rna

rna<-fread("data_mrna_seq_rpkm.txt")

# rna<-fread("data_mrna_seq_rpkm_zscores_ref_all_samples.txt")

rna$Entrez_Gene_Id<-NULL

rna_res$`#Patient Identifier` <-colnames(rna)


clin_names <- colnames(clin)
clin_names[clin_names != "#Patient Identifier"] <- paste0("meta.", clin_names[clin_names != "#Patient Identifier"])
colnames(clin) <- clin_names


total <- cbind(rna_res, clin[match(rna_res$`#Patient Identifier`, clin$`#Patient Identifier`), -1])

rm(clin,clin_pat,clin_samp,conv_names,rna,rna_res)

colnames(total)[colnames(total) == "#Patient Identifier"] <- "meta.#Patient Identifier"

for (col in colnames(total)) {
  if (!startsWith(col, "meta.")) {
    total[[col]] <- as.numeric(total[[col]])
  } else {
    total[[col]] <- as.factor(total[[col]])
  }
}

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

total$meta.vital_status <- ifelse(total$meta.vital_status == "0:LIVING", 0, 
                                  ifelse(total$meta.vital_status == "1:DECEASED", 1, total$meta.vital_status))


total$meta.days_to_event <- as.numeric(levels(total$meta.days_to_event)[total$meta.days_to_event])

total$meta.days_to_event <- as.numeric(total$meta.days_to_event) * 30


duplicated_columns <- total[, duplicated(names(total))]

for (col in names(total)[duplicated_columns]) {
  avg1 <- mean(total[[col]], na.rm = TRUE)
  avg2 <- mean(total[[paste0(col, ".1")]], na.rm = TRUE)  # dup cols
  if (!is.na(avg1) && !is.na(avg2) && avg1 > avg2) {
    total[[paste0(col, ".1")]] <- NULL  
  } else {
    total[[col]] <- NULL  
  }
}

for (col in colnames(total)) {
  if (!startsWith(col, "meta.")) {
    total[[col]] <- as.numeric(total[[col]])
  }
}

total$`meta.Mutation Load` <- as.numeric(levels(total$`meta.Mutation Load`)[total$`meta.Mutation Load`])
total$`meta.Neo-antigen Load` <- as.numeric(levels(total$`meta.Neo-antigen Load`)[total$`meta.Neo-antigen Load`])

total$meta.patient<-total$`meta.#Patient Identifier`
    
total[, meta.project_id := as.factor(sapply(1:.N, function(x) paste(sample(c(letters, 0:9), 7, replace = TRUE), collapse = "")))]


setnames(total, gsub(" ", ".", names(total)))


load("C:/Users/memre/OneDrive - Izmir Yuksek Teknoloji Enstitusu/IYTE/Masters Thesis/R/projects/TCGEx/shiny/genesets/msigdb_hallmark.rda")

ifng_genes<- msigdb_hallmark[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]

ifng_genes <- intersect(ifng_genes, colnames(total))

total[, meta_msigdb_ifng_score := rowMeans(.SD, na.rm = TRUE), .SDcols = ifng_genes]


ayers_ifng_genes <- c("CD3D", "IDO1", "CIITA", "CD3E","CCL5","GZMK","CD2","HLA-DRA","CXCL13","IL2RG","NKG7","HLA-E","CXCR6","LAG3","TAGAP","CXCL10","STAT1","GZMB")

### ayers_ifng_genes <- intersect(ayers_ifng_genes, colnames(total))

total[, meta_ayers_ifng_score := rowMeans(.SD, na.rm = TRUE), .SDcols = ayers_ifng_genes]

total[, meta_categorized_response := ifelse(total$meta.Durable.Clinical.Benefit %in% c("CR", "PR"), 1, ifelse(total$meta.Durable.Clinical.Benefit == "PD", 0, NA))]



saveRDS(total, "mel_ucla_2016_rpkm_all_pembro.rds")

# saveRDS(total, "mel_ucla_2016_zrpkm_all_pembro.rds")






#################################################################################################################################################################
uploaded_data <- readRDS("C:/Users/memre/OneDrive - Izmir Yuksek Teknoloji Enstitusu/IYTE/PhD/R/CbioPortal/mel_ucla_2016/mel_ucla_2016_rpkm_all_pembro.rds")

df_num = uploaded_data %>% select(where(is.numeric))

df_nonnum = uploaded_data %>% select(-where(is.numeric))

df_nongene = df_num %>% 
  select(starts_with("meta."),starts_with("meta_"))

df_gene <- df_num %>%
  select(-starts_with("meta."), -starts_with("meta_"))

df_gene <- df_gene[, apply(df_gene,2,var)!=0, with=F]

df_gene = log(df_gene +1, base=10)

uploaded_data = cbind(df_gene, df_nonnum,df_nongene)

saveRDS(uploaded_data, "mel_ucla_2016_rpkm_LOG_all_pembro.rds")

#################################################################################################################################################################

##"meta_" log


load("C:/Users/memre/OneDrive - Izmir Yuksek Teknoloji Enstitusu/IYTE/Masters Thesis/R/projects/TCGEx/shiny/genesets/msigdb_hallmark.rda")

ifng_genes<- msigdb_hallmark[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]

ifng_genes <- intersect(ifng_genes, colnames(total))

total[, meta_msigdb_ifng_score := rowMeans(.SD, na.rm = TRUE), .SDcols = ifng_genes]


ayers_ifng_genes <- c("CD3D", "IDO1", "CIITA", "CD3E","CCL5","GZMK","CD2","HLA-DRA","CXCL13","IL2RG","NKG7","HLA-E","CXCR6","LAG3","TAGAP","CXCL10","STAT1","GZMB")

ayers_ifng_genes <- intersect(ayers_ifng_genes, colnames(total))

total[, meta_ayers_ifng_score := rowMeans(.SD, na.rm = TRUE), .SDcols = ayers_ifng_genes]

saveRDS(total, "mel_ucla_2016_rpkm_LOG_meta_LOG_all_pembro.rds")

#################################################################################################################################################################




























#################################################################################################################################################################

cr_data <- total[meta.Durable.Clinical.Benefit == "CR"]
pd_data <- total[meta.Durable.Clinical.Benefit == "PD"]

cr_means <- cr_data[, lapply(.SD, mean), .SDcols = sapply(cr_data, is.numeric)]
pd_means <- pd_data[, lapply(.SD, mean), .SDcols = sapply(pd_data, is.numeric)]

differences <- cr_means - pd_means

sorted_differences <- sort(unlist(differences), decreasing = TRUE)

top_50 <- head(sorted_differences, 50)
bottom_50 <- tail(sorted_differences, 50)

cat("Top 50 Differences:\n")
print(top_50)

cat("\nBottom 50 Differences:\n")
print(bottom_50)

#################################################################################################################################################################

gene_expr <- total[, 1:21835]

# a<- total[, 1:21835]
# 
# a$meta.Durable.Clinical.Benefit<-total$meta.Durable.Clinical.Benefit
# 
# total<-a

# Limma iC'in bir tasarD1m matrisi oluEtur
design <- model.matrix(~0 + meta.Durable.Clinical.Benefit, data=total)

# rownames(design) <- colnames(gene_expr)

total<-t(total)

gene_expr<-t(gene_expr)

# Limma analizi
fit <- lmFit(gene_expr, design)
contrast_matrix <- makeContrasts(CRvsPD = meta.Durable.Clinical.BenefitCR - meta.Durable.Clinical.BenefitPD, levels = design)
fit_contrast <- contrasts.fit(fit, contrast_matrix)
fit_ebayes <- eBayes(fit_contrast)

# Log fold change'e gC6re sD1rala
top_genes <- topTable(fit_ebayes, coef = "CRvsPD", number = Inf, sort.by = "logFC")

# D0lk 10 geni gC6ster
print(top_genes[1:10, ])


#################################################################################################################################################################
# LASSO

library(data.table)
library(glmnet)

data <- readRDS("C:/IYTE/PhD/R/CbioPortal/mel_ucla_2016/mel_ucla_2016.rds")

data <- na.omit(data)

response_column <- "meta.Durable.Clinical.Benefit"

data[[response_column]] <- as.numeric(data[[response_column]])

data[, (response_column) := as.factor(get(response_column))]

predictors <- setdiff(names(data), response_column)

X <- as.matrix(data[, ..predictors])
y <- as.matrix(data[[response_column]])

model <- cv.glmnet(X, y, alpha = 1)

print(model)

plot(model, xvar = "lambda")

coef(model, s = "lambda.1se")

coefficients_df <- coefficients_df %>% 
  arrange(desc(Coefficient))

#################################################################################################################################################################


# total$meta.Durable.Clinical.Benefit <- ifelse(total$meta.Durable.Clinical.Benefit == "0:PR", 0,
#                                               ifelse(total$meta.Durable.Clinical.Benefit == "2:PD", 2,
#                                   ifelse(total$meta.Durable.Clinical.Benefit == "1:CR", 1, total$meta.Durable.Clinical.Benefit)))
# 
# total$meta.Durable.Clinical.Benefit<-as.numeric(total$meta.Durable.Clinical.Benefit)
# 
# a<- total[, 1:21835]
# 
# a$meta.Durable.Clinical.Benefit<-total$meta.Durable.Clinical.Benefit
# 
# total<-a
# 
# design <- model.matrix(~0 + meta.Durable.Clinical.Benefit, data =total)
# fit <- lmFit(gene_expr, design)
# contrast.matrix <- makeContrasts(CR - PD, levels=design)
# fit2 <- contrasts.fit(fit, contrast.matrix)
# fit2 <- eBayes(fit2)
# 
# lfc <- topTable(fit2, coef="CR - PD", number=n, sort.by="logFC")
# 
# top_columns <- lfc$Genes


#################################################################################################################################################################
library(openxlsx)

mir_columns <- names(total)[grepl("^MIR", names(total))]

mir_columns <- gsub("-", ".", mir_columns)

wb <- createWorkbook()

addWorksheet(wb, "MirColumns")

writeData(wb, sheet = "MirColumns", x = data.frame(MIR_Columns = mir_columns), startCol = 1, startRow = 1)

saveWorkbook(wb, file = "Mir_Columns_List.xlsx")


#prepare for tcgex

mel_ucla_2016_rpkm_LOG_meta_LOG_all_pembro$meta.project_id<-mel_ucla_2016_rpkm_LOG_meta_LOG_all_pembro$meta.Cancer.Type.Detailed

mel_ucla_2016_rpkm_LOG_meta_LOG_all_pembro$meta.definition<-mel_ucla_2016_rpkm_LOG_meta_LOG_all_pembro$meta.Biopsy.Time

mel_ucla_2016_rpkm_LOG_meta_LOG_all_pembro$meta.age_at_diagnosis<-mel_ucla_2016_rpkm_LOG_meta_LOG_all_pembro$meta.Age.at.Diagnosis

mel_ucla_2016_rpkm_LOG_meta_LOG_all_pembro$meta.Age.at.Diagnosis<-NULL

saveRDS(mel_ucla_2016_rpkm_LOG_meta_LOG_all_pembro,"mel_ucla_2016_tcgex.rds")
