BiocManager::install("AnnotationHub")

library(AnnotationHub)
hub <- AnnotationHub()
query(hub, c("homo sapiens","ensdb"))


ensdb <- hub[["AH113665"]]

require("ensembldb")
gns <- genes(ensdb)
gns


lncs <- gns[gns$gene_biotype %in% "lncRNA"]
lncs

lncRNAs<- lncs@elementMetadata@listData[["gene_name"]]

lncRNAs1<-as.vector(lncRNAs)

filtered_lncrnas <- lncRNAs1[lncRNAs1 != ""]

filtered_lncrnas <- gsub("-", ".", filtered_lncrnas)

saveRDS(filtered_lncrnas, file = "filtered_lncrnas.rds")

