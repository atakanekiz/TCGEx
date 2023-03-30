#' df_extractor
#'
#' @export
#'
#' @import Seurat
#'
#' @importFrom tibble add_column
#' @importFrom dplyr pull
#' @importFrom dplyr filter
#' @importFrom dplyr filter_
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#'
#' @description A function to obtain a simple data frame containing gene
#' expression, sample and cluster information from Seurat object for
#' programmatic plotting. This function is written for Seurat v3 objects.
#'
#' @param seurat_obj A Seurat object to extract from
#'
#' @param metadata_to_extract A named character vector to pull metadata
#' information from the Seurat object. The names provided here will constitute
#' column names in the extracted data frame. The items of the character vector
#' should match the column names of the `seurat_obj@metadata` data frame.
#'
#' @param assay Specify assay data that should be extracted from `seurat_obj`.
#' If `NULL`, the default assay will be extracted
#'
#' @param slot Specify data slot within the `Assay` that should be used.
#' If `NULL`, `data` (normalized data) slot will be extracted
#'
#' @param humanize Logical. Set to `TRUE` to convert gene names to human names.
#'  This is useful to generate GSEA-ready data frames to be used with MSigDB
#'  gene lists. See below for method for humanization.
#'
#' @param humanize_method Select how to perform gene conversion. Could be one of the following "uppercase", "biomart". When "uppercase" is selected, gene names are converted to uppercase letters which would be suitable for most of the gene conversions from mouse to human. This is a quick conversion that does not require biomart matching. When this parameter is set to "biomart", Ensembl-Biomart database is searched for annotated orthologs. This parameters does not matter when `humanize = FALSE`.
#'
#' @return A data frame containing genes in columns, and cells/observations in the rows. This data frame can be directly used for ggplot-friendly plotting
#'
#' @examples
#'
#'   exprs_df <- df_extractor(seurat_obj,
#' metadata_to_extract = c(Cluster="new.idents", Sample = "group"),
#' humanize = F)
#'
#'   exprs_hum <- df_extractor(seurat_obj,
#' metadata_to_extract = c(Cluster="new.idents", Sample = "group"),
#' humanize = T, humanize_method = "uppercase")
#'



df_extractor <- function(seurat_obj,
                         metadata_to_extract = c(cluster="orig_clusters", sample="group"),
                         assay = NULL,
                         slot = NULL,
                         humanize = F,
                         humanize_method = "uppercase"


){


    if(is.null(assay)) assay = DefaultAssay(seurat_obj)
    if(is.null(slot)) slot = "data"

    exprs <- as.data.frame(as.matrix(GetAssayData(object = seurat_obj,
                                                  assay = assay,
                                                  slot = slot)))






    if(humanize==T){

        if(humanize_method == "biomart"){

        exprs <- add_column(exprs, gene = rownames(exprs), .after = 0)

        require(biomaRt)
        require(data.table)


        mouse_genes <- exprs$gene

        human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
        mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="www.ensembl.org")

        convert_genes = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = mouse_genes , mart = mouse,
                               attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

        exprs <- merge(convert_genes, exprs, by.x = "MGI.symbol", by.y="gene", sort = F)

        exprs <- exprs[,colnames(exprs) != "MGI.symbol"]

        exprs <- setDT(exprs)[,lapply(.SD, mean), by=HGNC.symbol] # MUCH FASTER than aggregate and dplyr functions

        human_genes <- exprs %>% pull(HGNC.symbol)  # Capture gene names to give it back as row names

        exprs <- as.data.frame(exprs)

        exprs <- exprs[, colnames(exprs) != "HGNC.symbol"]

        rownames(exprs) <- human_genes


    } else if (humanize_method == "uppercase"){

        newnames <- toupper(rownames(exprs))

        if(anyDuplicated(newnames) == 0){

            rownames(exprs) <- newnames

        } else {

            dups <- duplicated(newnames)

            newnames <- newnames[!dups]

            exprs <- exprs[!dups,]

            rownames(exprs) <- newnames
        }




    } else stop("If humanize=T, humanize_method can be one of 'uppercase' or 'biomart'")

    }

        # Transpose exprs
    exprs <- as.data.frame(t(as.matrix(exprs)))

    if(is.null(names(metadata_to_extract))) {

        warning("'metadata_to_extract' argument is an unnamed vector. Column names are created automatically. To ensure correct functioning of other downstream functions, ensure the column names are appropriate.")
        names(metadata_to_extract) <- metadata_to_extract

    }

    for(i in names(metadata_to_extract)){

        exprs <- add_column(exprs, !!i := seurat_obj[[metadata_to_extract[[i]]]][[1]], .after = 0)
    }

    exprs

}

