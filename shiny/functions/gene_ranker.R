#' gene_ranker
#'
#' @description A function to rank genes between two subsets of scRNAseq data. It uses multiple user-selected statistical approaches. This function works with a simple expression data frame created using `SCseqtools::df_extractor()`` function
#'
#' @export
#'
#' @importFrom dplyr pull
#' @importFrom dplyr filter
#' @importFrom dplyr filter_
#' @importFrom dplyr select
#' @importFrom tibble add_column
#' @importFrom magrittr %>%
#' @importFrom BWStest bws_stat
#'
#'
#'
#' @param exprs Expression data frame (Single cells are in rows, genes and metadata are in columns)
#'
#' @param pos_marker A character vector of gene names to positively gate cells (cells expressing these genes will be included in sample and reference)
#'
#' @param neg_marker A character vector of gene names to negatively gate cells (cells expressing these genes will be excluded from sample and reference )
#'
#' @param sample_id String. Name of the samples to use for the ranking. As the genes will be ranked from sample to reference in decreasing order, samples will be on the left side of enrichment plots. Regex based string matching can be applied including logic operations such as `|`. Data for this parameter should match to entries under `sample` column from the `expr` data frame.
#'
#' @param sample_cluster Which clusters to include in the analysis. Data for this parameter should match to entried under `cluster` column from `expr` data frame.
#'
#' @param reference_id Same as `sample_id`, but for reference (right side of the enrichment plots i.e. tail end of the ranking)
#'
#' @param reference_cluster Same as sample_cluster but for reference (right side of the enrichment plots)
#'
#' @param method Specify method to use for gene ranking. It can be one of the following: 's2n' (default), 'ttest', 'difference', 'ratio', 'welch', 'mwt', 'bws'. See PMID: 18344518
#'
#' @param save_ranked_object Logical. Set to TRUE for keeping ranking results as a global variable
#'
#' @param verbose Tells you what is being analyzed.
#'
#'
#' @examples
#' s2n_ranked <- gene_ranker(exprs = exprs,
#'                           pos_marker = c('Cd3e', 'Cd8a'), neg_marker = 'Cd4',
#'                           sample_id = 'young|middle_age', reference_id = 'old_age',
#'                           method = 's2n')
#'
#'
#'


gene_ranker <- function(exprs = NULL, pos_marker = NULL, neg_marker = NULL, sample_id, sample_cluster = NULL, reference_id, reference_cluster = NULL,
    method = "s2n", save_ranked_object = F, verbose = T) {



    # Subset cells that express markers of interest at any level (>0)
    if (!is.null(pos_marker)) {
        exprs <- filter_(exprs, paste(pos_marker, "!= 0", collapse = " & "))
    }

    # Discard cells that express a marker gene which we want to negatively gate in our analyses
    if (!is.null(neg_marker)) {
        exprs <- filter_(exprs, paste(neg_marker, "== 0", collapse = " & "))
    }

    # Subset cells to serve as sample in comparison (sample vs ref comparison i.e. left side of the GSEA plots later on)
    sample_df <- filter(exprs, sample %in% sample_id)



    # Subset clusters of interest
    if (!is.null(sample_cluster)) {

        sample_df <- filter(sample_df, cluster %in% sample_cluster)

    }




    # Subset cells to serve as reference in comparison (sample vs ref comparison i.e. right side of the GSEA plots later on)
    reference_df <- filter(exprs, sample %in% reference_id)




    # Subset clusters of interest
    if (!is.null(reference_cluster)) {

        reference_df <- filter(reference_df, cluster %in% reference_cluster)

    }


    if (method == "mwt") {

        # require(mwt)
        # require(BiocGenerics)

        discard_from_sample <- as.logical(colSums(sample_df == 0, na.rm = T))
        discard_from_ref <- as.logical(colSums(reference_df == 0, na.rm = T))
        discard_from_both <- discard_from_sample | discard_from_ref
        sample_df <- sample_df[, !discard_from_both]
        reference_df <- reference_df[, !discard_from_both]

        sample_df <- add_column(sample_df, Comparison = "sample", .after = 0)

        reference_df <- add_column(reference_df, Comparison = "Reference", .after = 0)




        trimmed_df <- rbind(sample_df, reference_df)
        group_colnames <- as.factor(trimmed_df$Comparison)
        trimmed_df <- trimmed_df %>% dplyr::select(-sample, -cluster, -Comparison) %>% t()

        colnames(trimmed_df) <- group_colnames



        mwt_results <- mwt(trimmed_df, grp = group_colnames)

    } else {


        if (verbose) {

            # Report which cells are being analyzed
            print(paste(dim(sample_df)[1], "cells with the following annotations are included in the 'sample' group (left side of GSEA plot)"))
            print(paste("sample:", paste(levels(droplevels(as.factor(sample_df$sample))), collapse = ", ")))
            print(paste("sample clusters in analysis:", paste(levels(droplevels(as.factor(sample_df$cluster))), collapse = ", ")))


            print(paste(dim(reference_df)[1], "cells with the following annotations are included in the 'reference' group (right side of GSEA plot)"))
            print(paste("Reference:", paste(levels(droplevels(as.factor(reference_df$sample))), collapse = ", ")))
            print(paste("Reference clusters in analysis:", paste(levels(droplevels(as.factor(reference_df$cluster))), collapse = ", ")))

        }

        # Get rid of unnecessary columns
        sample_df <- sample_df[, !colnames(sample_df) %in% c("sample", "cluster")]
        reference_df <- reference_df[, !colnames(reference_df) %in% c("sample", "cluster")]

        sample_mean <- colMeans(sample_df)
        sample_sd <- apply(sample_df, MARGIN = 2, FUN = sd)
        sample_var <- apply(sample_df, MARGIN = 2, FUN = var)
        sample_n <- dim(sample_df)[1]

        reference_mean <- colMeans(reference_df)
        reference_sd <- apply(reference_df, MARGIN = 2, FUN = sd)
        reference_var <- apply(reference_df, MARGIN = 2, FUN = var)
        reference_n <- dim(reference_df)[1]

    }






    if (method == "s2n") {

        ranking <- (sample_mean - reference_mean)/(sample_sd + reference_sd)

    } else if (method == "ttest") {

        ranking <- (sample_mean - reference_mean)/sqrt(((sample_sd^2)/sample_n) + ((reference_sd^2)/reference_n))

    } else if (method == "difference") {

        ranking <- sample_mean - reference_mean

    } else if (method == "ratio") {

        ranking <- sample_mean/reference_mean

    } else if (method == "welch") {

        ranking <- (sample_mean - reference_mean)/sqrt(((sample_var^2)/sample_n) + ((reference_var^2)/reference_n))

    } else if (method == "mwt") {

        ranking <- mwt_results$MWT
        # names(ranking) <- rownames(trimmed_df) #not needed?

    } else if (method == "bws") {

        # require(BWStest)
        # require(dplyr)

        if (colnames(sample_df) != colnames(reference_df)) {

            print("Different column (gene) names between sample and reference. Using shared genes for ranking")

            shared <- intersect(colnames(sample_df), colnames(reference_df))

            sample_df <- sample_df[, shared]
            reference_df <- reference_df[, shared]

            # ensure match
            if (colnames(sample_df) != colnames(reference_df))
                error("Column (gene) names still do not match between sample and reference")

        }

        ranking <- numeric()

        for (i in colnames(sample_df)) {

            sample_gene <- pull(sample_df, i)
            reference_gene <- pull(reference_df, i)

            ranking[i] <- bws_stat(sample_gene, reference_gene)

        }

    } else {
        stop("Select one of the following statistical tests: s2n, ttest, difference, ratio, welch, mwt, bws")
    }


    ranking <- subset(ranking, c(!is.na(ranking) & !is.infinite(ranking)))
    ranking <- sort(ranking, decreasing = T)

    if (save_ranked_object == T)
        {
            assign("ranked_geneset", ranking, envir = .GlobalEnv)
        }  # Useful for iterative analysis when executed in pipeline

    ranking  # Return ranked gene set

}

