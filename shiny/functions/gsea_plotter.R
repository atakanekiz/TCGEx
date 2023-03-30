#' gsea_plotter
#'
#' @description A wrapper function to create gsea plots.
#'
#' @export
#'
#' @importFrom gridExtra grid.arrange
#' @importFrom stringr str_replace_all
#' @importFrom dplyr pull
#' @importFrom dplyr filter
#' @importFrom dplyr filter_
#' @importFrom dplyr select
#'
#' @import ggplot2
#' @import fgsea
#'
#' @param exprs A data frame containing the gene expression values, samples and cluster information (in columns) from single cells (rows). It must include lower case column name for `sample` and `cluster`. This can be generated from a Seurat object using `df_extractor` function.
#'
#' @param preranked_genes A named numeric vector for ranked gene expression values. This can be obtained using `gene_ranker` function. When this is provided, the arguments for subsetting and ranking `exprs` data frame is ignored. The idea behind this argument is to speed up GSEA analyses on the same ranked gene expression data using different gene sets without the need to recalculate the ranking in each iteration.
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
#' @param gene_set reference pathways. It should be a named list object with character vector of gene names. A custom list can be provided or one of the pre-loaded lists can be selected: msigdb_all, msigdb_curated, msigdb_go, msigdb_hallmark (default), msigdb_immune, msigdb_motif.
#'
#' @param minSize Minimum number of overlapping genes to consider a pathway in analysis
#'
#' @param maxSize Maximum number of overallping genes to consider a pathway in analysis
#'
#' @param top_n Specificy how many of the top enriched pathways are shown
#'
#' @param gseaParam fgsea parameter to change bar sizes
#'
#' @param plot_individual Select pathway to plot
#'
#' @param append_title Add informative titles to individual plots
#'
#' @param top_plots_title Add informative titles to summary plots
#'
#' @param seed Random seed
#'
#' @param keep_results Logical to store results globally
#'
#' @param save_png Save resulting plot as png
#'
#' @param png_units Png size units
#'
#' @param png_width width
#'
#' @param png_height height
#'
#' @param append_to_filename Custom string to add to png file
#'
#' @param verbose Prints the analysis details
#'
#' @param annot_text_color color of NES-p val annotation text
#'
#' @param annot_text_size size of NES-pval annotation text
#'
#' @param annot_text_fontface fontface of annotation text (1,2,3,4, plain-bold-italic-bold and italic)
#'
#' @param ... Additional parameters to pass into fgsea function call. This can include things like nperm (omit for the recommended `fgseaMultilevel` call, provide a value for `fgseaSimple` function call), minSize, maxSize etc.
#'
#'
#'


gsea_plotter <- function(exprs = NULL, preranked_genes = NULL, pos_marker = NULL, neg_marker = NULL, sample_id = NULL, sample_cluster = NULL,
    reference_id = NULL, reference_cluster = NULL, method = "s2n", gene_set = msigdb_hallmark, minSize = 1, maxSize = ncol(exprs)-3, top_n = 10,
    gseaParam = 1, plot_individual = NULL, append_title = F, top_plots_title = T, seed = 123, keep_results = T, save_png = F, png_units = "in",
    png_width = 4, png_height = 3, append_to_filename = "", verbose = T, annot_text_color = "black", annot_text_size = 4, annot_text_fontface = 2, ...) {

    set.seed(seed)
  
  if (gene_set == "hallmark") {
    
    load(file='genesets/msigdb_hallmark.rda')
    
    gene_set <-msigdb_hallmark
    
  } else if (gene_set == "go") {
    
    load(file='genesets/msigdb_go.rda')
    
    gene_set <- msigdb_go
    
  } else if (gene_set == "curated") {
    
    load(file='genesets/msigdb_curated.rda')
    
    gene_set <- msigdb_curated
    
  } else if (gene_set == "immune") {
    
    load(file='genesets/msigdb_immune.rda')
    
    gene_set <- msigdb_immune
    
  } else if (gene_set == "motif") {
    
    load(file='genesets/msigdb_motif.rda')
    
    gene_set <- msigdb_motif
    
  } else if (gene_set == "all") {
    
    load(file='genesets/msigdb_all.rda')
    
    gene_set <- msigdb_all
    
  }
  
    # Read molecular signatures database (MSigDB) gene lists. files must be stored
    # if (gene_set == "hallmark") {
    #
    #     gene_set <- msigdb_hallmark
    #
    # }

    if (is.null(reference_cluster))
        {
            reference_cluster <- sample_cluster
        }  # Specify only sample_cluster for subsetting on the same clusters in sample and reference

    if (!is.null(preranked_genes)) {

        ranked_genes = preranked_genes

    } else {

        ranked_genes <- gene_ranker(exprs = exprs, pos_marker = pos_marker, neg_marker = neg_marker, sample_id = sample_id, sample_cluster = sample_cluster,
            reference_id = reference_id, reference_cluster = reference_cluster, method = method, verbose = verbose)
    }

    if (keep_results) assign("ranked_genes", ranked_genes, .GlobalEnv)




  # browser()

    res <- fgsea(pathways = gene_set, stats = ranked_genes, minSize = minSize, maxSize = maxSize, ...)

    if (keep_results)
        assign("gsea_res", res, .GlobalEnv)

    if (is.null(plot_individual)) {

        if (top_plots_title == T) {

            arg_list <- list(samp_clu = sample_cluster, ref_clu = reference_cluster, pos = paste(pos_marker, collapse = "."), neg = paste(neg_marker,
                collapse = "."))
            select_non_null <- !sapply(arg_list, function(x) {
                identical(x, "")
            })



            main_title <- paste(sample_id, "vs", reference_id)
            plot_subtitle <- paste(names(arg_list[select_non_null]), arg_list[select_non_null], sep = ": ", collapse = "__")

            plot_title <- paste0(main_title, "\n", plot_subtitle)

        } else {
            plot_title = ""
        }

        plot_grob <- top_plotter(gsea_results = res, ranked_genes = ranked_genes, gene_set = gene_set, top_n = top_n, gseaParam = gseaParam,
            plot_title = plot_title, do.plot = F)
        grid.arrange(plot_grob)


    } else {


        hits <- c(grep(plot_individual, res$pathway, ignore.case = T, value = T))

        # assign('hits', hits, .GlobalEnv)


        if (length(hits) > 1) {
            multiple_hits <- t(t(hits))
            colnames(multiple_hits) <- "Multiple pathway matches"
            rownames(multiple_hits) <- c(1:length(hits))
            print(multiple_hits)
            num <- as.numeric(readline(prompt = "Multiple pathways are found. Select a number from the list above "))


            while (!num %in% 1:length(hits)) {
                num <- as.numeric(readline(prompt = paste0("Please pick a number between 1 and ", length(hits), ":    ")))
            }
            assign("num", num, .GlobalEnv)

            annot_padj <- signif(as.numeric(res[res$pathway == hits[num], "padj"]), digits = 2)
            annot_NES <- signif(as.numeric(res[res$pathway == hits[num], "NES"]), digits = 2)
            annot_ES <- signif(as.numeric(res[res$pathway == hits[num], "ES"]), digits = 2)

            # grob<- grobTree(textGrob(paste('adj.p: ', annot_padj, '\nNES:', annot_NES), x= 0.1, y=annot_ES, hjust = 0, gp = gpar(col='red',
            # fontsize=3, fontface='italic')))

            annot_text <- paste("adj.p: ", annot_padj, "\nNES:", annot_NES)

            if (append_title == F) {

                # plotEnrichment(pathway = gene_set[[hits[num]]], stats = ranked_genes) + labs(title = hits[num]) + annotation_custom(grob)

                plot_grob <- plotEnrichment(pathway = gene_set[[hits[num]]], stats = ranked_genes) + labs(title = hits[num]) + annotate("text",
                  x = x_pos, y = annot_ES/2, label = annot_text, colour = annot_text_color, size = annot_text_size, fontface = annot_text_fontface) +
                  theme(plot.title = element_text(size = 5, hjust = 0.5))
                print(plot_grob)

            } else {

                arg_list <- list(samp_clu = sample_cluster, ref_clu = reference_cluster, pos = pos_marker, neg = neg_marker)
                select_non_null <- !sapply(arg_list, function(x) {
                  identical(x, "")
                })
                select_non_null2 <- !sapply(arg_list, is.null)
                select_non_null <- as.logical(select_non_null * select_non_null2)

                plot_subtitle <- paste(names(arg_list[select_non_null]), arg_list[select_non_null], sep = ": ", collapse = "__")



                # plotEnrichment(pathway = gene_set[[hits[num]]], stats = ranked_genes) + labs(title = paste0(hits[num], sample_id, ' vs ', reference_id),
                # subtitle = plot_subtitle)+ annotation_custom(grob)

                plot_grob <- plotEnrichment(pathway = gene_set[[hits[num]]], stats = ranked_genes) + labs(title = hits[num], subtitle = paste(sample_id,
                  "vs", reference_id, plot_subtitle)) + annotate("text", x = x_pos, y = annot_ES/2, label = annot_text, colour = annot_text_color,
                  size = annot_text_size, fontface = annot_text_fontface) + theme(plot.title = element_text(size = 10, hjust = 0.5), plot.subtitle = element_text(size = 6,
                  hjust = 0.5))
                print(plot_grob)

            }



            ##########################

        } else {

            if (verbose)
                message(paste("Plotting", hits))

            annot_padj <- signif(as.numeric(res[res$pathway == hits, "padj"]), digits = 2)
            annot_NES <- signif(as.numeric(res[res$pathway == hits, "NES"]), digits = 2)
            annot_ES <- signif(as.numeric(res[res$pathway == hits, "ES"]), digits = 2)
            x_pos <- length(ranked_genes)/4

            # grob<- grobTree(textGrob(paste('adj.p: ', annot_padj, '\nNES:', annot_NES), x= 0.1, y=annot_ES, hjust = 0, gp = gpar(col='red',
            # fontsize=3, fontface='italic')))

            annot_text <- paste("adj.p: ", annot_padj, "\nNES:", annot_NES)


            if (append_title == F) {

                # plotEnrichment(pathway = gene_set[[hits]], stats = ranked_genes) + labs(title = hits) + annotation_custom(grob)

                plot_grob <- plotEnrichment(pathway = gene_set[[hits]], stats = ranked_genes) + labs(title = hits) + annotate("text", x = x_pos,
                  y = annot_ES/2, label = annot_text, colour = annot_text_color, size = annot_text_size, fontface = annot_text_fontface)
                print(plot_grob)

            } else {

                arg_list <- list(samp_clu = sample_cluster, ref_clu = reference_cluster, pos = paste(pos_marker, collapse = "."), neg = paste(neg_marker,
                  collapse = "."))

                select_non_null <- !sapply(arg_list, function(x) {
                  identical(x, "")
                })
                select_non_null2 <- !sapply(arg_list, is.null)
                select_non_null <- as.logical(select_non_null * select_non_null2)

                plot_subtitle <- paste(names(arg_list[select_non_null]), arg_list[select_non_null], sep = ": ", collapse = " ")

                # plotEnrichment(pathway = gene_set[[hits]], stats = ranked_genes) + labs(title = paste0(hits, sample_id, ' vs ', reference_id), subtitle =
                # plot_subtitle)+ annotation_custom(grob)

                plot_grob <- plotEnrichment(pathway = gene_set[[hits]], stats = ranked_genes) + labs(title = hits, subtitle = paste(sample_id,
                  "vs", reference_id, plot_subtitle)) + annotate("text", x = x_pos, y = annot_ES/2, label = annot_text, colour = annot_text_color,
                  size = annot_text_size, fontface = annot_text_fontface) + theme(plot.title = element_text(size = 10, hjust = 0.5), plot.subtitle = element_text(size = 6,
                  hjust = 0.5))
                print(plot_grob)
            }

        }

    }

    if (save_png == T)
        {

            sample_cluster <- str_replace_all(sample_cluster, "[:punct:]|[:space:]", "")
            reference_cluster <- str_replace_all(reference_cluster, "[:punct:]|[:space:]", "")
            sample_id <- str_replace_all(sample_id, "[:punct:]|[:space:]", "")
            reference_id <- str_replace_all(reference_id, "[:punct:]|[:space:]", "")

            arg_list <- list(SAMPclus = sample_cluster, REFclus = reference_cluster, SAMPid = sample_id, REFid = reference_id, POSmarker = paste(pos_marker,
                collapse = "."), NEGmarker = paste(neg_marker, collapse = "."))

            select_non_null <- !sapply(arg_list, function(x) {
                identical(x, "")
            })
            select_non_null2 <- !sapply(arg_list, is.null)
            select_non_null <- as.logical(select_non_null * select_non_null2)




            if (sum(select_non_null) == 0) {

                filename <- paste0("Unsubsetted_data_plots__", ".pdf")

            } else {

                filename <- paste(names(arg_list[select_non_null]), arg_list[select_non_null], sep = "-", collapse = "___")
                filename <- paste0(filename, "___", append_to_filename, ".png")
            }

            ggsave(plot = plot_grob, filename = filename, width = png_width, height = png_height, units = png_units)

        }  #else grid.arrange(plot_grob)

}

############################################################################################################## # Examples of using master_gsea function

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) exprs <- readRDS('aging_exprs.rds')

# exprs <- readRDS('aging_exprs.rds')

# gsea_plotter(exprs, sample_id = 'Young \\(WT\\)', reference_id = 'Aged \\(WT\\)', sample_cluster = 'NK', reference_cluster = 'NK',
# gene_set = 'hallmark')




