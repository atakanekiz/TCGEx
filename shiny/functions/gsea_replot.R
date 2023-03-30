#' gsea_replot
#'
#' @description # This is a function to speed up GSEA plotting in batch execution. A global gsea_res object will be created by using gsea_plotter function. This script is an excerpt from gsea_plotter function to recycle results file to create graphs quickly. Only the graphing engine is used from gsea_plotter function preventing recalculation of the results.
#'
#' @export
#'
#' @importFrom dplyr pull
#' @importFrom dplyr filter
#' @importFrom dplyr filter_
#' @importFrom dplyr select
#' @importFrom gridExtra grid.arrange
#'
#' @import fgsea
#' @import ggplot2

#'
#' @param res Previously calculated gsea results (global object created with `gsea_plotter()` works)
#'
#' @param ranked_genes Named numeric vector of ranked genes (automatically created after gsea_plotter function, or can be provided separately)
#'
#' @param pos_marker A character vector of gene names to positively gate cells (cells expressing these genes will be included in sample and reference)
#'
#' @param neg_marker A character vector of gene names to negatively gate cells (cells expressing these genes will be excluded from sample and reference )
#'
#' @param sample_id String. Name of the samples to use for the ranking. As the genes will be ranked from sample to reference in decreasing order, samples will be on the left side of enrichment plots. Regex based string matching can be applied including logic operations such as `|`. Data for this parameter should match to entries under `sample` column from the `expr` data frame. Make sure you escape special characters such as parantheses and periods.
#'
#' @param sample_cluster Which clusters to include in the analysis. Data for this parameter should match to entried under `cluster` column from `expr` data frame.
#'
#' @param reference_id Same as `sample_id`, but for reference (right side of the enrichment plots i.e. tail end of the ranking)
#'
#' @param reference_cluster Same as sample_cluster but for reference (right side of the enrichment plots)
#'
#' @param gene_set A geneset (list object) by which gsea is calculated.  Can be one of 'hallmark'go', 'curated', 'immune', 'motif', 'all'. Or the user can provide a list object containing reference genesets.
#'
#' @param gseaParam fgsea parameter for changing the bar size in top-plots
#'
#' @param plot_individual String to plot enrichment results for a desired pathway. Must be an exact match to gene list
#'
#' @param append_title Add informative titles to the enrichment plot
#'
#' @param top_plots_title Add informative titles to summary enrichment plots
#'
#' @param save_png Export image as png
#'
#' @param png_units Png size units. Defults to 'in'
#'
#' @param png_width Png width. Defaults to 4
#'
#' @param png_height Png height. Defaults to 3
#'
#' @param append_to_filename Add custom string to png file
#'
#' @param verbose Logical. Tells you what's going on
#'
#' @param annot_text_color Color of the NES-p val annotations on graph. You can use `adjustcolor('red', alpha.f = 0)` to hide annotation text (by making it transparent)
#'
#' @param annot_text_size Size of the NES-p val annotations on graph
#'
#' @param annot_text_font Fontface of the annotation text (1,2,3,4, plain-bold-italic-bold and italic)
#'
#'
#'


gsea_replot <- function(res = NULL, ranked_genes= ranked_genes, pos_marker = NULL, neg_marker = NULL, sample_id = NULL, sample_cluster = NULL, reference_id = NULL, reference_cluster = NULL,
    gene_set = "hallmark", gseaParam = 1, plot_individual = NULL, append_title = T, top_plots_title = T, save_png = F, png_units = "in", png_width = 4,
    png_height = 3, append_to_filename = "", verbose = T, annot_text_color = "black", annot_text_size = 4, annot_text_fontface = 2) {


    # Read molecular signatures database (MSigDB) gene lists. files must be stored
    # if (gene_set == "hallmark") {
    #
    #     gene_set <- msigdb_hallmark
    #
    # }


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
            x_pos <- length(ranked_genes)/5

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

