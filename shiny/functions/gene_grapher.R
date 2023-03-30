#' gene_grapher
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr pull
#' @importFrom dplyr filter
#' @importFrom dplyr filter_
#' @importFrom dplyr select
#' @importFrom stringr str_replace_all
#' @importFrom gtools mixedsort
#' @importFrom purrr set_names
#'
#' @import ggpubr
#' @import ggplot2
#'

#' @import colorspace
#' @import gtools
#'
#' @description A function to generate gene expression plots from data frames derived from single cell RNAseq data.
#'
#' @param exprs An expression data frame. Easily generated using HTtools::df_extractor() function. Two columns should be present named "sample", and "cluster". The other columns contain gene names and normalized expression values
#'
#' @param genes_to_plot Character vector of gene names
#'
#' @param x_variable What to show on axis. It can be "sample" (default), "cluster"
#'
#' @param clusters_to_plot Character vector of clusters to include in the plots
#'
#' @param pos_marker Character vector of gene names for positively gating cells. Cells expressing zero counts will be excluded
#'
#' @param neg_marker Character vector of gene names for negatively gating cells. Cells expressing these genes (>0) will be excluded.
#'
#' @param plot_type The appearance of plots. Can be one of the following: "box", "violin", "bar"
#'
#' @param add_jitter Show transparent jitter points (logical)
#'
#' @param add_point Shot aligned points (instead of jitter)
#'
#' @param point_size Size of the jitter points
#'
#' @param point_alpha Transparency of the jitter points
#'
#' @param add_mean Logical. Set it to TRUE if you would like to show mean per group
#'
#' @param mean_color Color of the mean point
#'
#' @param add_median Logical. Set it to TRUE if you would like to median per group
#'
#' @param median_color Color of the median point
#'
#' @param sort_plots Logical. Set it to TRUE to alphabetically sort graphs based on gene name
#'
#' @param colors_to_use Character vector of colors per sample. Default is rainbow palette.
#'
#' @param show_stats Calculate statistical tests and display on graph
#'
#' @param comparisons A list of pairs of character vectors to show comparisons for. list(c(wt, ko1), c("wt", "ko2")). By default, all pairwise comparisons are calculated and plotted.
#'
#' @param stat_method Which statistical test to use. Select one of "wilcox.text" (default) or "t.test".
#'
#' @param p_adj_method How to adjust for multiple comparisons. Possible values are "holm" (default), "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")
#'
#' @param y_expand_low Expand low y-limit as a multiplicative factor (see mult argument of expand_scale())
#'
#' @param y_expand_high Expand high y-limit as a multiplicative factor (see mult argument of expand_scale())
#'
#' @param pval_y_offset Offset factor for determining y-coordinate of p values shown on the plot. It is a multiplicative factor y_max
#'
#' @param save_pdf Save results as a pdf file in the work directory
#'
#' @param append_to_filename String to append to the pdf filename for organizational purposes
#'
#' @param output_plot Logical. Set it to FALSE if you would like to suppress graphical output. Could be helpful for creating pdf files.
#'
#' @param assign_global_plotlist Logical. Creates a global plot_list object to allow further graphical manipulations
#'
#' @param show_progress Logical. Creates a message output to the console to let you know which gene is being analyzed and plotted
#'
#' @param image_columns Controls the number of figure columns on one output page (default is 4)
#'
#' @param image_rows Controls the number of figure rows on one output page (defaul is 4)
#'
#' @param ... other parameters to pass to gg
#'
#' @return Graphical outputs showing gene expression per sample or cluster
#'
#' @seealso [df_extractor()]
#'
#' @examples
#'
#' # Run function on exprs_df data frame
#' gene_grapher(exprs = exprs_df,
#' x_variable = "sample",
#' genes_to_plot = c("Ifng", "Gzmb", "Pdcd1"),
#' clusters_to_plot = "T cells",
#' pos_marker = NULL,
#' neg_marker = NULL ,
#' plot_type = "box",
#' add_jitter = T,
#' add_mean = T,
#' mean_color = "black",
#' add_median = F ,
#' sort_plots = T,
#' colors_to_use = cols,
#' show_stats = T,
#' comparisons = list(c("Sample1", "Sample2"),
#'                   c("Sample1", "Sample3")),
#' stat_method = "wilcox" ,
#' p_adj_method = "holm" ,
#' save_pdf = T,
#' append_to_filename = "_Tcell_genes" ,
#' output_plot = T,
#' assign_global_plotlist = F,
#' show_progress = F,
#' image_columns = 4,
#' image_rows = 4,
#' pval_y_offset = 4/5,
#' y_expand_high = 0.2)
#'
#'
#'
#'



gene_grapher <- function(exprs,
                         genes_to_plot,
                         x_variable = "sample",
                         clusters_to_plot = NULL,
                         pos_marker = NULL,
                         neg_marker = NULL,
                         plot_type = "box",
                         add_jitter = T,
                         add_point = F,
                         point_size = 0.2,
                         point_alpha = 0.2,
                         add_mean = T,
                         mean_color = "red",
                         add_median = T,
                         median_color = "blue",
                         sort_plots = F,
                         colors_to_use = NULL,
                         show_stats = T,
                         comparisons = NULL,
                         stat_method = "wilcox.test",
                         p_adj_method = "holm",
                         y_expand_low = 0,
                         y_expand_high = 0.2,
                         pval_y_offset = 5/6,
                         save_pdf = F,
                         append_to_filename = "",
                         output_plot = T,
                         assign_global_plotlist = F,
                         show_progress = F,
                         image_columns = 4,
                         image_rows = 4,
                         ...
)


{

    colnames(exprs) <- stringr::str_replace_all(colnames(exprs), "[:punct:]|[:space:]", "_")
    genes_to_plot <- stringr::str_replace_all(genes_to_plot, "[:punct:]|[:space:]", "_")

    gene_not_found <- genes_to_plot[!match(tolower(genes_to_plot), tolower(colnames(exprs)), nomatch = 0)]


    if(length(gene_not_found) > 0) {message(paste("Following genes are not found in data set: ", paste(gene_not_found, collapse = ", ")))}


    # # Match and correct capitalization of gene list
    genes_to_plot <- colnames(exprs)[tolower(colnames(exprs)) %in% tolower(genes_to_plot)]

    genes_to_plot <- unique(genes_to_plot)


    if(sort_plots ==T) {genes_to_plot <- gtools::mixedsort(genes_to_plot)}

    # Subset clusters of interest
    if(!is.null(clusters_to_plot)) {

        exprs <- dplyr::filter(exprs, cluster %in% clusters_to_plot)

    }

    # Subset cells that express markers of interest at any level (>0)
    if(!is.null(pos_marker)) {

        # pos_marker_split <- trimws(strsplit(pos_marker, split = ",")[[1]])
        # exprs <- filter_(exprs, paste(pos_marker_split , "!= 0", collapse = " & "))   # Can be used to pass a comma separated long string "gene1 , gene2"

        exprs <- dplyr::filter_(exprs, paste(pos_marker , "!= 0", collapse = " & "))

    }


    # Discard cells that express a marker gene which we want to negatively gate in our analyses
    if(!is.null(neg_marker)) {

        # neg_marker_split <- trimws(strsplit(neg_marker, split = ",")[[1]])
        # exprs <- filter_(exprs, paste(neg_marker_split , "== 0", collapse = " & "))  # Can be used to pass a comma separated long string "gene1 , gene2"

        exprs <- dplyr::filter_(exprs, paste(neg_marker , "== 0", collapse = " & "))

    }


    # Report which cells are being analyzed
    print((paste(dim(exprs)[1], "cells with the following annotations will be plotted")))
    print(paste("Samples:", paste(levels(droplevels(as.factor(exprs$sample))), collapse = ", ")))
    print(paste("Clusters:", paste(levels(droplevels(as.factor(exprs$cluster))), collapse = ", ")))


    omitted_genes <- character()

    for (gene in genes_to_plot){

        if(sum(exprs[,gene]) == 0){

            omitted_genes <- c(omitted_genes, gene)

        }

    }

    genes_to_plot <- genes_to_plot[!genes_to_plot %in% omitted_genes]

    if(length(genes_to_plot) == 0) stop("Requested genes have no expression in the subsetted data frame")

    if(length(omitted_genes) != 0) message(paste("Following genes are omitted due to zero expression: ", paste(omitted_genes, collapse = ", ")))

    plot_list <- list() # Initialize plot_list to store plots to be printed

    for(gene in genes_to_plot){

        if(show_progress == T) {message(paste("Plotting", gene))}

        if(show_stats == T){

            tryCatch(error = function(x){warning(paste("Statistics cannot be computed for ", gene))},
                     {

                         # Compute p-values
                         comparison_formula <- paste0(gene, "~", x_variable) %>%
                             as.formula()
                         stat_test <- ggpubr::compare_means(
                             comparison_formula,  data = exprs,
                             method = stat_method,
                             p.adjust.method = p_adj_method
                         )

                         # If a comparison list is provided, extract the comparisons of interest for plotting
                         if(!is.null(comparisons)){

                             stat_test <- do.call(rbind, comparisons) %>% # rbind to a matrix
                                 as.data.frame %>% # convert to a data.frame
                                 purrr::set_names(c("group1", "group2")) %>% # change the column names
                                 dplyr::inner_join(stat_test) # and inner join

                         }

                         # P-value y coordinates
                         y_max <- exprs %>%
                             dplyr::pull(gene) %>% max(na.rm = TRUE)

                         y_min <- exprs %>%
                             dplyr::pull(gene) %>% min(na.rm = TRUE)


                         p_value_y_coord <- rep(y_max*pval_y_offset, nrow(stat_test))

                         step_increase <- (1:nrow(stat_test))*(y_max/5)
                         p_value_y_coord <- p_value_y_coord + step_increase

                         stat_test <- stat_test %>%
                             dplyr::mutate(
                                 y.position =  p_value_y_coord,
                                 p.adj = format.pval(p.adj, digits = 1)
                             )

                     }) # close tryCatch

        } # close if(show_stats=T)

        if(is.null(colors_to_use)){

            group_number <- length(levels(as.factor(dplyr::pull(exprs, x_variable))))
            colors_to_use <- colorspace::rainbow_hcl(group_number)

        }

        if(grepl("^[0-9]", gene)) gene <- paste0("`", gene, "`")  # Append backticks for colnames starting with number

        if(plot_type == "box"){



            p <- ggpubr::ggboxplot(exprs, x= x_variable, y = gene, fill = x_variable, palette = colors_to_use,
                           title = gene, ylab="Expression", outlier.shape=NA, xlab = F, ...)+
                ggpubr::rremove("legend")+
                ggpubr::rotate_x_text(angle=45)

        } else if(plot_type == "bar"){



            p <- ggpubr::ggbarplot(exprs, x= x_variable, y = gene, fill = x_variable, palette = colors_to_use,
                           title = gene, ylab="Expression",
                           add =  "mean_se" , xlab = F, ...)+
                ggpubr::rremove("legend")+
                ggpubr::rotate_x_text(angle=45)



        } else if(plot_type == "violin"){



            p <- ggpubr::ggviolin(exprs, x= x_variable, y = gene, fill = x_variable, palette = colors_to_use,
                          title = gene, ylab="Expression",  scale="width", trim = T, ...)+ #draw_quantiles = c(0.25, 0.5, 0.75) can be added for marking IQR and median
                ggpubr::rremove("legend")+
                ggpubr::rotate_x_text(angle=45)


        } else {stop("Select one of the following as graph type: 'bar', 'box', 'violin'")}

        if(add_jitter == T){ p <- ggpubr::ggadd(p, add = "jitter", alpha=point_alpha, size = point_size) }

        if(add_point == T){ p <- ggpubr::ggadd(p, add = "point", alpha=point_alpha, size = point_size) }

        if(add_mean == T){ p <- ggpubr::ggadd(p, add = "mean", color = mean_color, size = 0.3)}

        if(add_median == T){ p <- ggpubr::ggadd(p, add = "median", color = median_color, size = 0.3)}

        if(show_stats == T){ tryCatch(error=function(x){},

                                      {

                                          p <- p +
                                              ggplot2::scale_y_continuous(expand = ggplot2::expand_scale(mult = c(y_expand_low, y_expand_high))) +
                                              ggpubr::stat_pvalue_manual(stat_test, label = "p.signif", size = 3.5)

                                      })
        }

        plot_list[[gene]] <- p


    }

    if(assign_global_plotlist == T) { assign("plot_list", plot_list, envir = .GlobalEnv)}

    if(output_plot ==T) print(ggpubr::ggarrange(plotlist = plot_list, ncol = image_columns, nrow = image_rows))

    if(save_pdf == T){

        arg_list <- list(cluster = clusters_to_plot, pos = paste(pos_marker, collapse = "."), neg = paste(neg_marker, collapse = "."))
        select_non_null <- !sapply(arg_list, function(x) {identical(x, "")})
        select_non_null2 <- !sapply(arg_list, is.null)
        select_non_null <- as.logical(select_non_null * select_non_null2)



        if(sum(select_non_null) == 0) {

            filename <- paste0("Unsubsetted_data_plots__", stat_method, ".pdf")

        } else {


            filename <- paste(arg_list[select_non_null], names(arg_list[select_non_null]), sep="_", collapse = "  ")
            filename <- paste0(filename,"__", stat_method, "_", append_to_filename, ".pdf")
        }

        filename <- gsub("\\\\", "", filename)

        ggpubr::ggexport(plotlist = plot_list, filename = filename,
                 nrow = image_rows, ncol=image_columns,
                 height = image_rows*3, width = image_columns*1.8, res = 300)

    }

}


