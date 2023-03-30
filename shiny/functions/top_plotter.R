#' top_plotter
#'
#' @description A function to speed up plotting GSEA results. Use this function with previously generated 'gene rankings' and 'gsea results'
#'
#' @import data.table
#'
#' @export
#'
#' @param gsea_results A previously calculated data frame containing gsea results
#'
#' @param ranked_genes A named numeric vector of ranked genes.
#'
#' @param gene_set List object containing reference gene sets for analyses. Reference pathways correspond to name of the list item.
#'
#' @param top_n Specify the number of top pathways to be shown
#'
#' @param gseaParam Numeric value to increase the bar size of the genes in enrichment plot
#'
#' @param plot_title Add custom title to the summary plots. It can be useful for generating plots in a loop
#'
#' @param do.plot Logical. Set it to TRUE to create graphical outputs.
#'
#' @examples
#'
#' top_plotter(gsea_results = f480_d12_gsea, ranked_genes = f480_d12_rnk, gene_set = msig_HALLMARK, top_n = #' 5, gseaParam = 0.4)
#'
#' top_plotter(gsea_results = act_cd8_d12_gsea, ranked_genes = act_cd8_d12_rnk, gene_set = msig_HALLMARK, #' top_n = 5, gseaParam = 0.4)
#'
#'
#'


top_plotter <- function(gsea_results = NULL, ranked_genes = NULL, gene_set = NULL, top_n = 10, gseaParam = 1, plot_title = NULL, do.plot = T) {
    
    topPathwaysUp <- gsea_results[ES > 0][head(order(pval), n = top_n), pathway]
    
    topPathwaysDown <- gsea_results[ES < 0][head(order(pval), n = top_n), pathway]
    
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    
    plotGseaTable2(gene_set[topPathways], ranked_genes, gsea_results, gseaParam = gseaParam, plot_title = plot_title, do.plot = do.plot)
    
}


