#' plotGseaTable2
#'
#' @description This function plots GSEA results. The function was edited from fgsea package to enable to graph naming
#'
#' @param pathways Pathways to plot table, as in 'fgsea' function. List object containing reference gene sets for analyses. Reference pathways correspond to name of the list item.
#'
#' @param stats Gene-level stats, as in 'fgsea' function. A named numeric vector of ranked genes.
#'
#' @param fgseaRes Table with fgsea results. A previously calculated data frame containing gsea results
#'
#' @param gseaParam GSEA-like parameter. Adjusts displayed statistic values, values closer to 0 flatten plots. Default = 1, value of 0.5 is a good choice too. Numeric value to increase the bar size of the genes in enrichment plot.
#'
#' @param colwidths Vector of five elements corresponding to column width for grid.arrange. If column width is set to zero, the column is not drawn.
#'
#' @param plot_title Add custom title to the summary plots. It can be useful for generating plots in a loop
#'
#' @param do.plot Logical. Set it to TRUE to create graphical outputs.
#'
#' @import fgsea
#' @import gridExtra
#' @import ggplot2
#' @import grid
#'
#'
#'


plotGseaTable2 <- function(pathways, stats, fgseaRes, gseaParam = 1, colwidths = c(5, 3, 0.8, 1.2, 1.2), plot_title = NULL, do.plot = T) {

    environment(plotGseaTable2) <- environment(plotGseaTable)

    rnk <- rank(-stats)
    ord <- order(rnk)
    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
    statsAdj <- statsAdj/max(abs(statsAdj))
    pathways <- lapply(pathways, function(p) {
        unname(as.vector(na.omit(match(p, names(statsAdj)))))
    })
    ps <- lapply(names(pathways), function(pn) {
        p <- pathways[[pn]]
        annotation <- fgseaRes[match(pn, fgseaRes$pathway), ]
        list(textGrob(pn, just = "right", x = unit(0.95, "npc")), ggplot() + geom_segment(aes(x = p, xend = p, y = 0, yend = statsAdj[p]), size = 0.2) +
            scale_x_continuous(limits = c(0, length(statsAdj)), expand = c(0, 0)) + scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
            xlab(NULL) + ylab(NULL) + theme(panel.background = element_blank(), axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
            panel.grid = element_blank(), axis.title = element_blank(), plot.margin = rep(unit(0, "null"), 4), panel.spacing = rep(unit(0, "null"),
                4)), textGrob(sprintf("%.2f", annotation$NES)), textGrob(sprintf("%.1e", annotation$pval)), textGrob(sprintf("%.1e", annotation$padj)))
    })
    rankPlot <- ggplot() + geom_blank() + scale_x_continuous(limits = c(0, length(statsAdj)), expand = c(0, 0)) + scale_y_continuous(limits = c(-1,
        1), expand = c(0, 0)) + xlab(NULL) + ylab(NULL) + theme(panel.background = element_blank(), axis.line = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), panel.grid = element_blank(), axis.title = element_blank(), plot.margin = unit(c(0, 0, 0.5, 0), "npc"),
        panel.spacing = unit(c(0, 0, 0, 0), "npc"))
    grobs <- c(lapply(c("Pathway", "Gene ranks", "NES", "pval", "padj"), textGrob), unlist(ps, recursive = FALSE), list(nullGrob(), rankPlot,
        nullGrob(), nullGrob(), nullGrob()))
    grobsToDraw <- rep(colwidths != 0, length(grobs)/length(colwidths))
    plot_grob <- arrangeGrob(grobs = grobs[grobsToDraw], ncol = sum(colwidths != 0), widths = colwidths[colwidths != 0], top = plot_title)

    if (do.plot) {
        grid.draw(plot_grob)
    } else {
        plot_grob
    }

}
