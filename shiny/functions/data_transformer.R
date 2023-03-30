#' Apply SCtransform to a list of Seurat objects
#'
#' @export
#'
#' @import sctransform
#'
#' @description This function is a wrapper to apply SCtransform to a list of Seurat objects.  See Seurat vignette for more details.
#'
#' @param object_list A list of Seurat objects (easily prepared by data_importer function in this package)
#'
#' @param verbose Whether to print messages and progress bar for `Seurat::SCTransform()`
#'
#' @param min_nFeature Minimum feature number for filtering out cells
#'
#' @param max_nFeature Maximum feature number for filtering out cells
#'
#' @param max_percent.mt Maximum mitochondrial gene percentage for filtering out cells
#'
#' @return SCtransformed object list
#'
#' @examples
#'
#' obj_list <- data_transformer(obj_list)
#'
#'
#'

data_transformer <- function(object_list,
                             verbose=T,
                             min_nFeature = 750,
                             max_nFeature = 5000,
                             max_percent.mt = 7.5){

    message("Loading sctransform package")
    suppressPackageStartupMessages(library(sctransform))


    options(future.globals.maxSize = 12000 * 1024^2)

    for(i in 1:length(object_list)){

        if(verbose) message(paste("Processing sample:", names(object_list)[i]))

        start_time <- Sys.time()

        object_list[[i]] <- subset(object_list[[i]], subset = nFeature_RNA > min_nFeature & nFeature_RNA < max_nFeature & percent.mt < max_percent.mt)

        object_list[[i]] <- SCTransform(object_list[[i]], vars.to.regress = "percent.mt", verbose = FALSE)

        end_time <- Sys.time()

        if(verbose) message(paste("Time elapsed:", round(difftime(end_time, start_time, units = "secs"), digits = 2)))
    }

    object_list


}




