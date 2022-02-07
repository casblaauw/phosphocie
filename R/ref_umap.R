#' `uwot::umap` pre-trained PhosphoSitePlus reference model
#'
#' An object for use with [uwot::umap_transform()], containing the 3D
#' representation of the PhosphoSitePlus phosphorylation sites dataset,
#' which contains almost 240000 phosphorylation sites. Originally trained
#' with spread = 5, using the code in `data_raw/create_reference_colourspace.R`.
#'
#'
#' @format A list, containing the elements of an UMAP model:
#' \describe{
#'   \item{embedding}{3D coordinates, as a (239419, 3) matrix}
#'   \item{...}{parameters of the model}
#' }
#' @source \url{https://www.phosphosite.org/staticDownloads}
"ref_umap"
