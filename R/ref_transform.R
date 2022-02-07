#' Reference transformation scalars
#'
#' A numeric vector for use with [kinase2cielab()], containing the transformations to
#' go from a 3D representation made with [uwot::umap_transform(your_data, ref_umap)]
#' to colours in the CIELAB coordinate space.
#' Based on the UMAP->UCIE transformation of the PhosphoSitePlus phosphorylation
#' sites dataset, which contains almost 240000 phosphorylation sites.
#' Created using the code in `data_raw/create_reference_colourspace.R`.
#'
#'
#' @format A numeric vector, containing the 7 values needed to transform:
#' \describe{
#'   \item{S}{The scaling parameter}
#'   \item{RL}{The rotation parameter around the L CIELAB axis}
#'   \item{Ra}{The rotation parameter around the a CIELAB axis}
#'   \item{Rb}{The rotation parameter around the b CIELAB axis}
#'   \item{TL}{The translation parameter on the L CIELAB axis}
#'   \item{Ta}{The translation parameter on the a CIELAB axis}
#'   \item{Tb}{The translation parameter on the b CIELAB axis}
#' }
#' @source \url{https://www.phosphosite.org/staticDownloads}
"ref_transform"
