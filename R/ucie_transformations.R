#' Calculate UCIE geometric transformation values
#'
#' Scales, rotates and transforms the data to best fit in the CIELAB colourspace.
#' Wrapper around ucie:::FitColorsFunction()
#'
#' @param dataset A data frame or matrix of 3D values.
#' @param rownames_col Optional: name of column to ignore during calculations.
#'
#' @details `dataset` can be a 2D or 3D matrix or data frame, with optional
#'   rownames in the first column or a `rownames_col`-defined column.
#'
#' @return Returns a numeric vector with a scaling value, 3 rotation values, and 3 translation values for use in [kinase2cielab()].
#' @export
#'
#' @examples
ucie_transformations <- function(dataset, rownames_col = NULL) {
  prep_ucie_data(dataset, rownames_col = rownames_col)

  movement <- ucie:::FitColorsFunction(dataset, WL = 1, Wa = 1, Wb = 1) %>%
    stats::setNames(c('S', 'RotL', 'Rota', 'Rotb', 'TrL', 'Tra', 'Trb'))
  return(movement)
}

#' Title
#'
#' @param dataset A data frame or matrix of 3D values.
#' @param transform_vals A numeric vector with 1 scaling, 3 rotation, and 3 translation values,
#' as calculated by [ucie_transformations()].
#' @param LAB_coordinates Optional: boolean, whether to return CIELAB coordinates instead of RGB codes.
#' @param rownames_col Optional: name of column to ignore during calculations.
#'
#' @details `dataset` can be a 2D or 3D matrix or data frame, with optional
#'   rownames in the first column or a `rownames_col`-defined column.
#'
#' @return
#' @export
#'
#' @examples
kinase2cielab <- function(dataset, transform_vals, LAB_coordinates = FALSE, rownames_col = NULL) {

  dataset <- prep_ucie_data(dataset, rownames_col = rownames_col)

  dataset <- ucie:::Scaling(dataset, transform_vals[1]*1)
  dataset <- ucie:::Rotation(as.matrix(dataset), transform_vals[2], transform_vals[3], transform_vals[4])
  dataset <- ucie:::Translation(as.matrix(dataset), transform_vals[5], transform_vals[6], transform_vals[7])

  Lab <- dataset
  Lab <- round(Lab, 2)
  rawdata = structure(
    list(
      Lstar = c(Lab[, 1]),
      Astar = c(Lab[, 2]),
      Bstar = c(Lab[, 3])
    ),
    .Names = c("Lstar", "Astar", "Bstar"),
    row.names = c(rownames(dataset)),
    class = "data.frame"
  )

  LABdata <- with(rawdata, colorspace::LAB(Lstar, Astar, Bstar))

  if(LAB_coordinates==FALSE){
    colors <- as.data.frame(cbind(rownames(dataset),colorspace::hex(LABdata, fix = TRUE)))
  } else {
    colors <- as.data.frame(cbind(rownames(dataset),as.data.frame(LABdata@coords)))
    colnames(colors) <- c("names", "L", "a", "b")
  }

  return(colors)
}

#' Check and reformat data for use with U-CIE functions
#'
#' Takes a 2D or 3D matrix or data frame and returns a 3D data frame, checked for non-numeric values.
#'
#' @param dataset A data frame or matrix of 2D or 3D values.
#' @param rownames_col Optional: name of column to ignore during calculations.
#'
#' @details
#'   The function does a number of processing steps to accept as many data formats as feasible:
#'   * Rownames in the first column of the matrix (as tested by unconvertability to numeric) or of the data frame are discarded.
#'   * Character matrices are converted to numeric.
#'   * Matrices are transformed to data frames.
#'   * If rownames_col is given, that column is discarded as well.
#'   * 2D matrices or data frames are expanded to 3D by padding with 1's.
#' @keywords internal
#'
#' @examples
prep_ucie_data <- function(dataset, rownames_col) {

  # If matrix, transform to data frame
  if (!inherits(dataset, "data.frame")) {
    # If rownames encoded in first matrix row, move to true rownames
    if (anyNA(as.numeric(dataset[,1]))) {
      val <- dataset[,1]
      dataset <- dataset[,-1]
      rownames(dataset) <- val
    }
    # If character matrix, transform into numeric
    if (!all(is.numeric(dataset))) {
      class(dataset) <- 'numeric'
    }
    # Check for missing/leftover character data
    if (anyNA(dataset)) {
      rlang::abort(paste("Your matrix data contains missing data or",
                         "non-numeric data outside the first column.",
                         "Please check your dataset."))
    }
    dataset <- as.data.frame(dataset)
    warning("The dataset has been transformed into a data frame.")
  }

  # If rownames column supplied or first column looks like rownames, move to actual rownames
  if (!is.null(rownames_col)) {
    rownames_temp <- dataset[[rownames_col]]
    rownames(dataset) <- rownames_temp
    dataset <- dataset[colnames(dataset) != rownames_col]
  } else if (is.character(dataset[,1])) {
    rownames(dataset) <- dataset[,1]
    dataset <- dataset[-1]
  }

  # Check for non-numeric data
  if (!all(purrr::map_lgl(dataset, is.numeric))) {
    rlang::abort(paste("The dataset contains non-numeric columns.",
                       "Please check your data and use the rownames_col argument if necessary."))
  }

  # Check data dimensionality
  if (ncol(dataset) == 2) {
    warning("Data expanded to 3D!")
    dataset <- cbind(dataset, rep(1, nrow(dataset)))
  } else if (ncol(dataset) > 3) {
    rlang::abort("The dataset should have 3 numeric columns!")
  }

  # Check for missing values
  if (anyNA(dataset)) {
    rlang::abort("The dataset has missing values. Check again!")
  }

  return(dataset)
}
