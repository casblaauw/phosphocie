#' Calculate UCIE geometric transformation values
#'
#' Scales, rotates and transforms the data to best fit in the CIELAB colourspace.
#' Wrapper around ucie:::FitColorsFunction()
#'
#' @param dataset A data frame or matrix of 3D values.
#' @param rownames_col Optional: name of column to ignore during calculations.
#' @param center A boolean indicating whether data should be centered during the fitting
#' process. Improves performance, but requires `center = TRUE` in [transform_data()] too.
#'
#' @details `dataset` can be a 2D or 3D matrix or data frame, with optional
#'   rownames in the first column or a `rownames_col`-defined column.
#'
#' @return Returns a numeric vector with a scaling value (S), 3 rotation values
#' (RotL, Rota, Rotb), and 3 translation values (TrL, Tra, Trb) for use in [transform_data()].
#' @export
#'
#' @examples
fit_transform <- function(dataset, rownames_col = NULL, center = TRUE) {
  # Prepare dataset
  dataset <- prep_ucie_data(dataset, rownames_col = rownames_col)

  # Fit all parameters from untransformed data
  fit_params <- ucie:::FitColorsFunction(dataset, WL = 1, Wa = 1, Wb = 1, center = center)
  return(fit_params)
}

#' Get colours from fitted kinase data
#'
#' @param dataset A data frame or matrix of 3D values.
#' @param transform_vals A numeric vector with 1 scaling, 3 rotation, and 3 translation values,
#' as calculated by [fit_transform()].
#' @param LAB_coordinates Optional: boolean, whether to return CIELAB coordinates instead of RGB codes.
#' @param rownames_col Optional: name of column to ignore during calculations.
#' @param fix Optional: boolean, whether to force points outside the colour space inside (TRUE)
#' or return NA (FALSE). Default is TRUE.
#'
#' @details `dataset` can be a 2D or 3D matrix or data frame, with optional
#'   rownames in the first column or a `rownames_col`-defined column.
#'
#' @return Returns a data frame with names and hex codes in `name` and `colour`
#' (or `name`/`L`/`a`/`b` if `LAB_coordinates = TRUE`).
#' @export
#'
#' @examples
transform_data <- function(dataset, transform_vals, LAB_coordinates = FALSE, rownames_col = NULL, fix = TRUE, center = TRUE) {

  # Prep dataset
  dataset <- prep_ucie_data(dataset, rownames_col = rownames_col)

  # Transform data
  dataset <- ucie:::Scaling(dataset, transform_vals[1])
  dataset <- ucie:::Rotation(dataset, transform_vals[2], transform_vals[3], transform_vals[4])
  dataset <- ucie:::Translation(dataset, transform_vals[5], transform_vals[6], transform_vals[7])

  if (center) {
    data_centroid <- ucie:::PolygonCentroid(dataset)
    cielab_centroid <- ucie:::PolygonCentroid(ucie:::RGB2Lab(ucie:::RGB_space))
    translation_set <- cielab_centroid - data_centroid
    dataset <- ucie:::Translation(dataset, translation_set[1], translation_set[2], translation_set[3])
  }

  colorspace_obj <- colorspace::LAB(round(dataset, 2))

  if (LAB_coordinates) {
    # Turn coordinates matrix into a data frame
    col_coords <- colorspace::coords(colorspace_obj) %>%
      as.data.frame() %>%
      tibble::rownames_to_column('name')
    return(col_coords)
  } else {
    # Turn hex character vector into a data frame
    col_hex <- colorspace::hex(colorspace_obj, fixup = fix) %>%
      data.frame(colour = .) %>%
      tibble::rownames_to_column('name')
    return(col_hex)
  }

  return(colors)
}

#' Get kinase score colours based on reference colourspace
#'
#' @param data A data frame or matrix, with the same columns as the reference
#' kinase scoring dataset (NetPhorest, 61 kinases), and optionally a rownames
#' column (see `rownames_col`)
#' @param k Number of neighbours to consider.
#' @param rownames_col Optional: name of a column containing rownames, set as
#' true rownames before continuing. Allows usage of tidy data in this matrix-based function.
#' @param fix Optional: boolean, whether to force points outside the colour space inside (TRUE)
#' or return NA (FALSE). Default is TRUE.
#'
#' @return A data frame with columns `name` and `colour`. Use colours with [ggplot2::scale_color_identity()].
#' @export
#'
#' @examples
kinase2cielab <- function(data, k = 10, LAB_coordinates = FALSE, rownames_col = NULL, fix = TRUE) {

  # Check and prepare data
  data <- prep_ucie_data(data, rownames_col = rownames_col, check_3D = FALSE)
  if (ncol(data) != ncol(ref_kinase)) {
    rlang::abort(paste('Your data must have the same amount of columns as the reference.',
                       glue::glue('Data columns: {ncol(data)}, ref columns: {ncol(ref_kinase)}')
                       ))
  }

  # Separate tyrosines
  tyrosine_kinases <- c('Eph_group', 'Src_group', 'Met_group', 'EGFR_group', 'FLT3_CSF1R_Kit_PDGFR_group', 'Tec_group', 'KDR_FLT1_group', 'InsR_group', 'Abl_group', 'JAK2', 'Trk_group', 'Syk_group', 'Tyk2')
  serine_threonine_kinases <- netphorest_known_kinases[!netphorest_known_kinases %in% tyrosine_kinases]
  tyrosine_scores <- rowMeans(data[,tyrosine_kinases])
  serine_threonine_scores <- rowMeans(data[,serine_threonine_kinases])
  tyrosine_indices <- which(tyrosine_scores > serine_threonine_scores)
  serine_threonine_indices <- which(serine_threonine_scores > tyrosine_scores)

  # Map nearest neighbours
  indices <- FNN::get.knnx(ref_kinase, data[serine_threonine_indices, ], k = k)$nn.index

  # Get and average colours of nearest neighbours
  if (LAB_coordinates) {
    colour_coords <- matrix(0, nrow = nrow(data), ncol = 3) # 0,0,0 is #000000 in UCIE space
    colnames(colour_coords) <- c('L', 'A', 'B')
    colour_coords[serine_threonine_indices,] <- t(apply(indices, 1, function(indices_vec) colMeans(ref_pca[indices_vec,])))
    colour_coords <- as.data.frame(colour_coords)
    if (!is.null(rownames(data))) {rownames(colour_coords) <- rownames(data)}
    colour_coords <- tibble::rownames_to_column(colour_coords, 'name')
    return(colour_coords)
  } else {
    colours <- rep('#000000', nrow(data))
    colour_coords <- t(apply(indices, 1, function(indices_vec) colMeans(ref_pca[indices_vec,])))
    colours[serine_threonine_indices] <- colorspace::hex(colorspace::LAB(colour_coords), fixup = TRUE)
    colours <- data.frame(colour = colours)
    if (!is.null(rownames(data))) {rownames(colours) <- rownames(data)}
    colours <- tibble::rownames_to_column(colours, 'name')
    return(colours)
  }
}


#' Check and reformat data for use with U-CIE functions
#'
#' Takes a 2D or 3D matrix or data frame and returns a 3D matrix, checked for non-numeric values.
#'
#' @param dataset A data frame or matrix of 2D or 3D values.
#' @param rownames_col Optional: name of column to ignore during calculations.
#' @param check_3D Optional: whether to check and force 3D output. Needed for reduced/cielab data,
#' not needed for KNN data.
#'
#' @details
#'   The function does a number of processing steps to accept as many data formats as feasible:
#'   * Rownames in the first column of the matrix (as tested by unconvertability to numeric) or of the data frame are discarded.
#'   * Character matrices are converted to numeric.
#'   * Data frames are converted to matrices.
#'   * If rownames_col is given, that column is discarded as well.
#'   * 2D matrices or data frames are expanded to 3D by padding with 1's.
#'
#' @keywords internal
prep_ucie_data <- function(dataset, rownames_col = NULL, check_3D = TRUE) {

  if (!is.null(rownames_col)) {
    if (ncol(dataset) < 3) {
      rlang::abort('Your data needs to be at least 2D to use `rownames_col`.')
    }
    if (inherits(dataset, "data.frame")) {
      dataset <- tibble::column_to_rownames(dataset, rownames_col)
    } else if (inherits(dataset, 'matrix')) {
      rownames_temp <- dataset[,rownames_col]
      dataset <- dataset[,which(colnames(dataset) != rownames_col), drop = FALSE]
      rownames(dataset) <- unname(rownames_temp)
    }
  }

  # If data frame, transform to matrix
  if (inherits(dataset, "data.frame")) {
    # If rownames column supplied or first column looks like rownames, move to actual rownames
    if (is.character(dataset[[1]]) & is.null(rownames_col)) {
      rownames(dataset) <- dataset[[1]]
      dataset <- dataset[[-1]]
    }

    # Check for non-numeric data
    if (!all(purrr::map_lgl(dataset, is.numeric))) {
      rlang::abort(paste("The dataset contains non-numeric columns.",
                         "Please check your data and use the rownames_col argument if necessary."))
    }
    dataset <- as.matrix(dataset)
  }

  # If rownames encoded in first matrix row, move to true rownames if not present or drop
  if (anyNA(as.numeric(dataset[,1]))) {
    val <- dataset[,1]
    dataset <- dataset[,-1]
    if (is.null(rownames(dataset))) {
      rownames(dataset) <- val
    }
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

  if (check_3D) {

    # Check data dimensionality
    if (ncol(dataset) == 2) {
      warning("Data expanded to 3D!")
      dataset <- cbind(dataset, 1)
    }
    if (ncol(dataset) != 3) {
      rlang::abort("The dataset should have 3 numeric columns!")
    }

    # Check for missing values
    if (anyNA(dataset)) {
      rlang::abort("The dataset has missing values. Check again!")
    }

    # Set column names
    colnames(dataset) <- c('X', 'Y', 'Z')
  }

  return(dataset)
}


