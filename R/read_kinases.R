#' Read in KSP kinase prediction output
#'
#' @param path Path to the `*.*_kspoutput` file
#'
#' @return A matrix with of sites x (detected) kinases.
#' @export
#'
#' @examples
read_ksp <- function(path) {
  # Read file as single string
  file_chr <- readr::read_file(path) %>%
    # Split on separator line into one string per site
    stringr::str_split('\n?\\d* In terms of KSP, the top kinases with KSPScores are as follows:\n', simplify = TRUE) %>%
    stringr::str_trim() %>%
    # Split each site into individual kinase-score entries
    stringr::str_split('\n') %>%
    magrittr::extract(-1)

  # Map over sites, convert each kinase-score entry in each site into tidy data
  file_list <-
    # map(file_chr, ~.x %>% stringr::str_sub(start = 2) %>% stringr::str_split_fixed(' ', n = 2) %>% magrittr::set_colnames(c('kinase', 'score')))
    purrr::map_dfr(file_chr,
                   ~.x %>%
                     stringr::str_sub(start = 2) %>%            # Drop pre-kinase underscore
                     stringr::str_split_fixed(' ', n = 2) %>%   # Split kinase and score apart
                     {set_names(as.numeric(.[,2]), .[,1])}      # Turn into named numeric vector
                   ) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~replace(., is.na(.), 0.0))) # Turn NAs (non-predicted scores) into 0

  file_mat <- as.matrix(file_list)
  return(file_mat)

  # add non-detected kinases
  # attach site info somehow
  # normalise?
  #return(data_frame)
}

#' Read in NetPhorest prediction output
#'
#' Read in NetPhorest input. Keeps only kinase scores.
#'
#' @param path Path to the netphorest output file
#' @param return_long Optional boolean: whether to return the data in long format
#' (kinases on rows, like raw data), or whether to widen data (kinases as columns,
#' for use with U-CIE). Default is false (return wide)
#' @param cols_to_keep Columns to keep when pivoting to wide data.
#'   By default, keeps 'fasta_id', 'position', 'residue', and 'fragment_11',
#'   and drops 'method', 'organism', and 'binder_type'.
#'
#' @return A tibble with site/protein data and kinase scores
#' @export
#'
#' @examples
#' kinsub_netphorest_path <- system.file('extdata', 'kinsub_human_netphorest', package = 'phosphocie')
#' kinsub_netphorest <- read_netphorest(kinsub_netphorest_path)
#'
read_netphorest <- function(path, return_long = FALSE, cols_to_keep = c('fasta_id', 'position', 'residue', 'fragment_11')) {
  # Load in knowledge about netphorest output
  netphorest_colnames <- c('fasta_id', 'position', 'residue', 'fragment_11', 'method', 'organism', 'binder_type', 'kinase_fam', 'posterior', 'prior')
  # netphorest_known_kinases is included as internal package object and can be generated with data-raw/extract_netphorest_kinases.R

  # Read in data and keep only kinase predictions
  data <- readr::read_tsv(path, col_names = netphorest_colnames, skip = 1, show_col_types = FALSE)
  data <- dplyr::filter(data, binder_type == 'KIN')
  data <- dplyr::distinct(data, .keep_all = TRUE)

  # Option: return raw-ish long-format data
  if (return_long) return(data)

  # Check whether all kinases match known netphorest kinase families
  unknown_kinases <- dplyr::filter(data, !kinase_fam %in% netphorest_known_kinases)
  if (nrow(unknown_kinases) > 0) rlang::abort(glue::glue('Unknown kinases found: {paste(unknown_kinases$kinase_fam, collapse = ", ")}'))

  # Reshape data into wide format
  message('Reshaping data into wide matrix-like format, this might take a while.')
  data <- tidyr::pivot_wider(
      data,
      id_cols = dplyr::all_of(cols_to_keep),
      names_from = kinase_fam,
      values_from = posterior,
      values_fill = 0
    )

  # Append empty columns for any kinases not listed
  unpredicted_kinases <- netphorest_known_kinases[!netphorest_known_kinases %in% colnames(data)]
  if (!rlang::is_empty(unpredicted_kinases)) {
    empty_cols <- matrix(0, nrow = nrow(data), ncol = length(unpredicted_kinases)) %>%
      magrittr::set_colnames(unpredicted_kinases) %>%
      dplyr::as_tibble()
    data <- dplyr::bind_cols(data, empty_cols)
  }

  # Sort kinase columns
  data <- dplyr::relocate(data, dplyr::all_of(sort(netphorest_known_kinases)), .after = dplyr::last_col())

  return(data)
}


#' Remove off-target NetPhorest sites
#'
#' @description NetPhorest will always scan the entire sequence for possible sites, even when
#' the sequence is a short fragment with the desired site in the middle, like
#' when using [build_fastas()] to build NetPhorest input.
#' `filter_netphorest()` detects and removes these unwanted sites, in favour of the
#' site in the middle of the original sequence.
#'
#' Sites at the end of proteins are accounted for (see details), but sites at the
#' beginning (within `source_window_size`/2) cannot reliably be detected without
#' data on the original position in the protein (see `protein_pos_col`), due to
#' limitations in NetPhorest's output. Set `keep_uncertain` to handle these.
#'
#' `name_col`, `seq_col`, and `pos_col` defaults are based on the [read_netphorest()] default column names.
#'
#' @param data Data frame with data with one possible site per row (wide format),
#' presumably from [read_netphorest()].
#' @param source_window_size Length of the original site window sequences fed to NetPhorest.
#' @param name_col Name of column containing 'true' site names.
#' Data will be filtered down to one row in every group denoted by this column.
#' Default is 'fasta_id'.
#' @param seq_col Name of column containing netphorest-outputted sequences.
#'   Default is 'fragment_11'.
#' @param pos_col Name of column containing position of detected site in
#'   NetPhorest FASTA sequence. Default is 'position'.
#' @param fragment_col Optional: name of column containing sequence fragments
#' surrounding the true site, as generally extracted from the header column.
#' If not supplied and `match_fragments` is TRUE (both default), will attempt
#' to extract the fragments from `name_col` by extracting the last 3-9 letter/underscore/dash
#' characters from the fasta header.
#' @param match_fragments: Optional: Whether to use fragments to match the detected site
#' to the true site. Fragments can be any uneven length and are either extracted
#' from the name column or directly used from `fragment_col`. Default is TRUE.
#' @param keep_uncertain One of TRUE/FALSE/NULL (default). Some choices can be uncertain,
#' especially if there is no extra information from fragments (see details).
#' If `keep_uncertain` = TRUE, all possible uncertain values are kept;
#' if `keep_uncertain` = FALSE, all uncertain groups are fully dropped;
#' if `keep_uncertain` = NULL (default), a best guess is made based on proximity
#' to the midpoint of the sequence.
#'
#' @details This function will work with just the default netphorest output data
#' of `name_col`, `seq_col` and `pos_col`. However, filtering can be improved by providing
#' the section around the true site (fragment), either extracted by default (`fragment_col`
#' kept as NULL) or from a separate column (`fragment_col` supplied).
#'
#' @return The original dataset, without incorrect sites. Global row order is not preserved.
#' @export
#'
#' @examples
#' # Default usage
#' kinsub_netphorest_path <- system.file('extdata', 'kinsub_human_netphorest', package = 'phosphocie')
#' kinsub_netphorest <- read_netphorest(kinsub_netphorest_path)
#' kinsub_filtered <- filter_netphorest(kinsub_netphorest, source_window_size = 15)
#'
#' # Handle ambiguous sites
#' ambiguous_data <- data.frame(id = c("P13796|LCP1|L-plastin|S5|RGSVS", "P13796|LCP1|L-plastin|S5|RGSVS"),
#'                              pos = c(5, 7),
#'                              seq = c("-MARGsVSDEE", "ARGSVtDEEMM"))
#'
#'
#' ## Without further info, filter_phosphosite should pick 7
#' ## because it looks like an erroneous non-central site.
#' filter_netphorest(ambiguous_data,
#'                   name_col = 'id',
#'                   seq_col = 'seq',
#'                   pos_col = 'pos',
#'                   match_fragments = FALSE)
#'
#' ## Return all or none instead with `keep_uncertain`
#' filter_netphorest(ambiguous_data, 'id', 'seq', 'pos', match_fragments = FALSE, keep_uncertain = TRUE)
#' filter_netphorest(ambiguous_data, 'id', 'seq', 'pos', match_fragments = FALSE, keep_uncertain = FALSE)
#'
#' ## Or return the true value by integrating data from the fasta header, manually or automatic:
#' ambiguous_data_extra <- tidyr::extract(ambiguous_data, id, 'fragment', '\\|([A-Za-z_]{3,9})$',
#'                                        remove = FALSE, convert = TRUE)
#'
#' filter_netphorest(ambiguous_data_extra, 15, 'id', 'seq', 'pos')
#' filter_netphorest(ambiguous_data_extra, 15, 'id', 'seq', 'pos', fragment_col = 'fragment')
#'
filter_netphorest <- function(data,
                              name_col = 'fasta_id',
                              seq_col = 'fragment_11',
                              pos_col = 'position',
                              fragment_col = NULL,
                              match_fragments = TRUE,
                              source_window_size,
                              match_middle = TRUE,
                              keep_uncertain = NULL) {

  # Prep: check parameters
  submitted_cols <- c('name_col' = name_col,
                      'seq_col' = seq_col,
                      'pos_col' = pos_col)

  if (!all(submitted_cols %in% colnames(data))) {
    missing_cols <- submitted_cols[!submitted_cols %in% colnames(data)]
    rlang::abort(glue::glue(
      "Not all submitted column names are in the data.",
      "Faulty: {paste(names(missing_cols), missing_cols, sep = ' = ', collapse = ', ')}"
      ))
  }

  if (match_middle & missing(source_window_size)) {
    stop('In order to choose most middle site, source_window_size needs to be set.')
  }

  if (match_fragments) {
    # Extract fragment from name with default pattern if not provided separately
    if (is.null(fragment_col)) {
      data$temp_fragment <- stringr::str_extract(data[[name_col]], '[a-zA-Z_-]{3,9}$')
      fragment_col <- 'temp_fragment'

      if (anyNA(data$temp_fragment)) {
        nonextracted_headers <- data[[name_col]][is.na(data$temp_fragment)]
        if (length(nonextracted_headers > 50)) {
          nonextracted_headers <- c(nonextracted_headers[1:50], glue::glue('and {length(nonextracted_headers)-50} more...'))
        }
          stop(paste('Attempted to extract fragments from headers, but failed for headers ',
                     glue::glue('{paste(nonextracted_headers, collapse = ", ")}'),
                     '\nProvide fragments manually with fragment_col, or set match_fragments to FALSE.'))
      }
    }

    # Check whether fragments look like sequences
    if (any(stringr::str_detect('[^ARNDCEQGHILKMFPSTWYV_-]', toupper(data[[fragment_col]])))) {
      warning(paste('Attempted to extract fragments from headers, but non-standard AA residues detected.',
                    '\nCheck your data, provide fragment_col manually, or set match_fragments to FALSE.'))
    }

    # Make sure fragment has dashes for empty positions, like netphorest output
    data[[fragment_col]] <- stringr::str_replace_all(data[[fragment_col]], '_', '-')


    # Extract site surroundings to match against
    fragment_size = nchar(data[[fragment_col]])
    fragment_edge_size = (fragment_size - 1)/2

    data$temp_site_surroundings <- stringr::str_sub(
      data[[seq_col]],
      start = ceiling(nchar(data[[seq_col]])/2) - fragment_edge_size,
      end = ceiling(nchar(data[[seq_col]])/2) + fragment_edge_size
    )
  }

  # Prep: group and separate data
  data <- data %>%
    dplyr::group_by(.data[[name_col]])

  unique_data <- data %>%
    dplyr::filter(dplyr::n() == 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-dplyr::any_of(c('temp_fragment', 'temp_site_surroundings')))

  nonunique_data <- data %>%
    dplyr::filter(dplyr::n() > 1)



  if (match_fragments) {
  nonunique_data$step1 <- toupper(nonunique_data$temp_site_surroundings) == toupper(nonunique_data[[fragment_col]])
  unique_data <- nonunique_data %>%
    dplyr::filter(step1 & sum(step1) == 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-dplyr::any_of(c('temp_fragment', 'temp_site_surroundings', 'step1'))) %>%
    dplyr::bind_rows(unique_data, .)

  nonunique_data <- nonunique_data %>%
    dplyr::filter(sum(step1) != 1)
  }




  # Step 2: Find middle site.
  # For regular sites, just get middle site (ceiling(source_window_size/2), so site 8 for standard 15-width from phosphositeplus)
  if (match_middle) {
    nonunique_data$step2 = nonunique_data[[pos_col]] == ceiling(source_window_size/2)
    unique_data <- nonunique_data %>%
      dplyr::filter(step2 & sum(step2) == 1) %>%
      dplyr::select(-dplyr::starts_with('step')) %>%
      dplyr::ungroup() %>%
      dplyr::bind_rows(unique_data, .)

    nonunique_data <- nonunique_data %>% dplyr::filter(sum(step2) != 1)
  }

  # Handle inconclusive sites
  if (nrow(nonunique_data) > 0) {
    # Prepare warning text
    nonunique_names <- unique(nonunique_data[[name_col]])
    if (length(nonunique_names) > 50) {
      nonunique_names <- c(nonunique_names[1:50], glue::glue('and {length(nonunique_names)-50} more.'))
    }
    base_warn <- paste0("Could not conclusively determine original site for IDs ",
                        "{paste(nonunique_names, collapse = ', ')}.")
    extra_cols_warn <- ifelse(
      !match_fragments,
      paste0("\n(This generally happens for sites at the start of proteins. ",
             "Consider using `match_fragments` to reduce ambiguity.)"),
      "")

    # Filter, keep all, or discard data based on keep_uncertain
    if (is.null(keep_uncertain)) {
      uncertain_warn <- paste0("\nThe site closest to the center will be chosen. ",
                               "Change `keep_uncertain` to adjust this behaviour.")
      unique_data <- nonunique_data %>%
        dplyr::filter(abs(.data[[pos_col]] - ceiling(source_window_size/2)) == min(abs(.data[[pos_col]] - ceiling(source_window_size/2)))) %>%
        dplyr::select(-dplyr::starts_with('step')) %>%
        dplyr::ungroup() %>%
        dplyr::bind_rows(unique_data, .)
    } else if (isTRUE(keep_uncertain)) {
      uncertain_warn <- paste0("\nAll options for these groups will be retained. ",
                               "Change `keep_uncertain` to adjust this behaviour.")
      unique_data <- nonunique_data %>%
        dplyr::select(-dplyr::starts_with('step')) %>%
        dplyr::ungroup() %>%
        dplyr::bind_rows(unique_data, .)
    } else if (isFALSE(keep_uncertain)) {
      uncertain_warn <- paste0("\nThese groups will be fully discarded. ",
                               "Change `keep_uncertain` to adjust this behaviour.")
    }
    message(glue::glue(base_warn, uncertain_warn, extra_cols_warn))
  }

  return(unique_data)
}
