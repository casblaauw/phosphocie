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
#' @param path Path to the netphorest output file
#' @param return_long Optional boolean: whether to return the data in long format
#' (kinases on rows, like raw data), or whether to widen data (kinases as columns,
#' for use with U-CIE). Default is false (return wide)
#' @param split_fasta_header Optional boolean: whether to try and extract data
#' from the fasta header, like accession ID, gene name, and original site residue.
#' Supports UniProt FASTA headers (acc_id, protein) and [build_fastas()] headers
#' (acc_id, gene, protein, protein_res_col, protein_pos_col).
#'
#' @return A tibble with site/protein data and kinase scores
#' @export
#'
#' @examples
#' kinsub_netphorest_path <- system.file('extdata', 'kinsub_human_netphorest', package = 'phosphocie')
#' kinsub_netphorest <- read_netphorest(kinsub_netphorest_path)
#'
read_netphorest <- function(path, return_long = FALSE, split_fasta_header = FALSE) {
  # Load in knowledge about netphorest output
  netphorest_colnames <- c('fasta_id', 'position', 'residue', 'fragment_11', 'method', 'organism', 'binder_type', 'kinase_fam', 'posterior', 'prior')
  netphorest_kinase_names <- c("AMPK_group", "CDK2_CDK3_CDK1_CDK5_group", "CK1_group", "DMPK_group",
                               "EGFR_group", "InsR_group", "Src_group", "p38_group", "MAPK3_MAPK1_MAPK7_NLK_group",
                               "PKD_group", "CLK_group", "DAPK_group", "RCK_group", "LKB1",
                               "MST_group", "YSK_group", "PAK_group", "Pim3_Pim1_group", "Pim2",
                               "SLK_group", "ACTR2_ACTR2B_TGFbR2_group", "TLK_group", "MSN_group",
                               "p70S6K_group", "Eph_group", "NEK1_NEK5_NEK3_NEK4_NEK11_NEK2_group",
                               "ATM_ATR_group", "Abl_group", "AuroraA", "CDK4_CDK6_group", "CDK7",
                               "CK2_group", "CaMKII_group", "CaMKIV", "CaMKI_group", "DNAPK",
                               "EIF2AK2", "FLT3_CSF1R_Kit_PDGFR_group", "GRK_group", "GSK3_group",
                               "HIPK1_HIPK2_group", "IKKalpha_IKKbeta_group", "JAK2", "JNK_group",
                               "KDR_FLT1_group", "MAP2K_group", "Met_group", "PDHK_group", "PKA_group",
                               "PKB_group", "PKC_group", "PKGcGK_group", "ROCK_group", "RSK_group",
                               "SGK_group", "Syk_group", "TTK", "Tec_group", "Trk_group", "Tyk2"
  )

  # Read in data and keep only kinase predictions
  long_data <- readr::read_tsv(path, col_names = netphorest_colnames, skip = 1, show_col_types = FALSE) %>%
    dplyr::filter(binder_type == 'KIN')

  # Option: split fasta header into constituent components
  if (split_fasta_header) {
    # Detect fasta header type
    ## Detect uniprot fasta if it starts with a database marker
    if (mean(stringr::str_detect(head(long_data$fasta_id, 100), 'sp|tr\\|')) > 0.75) { # 75+% of data matches pattern
      long_data <- long_data %>%
        tidyr::separate(fasta_id, c(NA, 'acc_id', 'uniprot_name'), sep = '\\|', remove = FALSE) %>%
        tidyr::separate(uniprot_name, c('protein', NA), sep = '_')
      ## Detect own fasta if it starts with a uniprot ID (6-10 word characters)
    } else if (mean(stringr::str_detect(head(long_data$fasta_id, 100), '^\\w{6,10}\\|')) > 0.75) {
      long_data <- long_data %>%
        tidyr::separate(fasta_id, c('acc_id', 'gene', 'protein', 'orig_ptm_residue'), sep = '\\|', remove = FALSE) %>%
        tidyr::separate(orig_ptm_residue, c('orig_res', 'orig_pos'), sep = 1)
    } else {
      rlang::abort(glue::glue('Could not detect fasta header format. Set split_fasta_header to FALSE to disable splitting. Header example: {long_data$fasta_id[1]}'))
    }
  }

  # Option: return raw-ish long-format data
  if (return_long) return(long_data)

  # Check whether all kinases match known netphorest kinase families
  unknown_kinases <- dplyr::filter(long_data, !kinase_fam %in% netphorest_kinase_names)
  if (nrow(unknown_kinases) > 0) rlang::abort(glue::glue('Unknown kinases found: {paste(unknown_kinases$kinase_fam, collapse = ", ")}'))

  # Reshape data into wide format
  print('Reshaping data into wide matrix-like format, this might take a while.')
  wide_data <- long_data %>%
    tidyr::pivot_wider(
      id_cols = c(dplyr::all_of(c('fasta_id', 'position', 'residue', 'fragment_11')), dplyr::any_of(c('acc_id', 'gene', 'protein', 'orig_ptm_residue', 'orig_res', 'orig_pos'))),
      names_from = kinase_fam,
      values_from = posterior,
      values_fill = 0
    )

  # Append empty columns for any kinases not listed
  unpredicted_kinases <- netphorest_kinase_names[!netphorest_kinase_names %in% unique(long_data$kinase_fam)]
  if (!rlang::is_empty(unpredicted_kinases)) {
    empty_cols <- matrix(0, nrow = nrow(wide_data), ncol = length(unpredicted_kinases)) %>%
      magrittr::set_colnames(unpredicted_kinases) %>%
      dplyr::as_tibble()
    wide_data <- dplyr::bind_cols(wide_data, empty_cols)
  }

  # Sort kinase columns
  wide_data <- wide_data %>%
    dplyr::relocate(dplyr::all_of(sort(netphorest_kinase_names)), .after = dplyr::last_col())

  return(wide_data)
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
#' @param protein_res_col Optional: name of column containing residue of true site
#'   in NetPhorest FASTA sequence. If provided, used to filter out obvious false
#'   detected sites by comparing with detected residue (see details).
#' @param detected_res_col Optional: name of column containing residue of detected site
#'   in NetPhorest FASTA sequence. If provided, used to filter out obvious false
#'   detected sites by checking against true residue from `protein_res_col` (see details).
#'   Requires `protein_res_col` too. If missing, inferred from `seq_col`.
#' @param protein_pos_col Optional: name of column containing original position
#'   in the protein. If provided, used to filter sites at start of protein (see details)
#' @param netphorest_window_size Width of sequence windows in seq_col.
#'   If missing, inferred from `seq_col` widths.
#' @param keep_uncertain One of TRUE/FALSE/NULL. Some choices can be uncertain,
#' especially if there is no extra information from `protein_pos_col` or `protein_res_col`,
#' (see details).
#' If `keep_uncertain` = TRUE, all possible uncertain values are kept;
#' if `keep_uncertain` = FALSE, all uncertain groups are fully dropped;
#' if `keep_uncertain` = NULL (default), a best guess is made based on proximity
#' to the midpoint of the sequence.
#'
#' @details This function will work with just the default netphorest output data
#' of `name_col`, `seq_col` and `pos_col`. However, filtering can be improved by providing
#' the true residue of the site (`protein_res_col`) (and optionally the detected site residue
#' `detected_res_col` manually) and/or the original position in the full protein (`protein_pos_col`).
#' \describe{
#'   \item{`protein_res_col`/`detected_res_col`}{The columns from these two arguments are
#'    compared against eachother, throwing out all detected sites that aren't the true residue.}
#'   \item{`protein_pos_col`}{The column from this argument is compared to the intra-sequence
#'    position from `pos_col`. If they match, the true site is from the start of the protein
#'    and the proper detected site can be matched. Otherwise, these early sites will be
#'    indistinguishable from erroneous detected sites at the left edge of the netphorest window,
#'    which are generally discarded (as we expect the site to be in the middle of the netphorest window).}
#' }
#'
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
#' ambiguous_data <- data.frame(id = c("P13796|LCP1|L-plastin|S5", "P13796|LCP1|L-plastin|S5"),
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
#'                   source_window_size = 15)
#'
#' ## Return all or none instead with `keep_uncertain`
#' filter_netphorest(ambiguous_data, 15, 'id', 'seq', 'pos', keep_uncertain = TRUE)
#' filter_netphorest(ambiguous_data, 15, 'id', 'seq', 'pos', keep_uncertain = FALSE)
#'
#' ## Or return the true value by integrating outside site data, like extracted from the fasta header:
#' ambiguous_data_extra <- tidyr::extract(ambiguous_data, id, c('orig_site', 'orig_res'), '\\|(\\w)(\\d{1,4})$',
#'                                        remove = FALSE, convert = TRUE)
#'
#' filter_netphorest(ambiguous_data_extra, 15, 'id', 'seq', 'pos', protein_res_col = 'orig_res')
#' filter_netphorest(ambiguous_data_extra, 15, 'id', 'seq', 'pos', protein_pos_col = 'orig_site')
#'
filter_netphorest <- function(data,
                              source_window_size,
                              name_col = 'fasta_id',
                              seq_col = 'fragment_11',
                              pos_col = 'position',
                              protein_res_col = NULL,
                              detected_res_col = NULL,
                              protein_pos_col = NULL,
                              netphorest_window_size,
                              keep_uncertain = NULL) {

  # Prep: check parameters
  submitted_cols <- c('name_col' = name_col,
                      'seq_col' = seq_col,
                      'pos_col' = pos_col,
                      'protein_res_col' = protein_res_col,
                      'detected_res_col' = detected_res_col,
                      'protein_pos_col' = protein_pos_col)

  if (!all(submitted_cols %in% colnames(data))) {
    missing_cols <- submitted_cols[!submitted_cols %in% colnames(data)]
    rlang::abort(glue::glue(
      "Not all submitted column names are in the data.",
      "Faulty: {paste(names(missing_cols), missing_cols, sep = ' = ', collapse = ', ')}"
      ))
  }

  if (missing(netphorest_window_size)) {
    netphorest_window_size <- round(mean(nchar(data[[seq_col]])), 0)
  }

  # Optional step 1: Match proposed site residues to true site residues
  if (!is.null(protein_res_col)) {
    if (is.null(detected_res_col)) {
      detected_res_col <- 'res_col_extract'
      data <- dplyr::mutate(data, res_col_extract = toupper(
        stringr::str_sub(.data[[seq_col]],
                         start = ceiling(nchar(.data[[seq_col]])/2),
                         end = ceiling(nchar(.data[[seq_col]])/2))))}
    data <- dplyr::filter(data, toupper(.data[[detected_res_col]]) == toupper(.data[[protein_res_col]]))
    print(data)
  }

  # Prep: group and separate data
  data <- data %>%
    dplyr::group_by(.data[[name_col]])

  unique_data <- data %>%
    dplyr::filter(dplyr::n() == 1) %>%
    dplyr::ungroup()

  nonunique_data <- data %>%
    dplyr::filter(dplyr::n() > 1)

  # Optional step 2: For start-truncated sites, match full position to res position
  if (!is.null(protein_pos_col)) {
    step1 <- nonunique_data %>%
      dplyr::mutate(step1 = .data[[pos_col]] == .data[[protein_pos_col]])

    unique_data <- step1 %>%
      dplyr::filter(step1) %>%
      dplyr::ungroup() %>%
      dplyr::select(-dplyr::starts_with('step')) %>%
      dplyr::bind_rows(unique_data, .)

    nonunique_data <- step1 %>%
      dplyr::filter(!any(step1))
  }


  # Step 3/4: Find middle site.
  # For end-truncated sites (i.e shorter sequences), use expected # of dashes to shift midpoint backwards
  # For regular sites, just get middle site (ceiling(source_window_size/2), so site 8 for standard 15-width from phosphositeplus)

  # If small view window extends beyond big window, expect filling dashes:
  # expected_n_dashes = position + half_window_width - source_window_size_size
  # Dashes cannot be negative, so
  # expected_n_dashes = max(0, position + half_netphorest_window_size_width - source_window_size_size)
  # If original sequence (source_window_size) was less than 15 seq, there will be more dashes than expected
  # So shift the midpoint by difference between expected and observed:
  # new_mid = midpoint - (n_dash - max(0, position + half_window_width - source_window_size_size))
  step2 <- nonunique_data %>%
    dplyr::mutate(
      step2_shift = dplyr::if_else(
        stringr::str_ends(.data[[seq_col]], '-'),
        stringr::str_count(.data[[seq_col]], '-') - max(c(0, .data[[pos_col]] + floor(netphorest_window_size/2) - source_window_size)),
        0
      ),
      step2_mid = ceiling(source_window_size/2) - .data[['step2_shift']],
      step2 = .data[[pos_col]] == .data[['step2_mid']]
    )
  unique_data <- step2 %>%
    dplyr::filter(step2) %>%
    dplyr::select(-dplyr::starts_with('step')) %>%
    dplyr::ungroup() %>%
    dplyr::bind_rows(unique_data, .)

  nonunique_data <- step2 %>% dplyr::filter(!any(step2))

  # Handle inconclusive sites
  if (nrow(nonunique_data) > 0) {
    # Prepare warning text
    base_warn <- paste0("Could not conclusively determine original site for IDs ",
                        "{paste(unique(nonunique_data[[name_col]]), sep = ', ')}.")
    extra_cols_warn <- ifelse(
      is.null(protein_res_col) | is.null(protein_pos_col),
      paste0("\n(This generally happens for sites at the start of proteins. ",
             "Consider providing `protein_res_col`, `detected_res_col` and/or `protein_pos_col`",
             "to reduce ambiguity.)"),
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
      uncertain_warn <- paste0("\nThese groups will be fully discarded.",
                               "Change `keep_uncertain` to adjust this behaviour.")
    }
    message(glue::glue(base_warn, uncertain_warn, extra_cols_warn))
  }

  return(unique_data)
}
