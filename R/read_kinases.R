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
    dplyr::mutate(dplyr::across(everything(), ~replace(., is.na(.), 0.0))) # Turn NAs (non-predicted scores) into 0

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
#' for use with U-CIE)
#'
#' @return A tibble with site/protein data and kinase scores
#' @export
#'
#' @examples
#' kinsub_netphorest_path <- system.file('extdata', 'kinsub_human_netphorest', package = 'phosphocie')
#' kinsub_netphorest <- read_netphorest(kinsub_netphorest_path)
#'
read_netphorest <- function(path, return_long = FALSE) {
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

  # Detect fasta header type

  ## Detect uniprot fasta if it starts with a database marker
  if (mean(stringr::str_detect(head(long_data$fasta_id, 100), 'sp|tr\\|')) > 0.75) { # 75+% of data matches pattern
    long_data <- long_data %>%
      tidyr::separate(fasta_id, c(NA, 'acc_id', 'uniprot_name'), sep = '\\|', remove = FALSE) %>%
      tidyr::separate(uniprot_name, c('protein', NA), sep = '_')
  ## Detect own fasta if it starts with a uniprot ID (6-10 word characters)
  } else if (mean(stringr::str_detect(head(long_data$fasta_id, 100), '^\\w{6,10}\\|')) > 0.75) {
    long_data <- tidyr::separate(long_data, fasta_id, c('acc_id', 'gene', 'protein', 'orig_ptm_residue'), sep = '\\|', remove = FALSE)
  } else {
    warning(glue::glue('Could not detect fasta header format. Returning without separating out header info. Header example: {long_data$fasta_id[1]}'))
  }


  # Option: return raw-ish long-format data
  if (return_long) return(long_data)

  # Check whether all kinases match known netphorest kinase families
  unknown_kinases <- dplyr::filter(long_data, !kinase_fam %in% netphorest_kinase_names)
  if (nrow(unknown_kinases) > 0) stop(glue('Unknown kinases found: {paste(unknown_kinases$kinase_fam, collapse = ", ")}'))

  # Reshape data into wide format
  print('Reshaping data into wide matrix-like format, this might take a while.')
  wide_data <- long_data %>%
    tidyr::pivot_wider(
      id_cols = c(dplyr::all_of(c('fasta_id', 'position', 'residue', 'fragment_11')), dplyr::any_of(c('acc_id', 'gene', 'protein', 'orig_ptm_residue'))),
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
    wide_data <- bind_cols(wide_data, empty_cols)
  }

  # Sort kinase columns
  wide_data <- wide_data %>%
    dplyr::relocate(dplyr::all_of(sort(netphorest_kinase_names)), .after = dplyr::last_col())

  return(wide_data)
}
