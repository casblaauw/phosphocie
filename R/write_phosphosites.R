#' Create an input file for KSP
#'
#' Creates an input file to use for prediction with [KSP](https://doi.org/10.1186/s12864-020-06895-2).
#' Wrapper around [readr::write_tsv()].
#'
#' @param data Data frame in long format, containing gene names and sequence windows.
#' @param path Path to write the file to, including file name and extension.
#' @param name_col Name of the column that contains gene names.
#' @param seq_col Name of the column that contains sequence windows.
#'
#' @return Returns path that file was written to, invisibly.
#' @export
#'
#' @examples
#' kinsub_path <- system.file('extdata', 'Kinase_Substrate_Dataset_head', package = 'phosphocie')
#' kinsub <- phosphocie:::read_kinsub(kinsub_path)
#' tmp <- tempfile()
#'
#' write_ksp_input(kinsub, tmp, name_col = 'gene', seq_col = 'fragment_15')

write_ksp_input <- function(data, path, name_col, seq_col) {
  if (any(!c(name_col, seq_col) %in% colnames(data))) {
    rlang::abort(glue::glue("Please ensure your data contains the supplied column names. \nSupplied names: {paste(name_col, seq_col, sep = ', ')}\nDetected colnames: {paste(colnames(data), collapse = ', ')}"))
  }

  data %>%
    dplyr::select(dplyr::all_of(c(name_col, seq_col))) %>%
    readr::write_tsv(path, col_names = FALSE)
  print(glue::glue('file created at {path}'))
  return(invisible(path))

}

#' Create an input file for KSP's PWM mode
#'
#' Creates an input file to use for prediction with
#' [KSP](https://doi.org/10.1186/s12864-020-06895-2)'s position-weight matrix-based scoring tool.
#' Wrapper around [readr::write_tsv()].
#'
#' @param data Data frame in long format, containing gene names and sequence windows.
#' @param path Path to write the file to, including file name and extension.
#' @param name_col Name of the column that contains gene names.
#' @param seq_col Name of the column that contains sequence windows.
#'
#' @return Returns path that file was written to, invisibly.
#' @export
#'
#' @examples
#' kinsub_path <- system.file('extdata', 'Kinase_Substrate_Dataset_head', package = 'phosphocie')
#' kinsub <- phosphocie:::read_kinsub(kinsub_path)
#' tmp <- tempfile()
#'
#' write_pwm_input(kinsub, tmp, name_col = 'gene', seq_col = 'fragment_15')

write_pwm_input <- function(data, path, name_col, seq_col) {
  if (any(!c(name_col, seq_col) %in% colnames(data))) {
    rlang::abort(glue::glue("Please ensure your data contains the supplied column names. \nSupplied names: {paste(name_col, seq_col, sep = ', ')}\nDetected colnames: {paste(colnames(data), collapse = ', ')}"))
  }

  data %>%
    dplyr::select(dplyr::all_of(c(name_col, seq_col))) %>%
    dplyr::bind_cols(placeholder = ' ', .) %>%
    readr::write_tsv(path, col_names = FALSE)
  print(glue::glue('file created at {path}'))
  return(invisible(path))
}

#' Retrieve full protein sequences from UniProt
#'
#' From a vector of uniprot accession IDs, retrieves FASTAs for each protein from
#' the [EMBL Proteins API](https://www.ebi.ac.uk/proteins/api/doc/index.html) and
#' concatenates these into one massive FASTA file.
#' Can serve as NetPhorest input, but [build_fastas()] is the recommended alternative
#' to predict for known sites rather than entire proteins.
#'
#' @param uniprot_acc Character vector of UniProt accession IDs. Supports isoform IDs like `Q02297-6`.
#' @param path Path to write the file to, including file name and extension.
#'
#' @details
#' * Filters out duplicated IDs and invalid IDs.
#' * Requests in batches of 100, which is the max for the API.
#' * If any batch fails, the process is stopped and the file is left as is, containing all
#'   successful batches up to then.
#' * Between each batch request, waits for 0.75 seconds to stay below rate limits.
#'
#' @return Returns path of the output FASTA, invisibly.
#' @export
#'
#' @examples
#' kinsub_path <- system.file('extdata', 'Kinase_Substrate_Dataset_head', package = 'phosphocie')
#' kinsub <- phosphocie:::read_kinsub(kinsub_path)
#' tmp <- tempfile()
#'
#' \dontrun{
#'   retrieve_fastas(kinsub$acc_id, tmp)
#' }

retrieve_fastas <- function(uniprot_acc, path) {

  # Filter uniprot IDs
  uniprot_acc <- unique(uniprot_acc)

  uniprot_regex <- '[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}'
  if (any(!stringr::str_detect(uniprot_acc, uniprot_regex))) {
    warning(glue::glue("Accession IDs with invalid uniprot ID formats detected. These will be ignored. {paste(uniprot_acc[!str_detect(uniprot_acc, uniprot_regex)], collapse = ',')}"))
    uniprot_acc <- uniprot_acc[stringr::str_detect(uniprot_acc, uniprot_regex)]
  }

  if (any(nchar(uniprot_acc) > 10)) {
    warning(glue::glue("Accession IDs with more than 10 characters detected, which is not allowed. These will be ignored. {paste(uniprot_acc[nchar(uniprot_acc) > 10])}" ))
    uniprot_acc <- uniprot_acc[nchar(uniprot_acc) <= 10]
  }

  # Loop over and retrieve for every batch of 100 IDs
  sections <- split(uniprot_acc, ceiling(seq_along(uniprot_acc)/100))
  counter <- 0
  open_file <- file(path, open = 'w+')
  for (section in sections) {
    counter <- counter + 1
    print(glue::glue('Requesting section {counter}/{length(sections)}'))
    requestURL <- glue::glue("https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=-1&accession={paste(section, collapse = ',')}")
    r <- httr::GET(requestURL, httr::accept("text/x-fasta"))

    if (httr::http_error(r)) {
      close(open_file)
      print(glue::glue("Problem with section {counter}. Accession IDs of current section: \n {paste(section, collapse = ',')}"))
      httr::stop_for_status(r)
    }
    writeLines(httr::content(r, as = 'text', encoding = "UTF-8"), open_file, sep = '')
    Sys.sleep(0.75)
  }
  close(open_file)
  print(glue::glue('Successfully retrieved fasta, saved at {path}'))
  return(invisible(path))
}

#' Create a FASTA file with sequence fragments
#'
#' From (potential) sites and their surrounding amino acids, create a FASTA-conforming file.
#' Requires a column of unique names for each site.
#'
#'
#' @param data Data frame in long format, containing unique names and sequence windows.
#' @param path Path to write the file to, including file name and extension.
#' @param name_col Name of the column that contains unique fasta IDs.
#'   If missing, uses name_pattern to construct unique IDs.
#' @param seq_col Name of the column that contains the sequence windows.
#' @param name_pattern Optional, only used if `name_col` is missing.
#'   Pattern for use in [glue::glue()], used to construct a new unique name column
#'   based on current columns.
#'
#' @return Returns path of the output FASTA, invisibly.
#' @export
#'
#' @examples
#' kinsub_path <- system.file('extdata', 'Kinase_Substrate_Dataset_head', package = 'phosphocie')
#' kinsub <- phosphocie:::read_kinsub(kinsub_path)
#' tmp <- tempfile()
#'
#' build_fastas(kinsub, tmp, name_col = 'unique_id', seq_col = 'fragment_15')
#'
#' build_fastas(kinsub, tmp, seq_col = 'fragment_15', name_pattern = '{acc_id}|{gene}|{substrate}|{residue}{position}')

build_fastas <- function(data, path, name_col, seq_col, name_pattern = NULL) {

  # Create name column if missing
  if (missing(name_col)) {
    if (is.null(name_pattern)) {
      rlang::abort('Either name_col or name_pattern must be supplied.')
    }
    data <- dplyr::mutate(data, build_fasta_id = glue::glue(name_pattern))
    name_col <- 'build_fasta_id'
    if (anyDuplicated(data[[name_col]]) > 0) {
      rlang::abort('Name column created with name_pattern is not unique, please supply a different pattern.')
    }
  }

  # Check column names and data
  if (any(!c(name_col, seq_col) %in% colnames(data))) {
    rlang::abort(glue::glue("Please ensure your data contains the supplied column names. \nSupplied names: {paste(name_col, seq_col, sep = ', ')}\nDetected colnames: {paste(colnames(data), collapse = ', ')}"))
  }

  if (any(nchar(data[[seq_col]]) > 70)) {
    rlang::abort("Peptide of 70+ AA detected. Please check your data.")
  }

  # Remove any filling underscores and other junk by extracting first alphanumeric group
  data[[seq_col]] <- stringr::str_extract(data[[seq_col]], '[:alpha:]+')

  # Walk over name and sequence at once and write to file
  open_file <- file(path, open = 'w+')
  purrr::pwalk(
    dplyr::select(data, dplyr::all_of(c(name_col, seq_col))),
    ~writeLines(
      text = c(
        paste0('>', ..1), # Write start line
        toupper(..2)     # Write sequence line, cut off any accidental
      ),
      con = open_file)
  )
  close(open_file)

  print(glue::glue('Successfully built fasta, saved at {path}'))
  return(invisible(path))
}
