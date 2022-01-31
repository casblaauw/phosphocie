#' Read in PhosphoSitePlus phosphorylation site data
#'
#' Wrapper around [readr::read_tsv()] that cleans up the default phosphorylation
#' site data to have nicer columns, cleaner values, and a prepared unique site ID.
#'
#' @param path The path to the unzipped kinase substrate file.
#'
#' @keywords internal
#'
read_phosphosite <- function(path) {
  phosphosite_colnames <-
    c('gene',
      'protein',
      'acc_id',
      'human_chr_loc',
      'ptm_info',
      'site_grp_id',
      'organism',
      'mw_kd',
      'domain',
      'fragment_15',
      'lt_lit',
      'ms_lit',
      'ms_cst',
      'cst_cat'
    )
    # Read in data
  readr::read_tsv(path, skip = 4, col_names = phosphosite_colnames) %>%
    # Drop useless columns
    dplyr::select(-c('lt_lit', 'ms_lit', 'ms_cst', 'cst_cat')) %>%
    # Separate ptm column into residue, position, and type info
    tidyr::separate(ptm_info,
                    into = c('ptm_residue', 'ptm_type'),
                    sep = '-',
                    remove = TRUE) %>%
    tidyr::separate(ptm_residue,
                    into = c('residue', 'position'),
                    sep = 1,
                    remove = FALSE,
                    convert = TRUE
    ) %>%
    # Make protein names output-safe and generate unique site ID for FASTA
    dplyr::mutate(protein = stringr::str_replace_all(protein, ' ', '_'),
                  unique_id = paste(acc_id, gene, protein, ptm_residue, sep = '|')
    ) %>%
    # Reorder columns into consistent order
    dplyr::relocate(
      dplyr::all_of(c(
        'unique_id',
        'gene',
        'protein',
        'acc_id',
        'residue',
        'position',
        'fragment_15',
        'ptm_type',
        'organism',
        'human_chr_loc',
        'ptm_residue',
        'mw_kd',
        'domain'
      )
    ))
}

#' Read in PhosphoSitePlus kinase substrate data
#'
#' Wrapper around [readr::read_tsv()] that cleans up the default kinase-substrate
#' data to have nicer columns, cleaner values, and a prepared unique site ID.
#'
#' @param path The path to the unzipped kinase substrate file.
#'
#' @keywords internal
#'
read_kinsub <- function(path) {
  kinsub_colnames <-
    c('kin_gene',
      'kin_prot',
      'kin_acc_id',
      'kin_organism',
      'substrate',
      'gene_id',
      'acc_id',
      'gene',
      'organism',
      'ptm_residue',
      'site_grp_id',
      'fragment_15',
      'domain',
      'in_vivo',
      'in_vitro',
      'cst_cat'
    )

    # Read in data
  readr::read_tsv(path, skip = 4, col_names = kinsub_colnames) %>%
    # Drop useless columns
    dplyr::select(-c('cst_cat')) %>%
    # Separate ptm column into residue and position columns
    tidyr::separate(ptm_residue,
                    into = c('residue', 'position'),
                    sep = 1,
                    remove = FALSE,
                    convert = TRUE
    ) %>%
    # Turn 'X'/NA columns into boolean
    dplyr::mutate(
      in_vivo = if_else(in_vivo == 'X', TRUE, FALSE) %>% tidyr::replace_na(FALSE),
      in_vitro = if_else(in_vitro == 'X', TRUE, FALSE) %>% tidyr::replace_na(FALSE)
    ) %>%
    # Make protein names output-safe and generate unique site ID for FASTA
    dplyr::mutate(
      substrate = stringr::str_replace_all(substrate, ' ', '_'),
      unique_id = paste(acc_id, gene, substrate, ptm_residue, sep = '|')
    ) %>%
    dplyr::relocate(
      dplyr::all_of(c(
        'unique_id',
        'gene',
        'substrate',
        'acc_id',
        'residue',
        'position',
        'fragment_15',
        'kin_gene',
        'kin_prot',
        'kin_acc_id',
        'ptm_residue',
        'organism',
        'kin_organism',
        'site_grp_id',
        'gene_id',
        'domain',
        'in_vivo',
        'in_vitro'
      )
    ))
}
