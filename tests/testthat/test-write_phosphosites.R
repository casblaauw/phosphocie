kinsub_path <- system.file('extdata', 'kinase_substrate_dataset_head', package = 'phosphocie')
kinsub <- read_kinsub(kinsub_path)

test_that("build_fastas works with supplied names", {
  tmp1 <- tempfile()
  build_fastas(head(kinsub), tmp1, name_col = 'gene', seq_col = 'fragment_15')
  expect_true(file.exists(tmp1))
})

test_that("build_fastas works without supplied names", {
  tmp1 <- tempfile()
  build_fastas(head(kinsub), tmp1, seq_col = 'fragment_15', header_pattern = '{acc_id}|{gene}|{substrate}|{residue}{position}')
  expect_true(file.exists(tmp1))
})

test_that("build_fastas creates the same file with header_pattern", {
  tmp1 <- tempfile()
  tmp2 <- tempfile()

  build_fastas(head(kinsub), tmp1, name_col = 'unique_id', seq_col = 'fragment_15')
  build_fastas(head(kinsub), tmp2, seq_col = 'fragment_15', header_pattern = '{unique_id}|{toupper(stringr::str_sub(.data[[seq_col]], ceiling(nchar(.data[[seq_col]])/2)-3, ceiling(nchar(.data[[seq_col]])/2)+3))}')
  expect_equal(readr::read_file(tmp1), readr::read_file(tmp2))
})

test_that("build_fastas headers have the site at the center of the sequence fragment", {
  tmp1 <- tempfile()

  build_fastas(head(kinsub), tmp1, name_col = 'unique_id', seq_col = 'fragment_15')
  center_aas <- readr::read_lines(tmp1) %>%
    magrittr::extract(seq(1, length(.), 2)) %>%
    stringr::str_extract('[a-zA-Z_]{3,9}$') %>%
    stringr::str_sub(ceiling(nchar(.)/2), ceiling(nchar(.)/2))
  site_aas <- head(kinsub)$ptm_residue %>% stringr::str_sub(1,1)

  expect_equal(toupper(center_aas), toupper(site_aas))
})
