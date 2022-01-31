kinsub_path <- system.file('extdata', 'Kinase_Substrate_Dataset_head', package = 'phosphocie')
kinsub <- read_kinsub(kinsub_path)

test_that("build_fastas works with supplied names", {
  tmp1 <- tempfile()
  build_fastas(head(kinsub), tmp1, name_col = 'gene', seq_col = 'fragment_15')
  expect_true(file.exists(tmp1))
})

test_that("build_fastas creates the same file with name_pattern", {
  tmp1 <- tempfile()
  tmp2 <- tempfile()

  build_fastas(head(kinsub), tmp1, name_col = 'unique_id', seq_col = 'fragment_15')
  build_fastas(head(kinsub), tmp2, seq_col = 'fragment_15', name_pattern = '{acc_id}|{gene}|{substrate}|{residue}{position}')
  expect_equal(readr::read_file(tmp1), readr::read_file(tmp2))
})
