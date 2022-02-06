phosphosite_path <- system.file('extdata', 'phosphorylation_site_dataset_head', package = 'phosphocie')
kinsub_path <- system.file('extdata', 'kinase_substrate_dataset_head', package = 'phosphocie')
# phosphosite <- read_phosphosite(phosphosite_path)

test_that("read_phosphosite reads successfully", {
  phosphosite_read <- read_phosphosite(phosphosite_path)
  phosphosite_human_read <- phosphosite_read %>% dplyr::filter(organism == 'human')
  testthat::expect_equal(dim(phosphosite_human_read), c(25L, 14L))
})

test_that("read_kinsub reads successfully", {
  kinsub_read <- read_kinsub(kinsub_path)
  kinsub_human_read <- kinsub_read %>% dplyr::filter(organism == 'human' & kin_organism == 'human')
  testthat::expect_equal(dim(kinsub_human_read), c(30L, 18L))
})
