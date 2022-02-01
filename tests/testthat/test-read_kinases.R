kinsub_netphorest_path <- system.file('extdata', 'kinsub_human_netphorest', package = 'phosphocie')

test_that("read_netphorest can read in data as long format", {
  kinsub_netphorest_path <- system.file('extdata', 'kinsub_human_netphorest', package = 'phosphocie')
  kinsub_netphorest <- read_netphorest(kinsub_netphorest_path, return_long = TRUE)
  expect_equal(dim(kinsub_netphorest), c(3458L, 14L))
})

test_that("read_netphorest can read in data as wide format", {
  kinsub_netphorest_path <- system.file('extdata', 'kinsub_human_netphorest', package = 'phosphocie')
  kinsub_netphorest <- read_netphorest(kinsub_netphorest_path, return_long = FALSE)
  expect_equal(dim(kinsub_netphorest), c(87L, 68L))
})
