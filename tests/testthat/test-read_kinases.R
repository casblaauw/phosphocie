# Test read_netphorest
kinsub_netphorest_path <- system.file('extdata', 'kinsub_human_netphorest', package = 'phosphocie')

test_that("read_netphorest can read in data as wide format", {
  kinsub_netphorest <- read_netphorest(kinsub_netphorest_path)
  expect_equal(dim(kinsub_netphorest), c(87L, 64L))
})

test_that("read_netphorest can read in data as long format", {
  kinsub_netphorest <- read_netphorest(kinsub_netphorest_path, return_long = TRUE)
  expect_equal(dim(kinsub_netphorest), c(3458L, 10L))
})

test_that("read_netphorest can split fasta header", {
  kinsub_netphorest <- read_netphorest(kinsub_netphorest_path, split_fasta_header = TRUE)
  expect_equal(dim(kinsub_netphorest), c(87L, 64L))
})


# Test filter_netphorest
ambiguous_data <- data.frame( # True value is 5, but ambiguous; w/o info, picks 7
    fasta_id = c("P13796|LCP1|L-plastin|S5", "P13796|LCP1|L-plastin|S5"),
    position = c(5, 7),
    residue = c("S", "T"),
    fragment_11 = c("-MARGsVSDEE", "ARGSVtDEEMM"),
    protein_res_col = c("S", "S"),
    orig_pos_col = c(5, 5),
    Abl_group = c(0, 0)
  )


kinsub_netphorest <- kinsub_netphorest_path %>%
  read_netphorest()

kinsub_netphorest_small <- kinsub_netphorest %>%
  dplyr::filter(fasta_id %in% head(unique(.$fasta_id), 10))


test_that("filter_netphorest can filter read_netphorest data by default", {
  expect_equal(
    nrow(filter_netphorest(kinsub_netphorest_small, source_window_size = 15)),
    10
  )
})

test_that("filter_netphorest output has only unique IDs", {
  expect_equal(
    anyDuplicated(filter_netphorest(kinsub_netphorest, source_window_size = 15)$fasta_id),
    0
  )
})

test_that("filter_netphorest output has only unique IDs if keep_uncertain = FALSE", {
  expect_equal(
    anyDuplicated(filter_netphorest(kinsub_netphorest, source_window_size = 15)$fasta_id),
    0
  )
})

test_that("filter_netphorest returns warnings on ambiguous data", {
  expect_message(
    filter_netphorest(ambiguous_data, source_window_size = 15),
    ".*P13796\\|LCP1\\|L-plastin\\|S5.*"
    )
})

test_that("filter_netphorest filters ambiguous data by default", {
  expect_equal(
    nrow(filter_netphorest(ambiguous_data, source_window_size = 15)),
    1
  )
})

test_that("filter_netphorest keeps ambiguous data if keep_uncertain = NULL", {
  expect_equal(
    nrow(filter_netphorest(ambiguous_data, source_window_size = 15, keep_uncertain = NULL)),
    1
  )
})

test_that("filter_netphorest keeps ambiguous data if keep_uncertain = TRUE", {
  expect_equal(
    nrow(filter_netphorest(ambiguous_data, source_window_size = 15, keep_uncertain = TRUE)),
    2
  )
})

test_that("filter_netphorest keeps ambiguous data if keep_uncertain = TRUE", {
  expect_equal(
    nrow(filter_netphorest(ambiguous_data, source_window_size = 15, keep_uncertain = FALSE)),
    0
  )
})

test_that("filter_netphorest can fix ambiguous data with protein_res_col", {
  expect_equal(
    filter_netphorest(ambiguous_data,
                      source_window_size = 15,
                      protein_res_col = 'protein_res_col'
                      )$position,
    5
  )
})

test_that("filter_netphorest can fix ambiguous data with protein_res_col and detected_res_col", {
  expect_equal(
    filter_netphorest(ambiguous_data,
                      source_window_size = 15,
                      protein_res_col = 'protein_res_col',
                      detected_res_col = 'residue'
                      )$position,
    5
  )
})

test_that("filter_netphorest can fix ambiguous data with protein_pos_col", {
  expect_equal(
    filter_netphorest(ambiguous_data,
                      source_window_size = 15,
                      protein_pos_col = 'orig_pos_col'
                      )$position,
    5
  )
})
