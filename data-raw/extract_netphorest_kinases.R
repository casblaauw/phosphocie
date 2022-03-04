library(stringr)
library(magrittr)

tmp_dir <- tempdir()
download.file('http://netphorest.info/download/NetPhorest_human_2.1.zip', file.path(tmp_dir, '/NetPhorest_human_2.1.zip'))
unzip(file.path(tmp_dir, '/NetPhorest_human_2.1.zip'), exdir = file.path(tmp_dir, 'NetPhorest'))

pssm_kinases <- readLines(file.path(tmp_dir, 'NetPhorest/pssm_code.h')) %>%
  grep(pattern = '\\thuman\\tKIN', x = ., value = TRUE, fixed = TRUE) %>%
  stringr::str_match('KIN\\\\t([\\w\\d]*)\\\\t') %>%
  .[,2]

nn_kinases <- readLines(file.path(tmp_dir, 'NetPhorest/nn_code.h')) %>%
  grep(pattern = '\\thuman\\tKIN', x = ., value = TRUE, fixed = TRUE) %>%
  stringr::str_match('KIN\\\\t([\\w\\d]*)\\\\t') %>%
  .[,2]
netphorest_known_kinases <- c(pssm_kinases, nn_kinases)
# dput(netphorest_known_kinases)


# If internal sysdata doesn't exist yet, initialise with use_data
if (!file.exists(system.file('R', 'sysdata.rda', package = 'phosphocie'))) {
  usethis::use_data(netphorest_known_kinases, internal = TRUE)
} else {
  # Else: append by loading and re-saving
  sysdata_filenames <- load("R/sysdata.rda")
  save(
    list = c(sysdata_filenames, "netphorest_known_kinases"),
    file = file.path(system.file('R', package = 'phosphocie'), 'sysdata.rda')
  )
}
