library(stringr)
library(magrittr)

download.file('http://netphorest.info/download/NetPhorest_human_2.1.zip', file.path('./NetPhorest_human_2.1.zip'))
unzip('NetPhorest_human_2.1.zip', exdir = './NetPhorest')

pssm_kinases <- readLines('NetPhorest/pssm_code.h') %>%
  grep(pattern = '\\thuman\\tKIN', x = ., value = TRUE, fixed = TRUE) %>%
  stringr::str_match('KIN\\\\t([\\w\\d]*)\\\\t') %>%
  .[,2]

nn_kinases <- readLines('NetPhorest/nn_code.h') %>%
  grep(pattern = '\\thuman\\tKIN', x = ., value = TRUE, fixed = TRUE) %>%
  stringr::str_match('KIN\\\\t([\\w\\d]*)\\\\t') %>%
  .[,2]
netphorest_known_kinases <- c(pssm_kinases, nn_kinases)
# dput(netphorest_known_kinases)

usethis::use_data(netphorest_known_kinases, overwrite = TRUE, internal = TRUE)
