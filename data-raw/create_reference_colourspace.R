# Creates a phosphocie reference colourspace, by predicting and representing
# kinase sites from the entire PhosphoSitePlus phosphorylation site dataset

library(phosphocie)
library(uwot)
library(magrittr)


# Download phospho data after logging in:
# https://www.phosphosite.org/staticDownloads -> Phosphorylation_site_dataset.gz
untar('Phosphorylation_site_dataset.gz', exdir = './PhosphoSitePlus')

# Read in massive phospho data
phosphosite <- phosphocie:::read_phosphosite('PhosphoSitePlus/Phosphorylation_site_dataset')
phosphosite_human <- phosphocie:::filter_phosphosite(phosphosite)

# Build fastas of site windows
build_fastas(
  phosphosite_human,
  name_col = 'unique_id',
  seq_col = 'fragment_15',
  path = 'PhosphoSitePlus/phosphosite_human_peptides.fasta'
)

# Predict kinases

# Requires downloading and compiling netphorest from netphorest.info:
# download.file('http://netphorest.info/download/NetPhorest_human_2.1.zip', file.path('./PhosphoSitePlus/NetPhorest_human_2.1.zip'))
# unzip('NetPhorest_human_2.1.zip', exdir = '.')

# system('cat PhosphoSitePlus/phosphosite_human_peptides.fasta | ./netphorest > PhosphoSitePlus/phosphosite_human_netphorest')

# Read in predictions
netphorest_kinase <- read_netphorest('PhosphoSitePlus/phosphosite_human_netphorest')
netphorest_kinase <- filter_netphorest(netphorest_kinase,
                                       source_window_size = 15,
                                       keep_uncertain = FALSE)

# Create reference UMAP
ref_kinase <- netphorest_kinase %>%
  tibble::column_to_rownames('fasta_id') %>%
  dplyr::select(Abl_group:YSK_group) %>%
  as.matrix()

ref_umap <- uwot::umap(netphorest_kinase, n_components = 3, ret_model = TRUE, verbose = TRUE, spread = 5)

# Create reference UCIE
ref_transform <- ucie_transformations(ref_umap$embedding)
ref_ucie <- kinase2cielab(ref_umap$embedding, ref_transform, LAB_coordinates = TRUE) %>%
  tibble::column_to_rownames('rownames') %>%
  as.matrix()

# If internal sysdata doesn't exist yet, initialise with use_data
if (!file.exists(system.file('R', 'sysdata.rda', package = 'phosphocie'))) {
  usethis::use_data(ref_kinase, ref_ucie, internal = TRUE)
} else {
  # Else: append by loading and re-saving
  sysdata_filenames <- load("R/sysdata.rda")
  save(
    list = c(sysdata_filenames, "ref_kinase", "ref_ucie"),
    file = file.path(system.file('R', package = 'phosphocie'), 'sysdata.rda')
  )
}



