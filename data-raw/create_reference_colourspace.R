# Creates a phosphocie reference colourspace, by predicting and representing
# kinase sites from the entire PhosphoSitePlus phosphorylation site dataset

library(phosphocie)
library(umap)
library(magrittr)


# Download phospho data
download.file('https://www.phosphosite.org/downloads/Phosphorylation_site_dataset.gz', file.path('./Phosphorylation_site_dataset.gz'))
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

# system('cat PhosphoSitePlus/phosphosite_human_peptides.fasta | .PhosphoSitePlus/netphorest > PhosphoSitePlus/phosphosite_human_netphorest')

# Read in predictions
netphorest_kinase <- read_netphorest('PhosphoSitePlus/phosphosite_human_netphorest')
netphorest_kinase <- filter_netphorest(netphorest_kinase, source_window_size = 15, keep_uncertain = FALSE)

# Create reference UMAP
umap_settings <- umap.defaults
umap_settings$n_components <- 3
netphorest_kinase_f <- netphorest_kinase_f %>%
  tibble::column_to_rownames('fasta_id') %>%
  dplyr::select(-c(dplyr::all_of(c('position', 'residue', 'fragment_11')), dplyr::any_of('orig_res', 'orig_pos')))
ref_umap <- umap(netphorest_kinase_f, umap_settings)

# Create reference UCIE
ref_transform <- ucie_transformations(ref_umap$layout)
ref_ucie <- kinase2cielab(netphorest_kinase, ref_transform, LAB_coordinates = TRUE)


usethis::use_data(ref_umap, overwrite = TRUE)
usethis::use_data(ref_transform, overwrite = TRUE)


