# SETUP ####

## packages ####
library(tidyverse)
library(dada2)
library(phyloseq)

## functions ####
source("./R/functions.R")

## paths ####
asv_table_dir <- "./data/ASV_Tables"

# find asv table files
asv_tables <- list.files(asv_table_dir,full.names = TRUE, pattern = "_ASV_Table.RDS")

# find taxonomy table files
tax_tables <- list.files(asv_table_dir,full.names = TRUE, pattern = "_Taxonomy_Table.RDS")

# out path
out_paths <- paste0(asv_tables %>% str_remove("_ASV_Table.RDS"),"_physeq_object_",
                    c("genus","species","assigned_species"),
                    ".RDS")
dir.create("./data/physeq_objects")
out_paths <- out_paths %>% str_replace("ASV_Tables","physeq_objects")
# DATA ####

# get metadata
meta <- read_csv("./data/metadata/sporosphere_metadata.csv")

# filter metadata to match remaining samples
# we removed 1,214 potential contaminant taxa during previous script
# we also removed a number of taxa that had low (<100) total observations
asv <- readRDS(asv_tables)
dim(asv)
lost_samples_metadata <- meta %>% dplyr::filter(sample_id %ni% row.names(asv))
saveRDS(lost_samples_metadata,"./output/lost_samples_metadata.RDS")
meta$sample_id %in% row.names(asv)
meta <- meta %>% dplyr::filter(sample_id %in% row.names(asv))

# BUILD PS OBJECTS ####
for(i in seq_along(tax_tables)){
  
  # load asv table
  asv <- readRDS(asv_tables)
  
  # load taxonomy table
  taxa <- readRDS(tax_tables[i])
  
  # set up physeq components
  otu <- otu_table(asv,taxa_are_rows = FALSE)
  met <- sample_data(meta)
  sample_names(met) <- meta$sample_id
  sample_names(otu) <- meta$sample_id[1:nsamples(otu)]
  tax <- tax_table(taxa)
  
  # build physeq object
  physeq <- phyloseq(otu,
                     met,
                     tax
  )
  
  # save physeq object
  saveRDS(physeq,out_paths[i])
  
}


