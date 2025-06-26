# SETUP ####

## packages ####
library(tidyverse)
library(dada2)
library(ShortRead)
library(Biostrings)

## functions ####
source("./R/functions.R")

# get taxonomy database files
genus_db <- "./taxonomy/silva_nr99_v138.2_toGenus_trainset.fa.gz"
species_db <- "./taxonomy/silva_nr99_v138.2_toSpecies_trainset.fa.gz"
assign_species_db <- "./taxonomy/silva_v138.2_assignSpecies.fa.gz"

# get possible asv table files
asv_tables <- list.files("./data/ASV_Tables",full.names = TRUE,pattern = "_ASV_Table.RDS")
asv <- readRDS(asv_tables)

# remove ASVs that don't have at least 100 observations
asv_thinned <- asv[,which(colSums(asv) >= 100)]


if(file.exists(genus_db)){
      outfile <- str_replace(asv_tables,"_ASV","_genus_Taxonomy")
      
      tax_genus <- assign_taxonomy_to_asv_table(asv.table=asv_thinned,
                                          tax.database=genus_db,
                                          multithread=parallel::detectCores(),
                                          random.seed=666,
                                          try.rc = TRUE,
                                          min.boot=70)
      # export as RDS
      saveRDS(tax_genus,outfile)
      
    } else {warning("Genus Taxonomic databases do not exist as specified.")}

if(file.exists(species_db)){
  outfile <- str_replace(asv_tables,"_ASV","_species_Taxonomy")
  
  tax_species <- assign_taxonomy_to_asv_table(asv.table=asv_thinned,
                                            tax.database=species_db,
                                            multithread=parallel::detectCores(),
                                            random.seed=666,
                                            try.rc = TRUE,
                                            min.boot=70)
  # export as RDS
  saveRDS(tax_species,outfile)
  
} else {warning("Species Taxonomic databases do not exist as specified.")}

if(file.exists(assign_species_db)){
  outfile <- str_replace(asv_tables,"_ASV","_assign_species_Taxonomy")
  
  # run assignSpecies algorithm and add to previous tax table
  tax_assign_species <-addSpecies(taxtab = tax_species,refFasta = assign_species_db,tryRC = TRUE)
   
  # export as RDS
  saveRDS(tax_assign_species,outfile)
  
} else {warning("assign_species Taxonomic databases do not exist as specified. Do not trust species-level assignments.")}


