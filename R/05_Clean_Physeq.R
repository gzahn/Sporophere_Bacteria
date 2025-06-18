# SETUP ####

## packages ####
library(tidyverse)
library(phyloseq)

## functions ####
source("./R/functions.R")

## load data ####
ps <- readRDS("./data/physeq_objects/16S_V3V4_physeq_object_species.RDS")

## sanity checks
ps
ps@tax_table[,7] %>% unname %>% unique
sample_names(ps)
ps@tax_table[1,] %>% unname
taxa_sums(ps) %>% summary 
sample_sums(ps) %>% plot

## remove ambiguous taxa ####
# remove taxa with missing phylum
ps <- 
  ps %>% 
  subset_taxa(!is.na(Phylum))
# remove chloroplasts
ps <- subset_taxa(ps,Class != "Chloroplast")
# remove mitochondria
ps <- subset_taxa(ps,Family != "Mitochondria")


## remove bad samples ####
# remove samples with fewer than 100 ASV observations
ps <- 
  ps %>% 
  subset_samples(sample_sums(ps) > 100)

# remove negative control samples still remaining (N=4)
ps <- 
  ps %>% 
  subset_samples(treatment != "negative_ctl")

## data wrangling in metadata ####

# check G-spp names
as.character(na.omit(ps@sam_data$g_spp %>% unique))
# make them a factor (by ornamentation level, for now)
ps@sam_data$g_spp <- factor(ps@sam_data$g_spp,
                            levels = ps@sam_data$g_spp[order(ps@sam_data$trait_orn_height_mean)] %>% 
                              unique %>% 
                              na.omit %>% 
                              as.character())
# same for family
ps@sam_data$family <- factor(ps@sam_data$family,
                            levels = ps@sam_data$family[order(ps@sam_data$trait_orn_height_mean)] %>% 
                              unique %>% 
                              na.omit %>% 
                              as.character())
# same for genus
ps@sam_data$genus <- factor(ps@sam_data$genus,
                             levels = ps@sam_data$genus[order(ps@sam_data$trait_orn_height_mean)] %>% 
                               unique %>% 
                               na.omit %>% 
                               as.character())
# same for order
ps@sam_data$order <- factor(ps@sam_data$order,
                             levels = ps@sam_data$order[order(ps@sam_data$trait_orn_height_mean)] %>% 
                               unique %>% 
                               na.omit %>% 
                               as.character())
# same for isolate
ps@sam_data$isolate <- factor(ps@sam_data$isolate,
                             levels = ps@sam_data$isolate[order(ps@sam_data$trait_orn_height_mean)] %>% 
                               unique %>% 
                               na.omit %>% 
                               as.character())





## split physeq ####

# INVAM bulk soil
invam <- 
  ps %>% 
  subset_samples(treatment == "invam_bulk_soil")
invam <- 
  invam %>% 
  subset_taxa(taxa_sums(invam) > 0)

# NEON soil for checking later
neon <- 
  ps %>% 
  subset_samples(treatment == "NEON")
neon <- 
  neon %>% 
  subset_taxa(taxa_sums(neon) > 0)

# Main study samples
ps <- 
  ps %>% 
  subset_samples(treatment %in% c("unsterile","sterile"))
ps <- 
  ps %>% 
  subset_taxa(taxa_sums(ps) > 0)

# EXPORT ####
saveRDS(ps,"./data/physeq_objects/clean_physeq_object.RDS")
saveRDS(invam,"./data/physeq_objects/clean_invam_bulk_soil_physeq_object.RDS")
saveRDS(neon,"./data/physeq_objects/clean_neon_soil_physeq_object.RDS")
