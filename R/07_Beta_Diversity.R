# SETUP ####

## packages ####
library(tidyverse)
library(phyloseq)
library(microbiome)
library(vegan)

## functions ####
source("./R/functions.R")
set.seed(666)

## load data ####
full <- readRDS("./data/physeq_objects/endo_ecto_microbiomes_ps.RDS")
endo <- readRDS("./data/physeq_objects/endo_microbiome_ps.RDS")
ecto <- readRDS("./data/physeq_objects/ecto_microbiome_ps.RDS")

# ORDINATE ####

# ordinate and plot several ways
for(i in c("DCA", "CCA", "RDA", "NMDS", "MDS", "PCoA")){
  x <- endo %>% microbiome::transform('pa')
  ord <- ordinate(x,method = i,distance = distance(x, method = "jaccard", binary = TRUE))
  p <- plot_ordination(endo,ord,color="trait_orn_height_mean") + ggtitle(i)
  print(p)
}

# PCoA is best
ord <- ordinate(ecto,method = "PCoA",distance = distance(x, method = "jaccard", binary = TRUE))
ord2 <- ordinate(endo,method = "PCoA",distance = distance(x, method = "jaccard", binary = TRUE))

# fit trait data to PCoA
env <- envfit(ord$values, 
              as(ecto@sam_data,'data.frame') %>% select(starts_with("trait_")),
              permutations = 999, 
              na.rm = FALSE)
# export to file
sink("./output/envfit_results_table.txt")
print("Ectosporic")
env
sink(NULL)

env2_data <- as(endo@sam_data,'data.frame') %>% 
  dplyr::filter(sample_id %in% row.names(ord2$vectors)[1:nrow(ord2$values)]) %>% 
  select(starts_with("trait_"))
env2 <- envfit(ord2$values, 
               env2_data,
               permutations = 999, 
               na.rm = FALSE)
# export to file
sink("./output/envfit_results_table.txt",append = TRUE)
print("Endosporic")
env2
sink(NULL)


# PERMANOVA ####
mat <- ecto %>% microbiome::transform('compositional') %>% otu_table() %>% as('matrix')
perm_ecto <- adonis2(formula = mat ~ trait_orn_height_mean, data = ecto@sam_data %>% as("data.frame"))
perm_ecto
mat <- endo %>% microbiome::transform('compositional') %>% otu_table() %>% as('matrix')
perm_endo <- adonis2(formula = mat ~ trait_orn_height_mean, data = endo@sam_data %>% as("data.frame"))
perm_endo



