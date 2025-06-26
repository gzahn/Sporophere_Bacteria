# SETUP ####

## packages ####
library(tidyverse)
library(phyloseq)
library(microbiome)
library(vegan)
library(corncob)
library(ranger)
library(vip)

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
  ord <- ordinate(x,method = i,distance = phyloseq::distance(x, method = "jaccard", binary = TRUE))
  p <- plot_ordination(endo,ord,color="trait_orn_height_mean") + ggtitle(i)
  print(p)
}

# PCoA is best
ord <- ordinate(ecto,method = "PCoA",distance = phyloseq::distance(x, method = "jaccard", binary = TRUE))
ord2 <- ordinate(endo,method = "PCoA",distance = phyloseq::distance(x, method = "jaccard", binary = TRUE))

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

ecto@sam_data$trait_investment_mean
# PERMANOVA ####
mat <- ecto %>% microbiome::transform('compositional') %>% otu_table() %>% as('matrix')
perm_ecto <- adonis2(formula = mat ~  trait_orn_height_mean + trait_shape_median + trait_vol_mean + trait_investment_mean,
                     data = ecto@sam_data %>% as("data.frame"),by = "margin")
perm_ecto
mat <- endo %>% microbiome::transform('compositional') %>% otu_table() %>% as('matrix')
perm_endo <- adonis2(formula = mat ~ trait_orn_height_mean + trait_shape_median + trait_vol_mean + trait_investment_mean,
                     data = endo@sam_data %>% as("data.frame"),by = "margin")
perm_endo

broom::tidy(perm_ecto) %>% 
  mutate(location="Ectosporic") %>% 
  full_join(
    broom::tidy(perm_endo) %>% 
      mutate(location="Endosporic")
  ) %>% 
  saveRDS("./output/PermANOVA_results_tables_endo-ecto_df.RDS")

# CORNCOB ####
# test for differential abundance by location (endo/ecto)
# watch for perfectly discriminatory taxa and remove them

ps_genus <- full %>% tax_glom("Genus",NArm = FALSE)


# how to account for AMF species or phylogeny here???
da_analysis <- 
  differentialTest(formula = ~ location, #abundance
                               phi.formula = ~ 1, #dispersion
                               formula_null = ~ 1, # control for __ abundance
                               phi.formula_null = ~ 1, # control for __ dispersion
                               test = "Wald", 
                               boot = FALSE,
                               B = 100,
                               data = full,
                               fdr_cutoff = 0.05,
                               full_output = FALSE)
beepr::beep()
# saveRDS(da_analysis,"./output/da_analysis_location-g_spp.RDS")
da_analysis$significant_taxa

# no taxa had significant differential abundance between endo and ecto locations!

# RANDOM FOREST ####

# build data frame
dat <- 
  data.frame(ectosporic = full@sam_data$location) %>% 
  bind_cols(full %>% microbiome::transform('compositional') %>% otu_table() %>% as("matrix")) %>% 
  mutate(ectosporic = case_when(ectosporic == "Ectosporic" ~ TRUE, TRUE ~ FALSE))

# fit model
rf_mod <- 
  ranger::ranger(data=dat,
                 formula = ectosporic ~ .,importance = 'permutation')
# get most important variables (taxa)
vi(rf_mod)$Importance %>% plot
important_taxa <- vi(rf_mod)[vi(rf_mod)$Importance > 0.002,]
important_taxa$taxonomy <- otu_to_taxonomy(important_taxa$Variable,full)

# plot those between the locations

full %>% 
  microbiome::transform('compositional') %>% 
  subset_taxa(taxa_names(full) %in% important_taxa$Variable) %>% 
  merge_samples("location",fun = 'sum') %>% 
  plot_bar2(fill="Family") +
  scale_fill_viridis_d(option = 'turbo') +
  labs(y="Composition",x="",caption="These taxa were identified by RF as highly predictive of microbiome location.") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face='bold.italic',size=12),
        legend.title = element_text(face='bold',size=16),
        legend.text = element_text(face='bold',size=12),
        axis.title = element_text(face='bold',size=16),
        axis.text.x = element_text(face='bold',size=14,hjust = 0,vjust=.5,angle=90))
ggsave("./output/figs/important_taxa_relabund_by_location-family.png")
