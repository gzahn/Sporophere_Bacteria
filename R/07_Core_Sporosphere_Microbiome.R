# SETUP ####

## packages ####
library(tidyverse)
library(phyloseq)
library(microbiome)
library(vegan)

## functions ####
source("./R/functions.R")

## load data ####
ps <- readRDS("./data/physeq_objects/clean_physeq_object.RDS")
invam <- readRDS("./data/physeq_objects/clean_invam_bulk_soil_physeq_object.RDS")


# FIND CORE MICROBIOMES ####
# remove taxa found in  bulk soil
# check overlap
taxa_names(ps) %in% taxa_names(invam) %>% summary


## Get taxa lists, unique to isolate and physical location ####

# need to find (in each isolate) which taxa are present in the sterile spores
# and remove them from the unsterile spores

# make list of physeq objects (one for each isolate)
# unsterile = taxa found on surface of spore, not found internally
# sterile = taxa found internally, not found in bulk soil or spore surface
isolate_unsterile_taxa_ps <- list()
isolate_sterile_taxa_ps <- list()

for(iso in levels(ps@sam_data$isolate)){
  
  x <- 
    ps %>% 
    subset_samples(isolate == iso)
  
  # if no invam bulk soil present for a given isolate, skip this step
  if(any(invam@sam_data$isolate == iso)){
    y <- 
      invam %>% 
      subset_samples(isolate == iso)
    y <- 
      y %>% 
      subset_taxa(taxa_sums(y) > 1)
    x.sterile <- 
      x %>% 
      subset_samples(treatment == "sterile")
    x.sterile <- 
      x.sterile %>% 
      subset_taxa(taxa_sums(x.sterile) > 0 & taxa_names(x.sterile) %ni% taxa_names(y)) # remove soil taxa from internal
    x.sterile.taxa <- taxa_names(x.sterile)
    
    x.unsterile <- 
      x %>% 
      subset_samples(treatment == "unsterile")
    x.unsterile <- 
      x.unsterile %>% 
      subset_taxa(taxa_sums(x.unsterile) > 0)
    x.unsterile <- 
      x.unsterile %>% 
      subset_taxa(taxa_names(x.unsterile) %ni% x.sterile.taxa)
    
    # now, remove any surface taxa from the sterile treatment, to be conservative and find
    # taxa that ONLY exist inside the spore
    # (this will remove internal taxa that are ALSO present on the surface)
    
    # THIS IS WRONG!!! comment it out!
    # x.sterile <- 
    #   x.sterile %>% 
    #   subset_taxa(taxa_names(x.sterile) %ni% taxa_names(x.unsterile))
    
    # save individual physeq for each group
    isolate_unsterile_taxa_ps[[iso]] <- x.unsterile
    isolate_sterile_taxa_ps[[iso]] <- x.sterile
    
  } else {
    x.sterile <- 
      x %>% 
      subset_samples(treatment == "sterile")
    x.sterile <- 
      x.sterile %>% 
      subset_taxa(taxa_sums(x.sterile) > 0) # remove soil taxa from internal
    x.sterile.taxa <- taxa_names(x.sterile)
    
    x.unsterile <- 
      x %>% 
      subset_samples(treatment == "unsterile")
    x.unsterile <- 
      x.unsterile %>% 
      subset_taxa(taxa_sums(x.unsterile) > 0)
    x.unsterile <- 
      x.unsterile %>% 
      subset_taxa(taxa_names(x.unsterile) %ni% x.sterile.taxa)
    
    # now, remove any surface taxa from the sterile treatment, to be conservative and find
    # taxa that ONLY exist inside the spore
    # (this will remove internal taxa that are ALSO present on the surface)
    x.sterile <- 
      x.sterile %>% 
      subset_taxa(taxa_names(x.sterile) %ni% taxa_names(x.unsterile))
    
    # save individual physeq for each group
    isolate_unsterile_taxa_ps[[iso]] <- x.unsterile
    isolate_sterile_taxa_ps[[iso]] <- x.sterile
  }
}

## Define core microbiome parameters ####

# Get a feel for appropriate "core" definition
x <- isolate_unsterile_taxa_ps$`MG101B-12` %>% microbiome::transform("compositional")
y <- isolate_unsterile_taxa_ps$`MG101B-12` %>% microbiome::transform("compositional") %>% 
  microbiome::core(detection = .5,prevalence = .05)

det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)
plot_core(x, 
          prevalences = prevalences, 
          detections = det, 
          plot.type = "lineplot") + 
  xlab("Relative Abundance (%)")

# core heatmap
prevalences <- seq(.01, 1, .01)
detections <- round(10^seq(log10(0.01), log10(.2), length = 9), 3)
gray <- gray(seq(0,1,length=5))
plot_core(x,
               plot.type = "heatmap", 
               colours = gray,
               prevalences = prevalences, 
               detections = detections) +
  labs(x = "Detection Threshold\n(Relative Abundance (%))") +
  theme(axis.text.y= element_blank(),
        axis.text.x.bottom=element_text(size=8),
        axis.title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10))

## Get core microbiomes ####

# get list of core taxa for each isolate (detection=.2 'relabund',prevalence='more than half of samples')
# good conservative definition, based on plots
get_core_members <- function(x){
  x <- x %>% microbiome::transform("compositional")
  microbiome::core_members(x,detection = 0.2,prevalence = (nsamples(x)-1)/nsamples(x),include.lowest = TRUE)
}
isolate_core_members <- map(isolate_unsterile_taxa_ps,get_core_members)
isolate_core_ntaxa <- isolate_core_members %>% map_dbl(length) %>% unlist
isolate_core_members <- 
  data.frame(isolate = rep(names(isolate_core_ntaxa),isolate_core_ntaxa),
           asv = isolate_core_members %>% purrr::reduce(c)) %>% 
  mutate(taxonomy = corncob::otu_to_taxonomy(asv,ps)) %>% 
  left_join(as(ps@sam_data,'data.frame') %>% dplyr::select(isolate,starts_with("trait"))) %>% 
  unique.data.frame()
saveRDS(isolate_core_members,"./output/surface_bacteria_core_membership_by_isolate_df.RDS")

isolate_core_members2 <- map(isolate_sterile_taxa_ps,get_core_members)
isolate_core_ntaxa <- isolate_core_members2 %>% map_dbl(length) %>% unlist
isolate_core_members2 <- 
  data.frame(isolate = rep(names(isolate_core_ntaxa),isolate_core_ntaxa),
             asv = isolate_core_members2 %>% purrr::reduce(c)) %>% 
  mutate(taxonomy = corncob::otu_to_taxonomy(asv,ps)) %>% 
  left_join(as(ps@sam_data,'data.frame') %>% dplyr::select(isolate,starts_with("trait"))) %>% 
  unique.data.frame()
saveRDS(isolate_core_members2,"./output/internal_bacteria_core_membership_by_isolate_df.RDS")





## Regress core taxa numbers against traits ####

# N core membership vs traits (glm)
n_core_vs_traits_mod <- 
isolate_core_members %>% 
  group_by(isolate) %>% 
  summarize(N_core = n(),
            trait_orn_height_mean = mean(trait_orn_height_mean),
            trait_color_most = mean(trait_color_most),
            trait_vol_mean = mean(trait_vol_mean),
            trait_shape_median = mean(trait_shape_median),
            trait_investment_mean = mean(trait_investment_mean)) %>% 
  select(-isolate) %>% 
  glm(data = .,
      formula = N_core ~ .)
summary(n_core_vs_traits_mod)
saveRDS(n_core_vs_traits_mod,"./output/n_ecto_coretaxa_vs_traits_mod.RDS")

n_core_vs_traits_mod2 <- 
  isolate_core_members2 %>% 
  group_by(isolate) %>% 
  summarize(N_core = n(),
            trait_orn_height_mean = mean(trait_orn_height_mean),
            trait_color_most = mean(trait_color_most),
            trait_vol_mean = mean(trait_vol_mean),
            trait_shape_median = mean(trait_shape_median),
            trait_investment_mean = mean(trait_investment_mean)) %>% 
  select(-isolate) %>% 
  glm(data = .,
      formula = N_core ~ .)
summary(n_core_vs_traits_mod2)
saveRDS(n_core_vs_traits_mod2,"./output/n_endo_coretaxa_vs_traits_mod.RDS")


## PLOT DIVERSITY BY LOCATION ####

# bind them together, having removed bacteria found in the surface-sterile controls
# from each isolate grouping
ecto_ps <- 
  isolate_unsterile_taxa_ps %>% 
  purrr::reduce(merge_phyloseq)
ecto_ps@sam_data$location <- "Ectosporic"
ecto_ps <- 
  ecto_ps %>% 
  subset_samples(sample_sums(ecto_ps)>0)
ecto_ps <- 
  ecto_ps %>% 
  subset_taxa(taxa_sums(ecto_ps) >= 100)
taxa_sums(ecto_ps) %>% summary

# bind together endosporic bacteria
endo_ps <- 
  isolate_sterile_taxa_ps %>% 
  purrr::reduce(merge_phyloseq)
endo_ps@sam_data$location <- "Endosporic"
endo_ps <- 
  endo_ps %>% 
  subset_samples(sample_sums(endo_ps)>0)
endo_ps <- 
  endo_ps %>% 
  subset_taxa(taxa_sums(endo_ps) >= 100)
taxa_sums(endo_ps) %>% summary

# export individual physeqs
saveRDS(endo_ps,"./data/physeq_objects/endo_microbiome_ps.RDS")
saveRDS(ecto_ps,"./data/physeq_objects/ecto_microbiome_ps.RDS")

# merge and clean up for plotting
full <- merge_phyloseq(ecto_ps,endo_ps)
full <- 
  full %>% 
  subset_taxa(taxa_sums(full) >= 100) # remove very rare ASVs (N >= 100 observations in total)
full <- 
  full %>% 
  subset_samples(sample_sums(full) > 0)

# create new variable for merging
full@sam_data$mergevar <- paste(full@sam_data$g_spp,full@sam_data$location,sep="|")
# Export 
full %>% saveRDS("./data/physeq_objects/endo_ecto_microbiomes_ps.RDS")



# merge samples
full_m <- full %>% merge_samples('mergevar')
# repair metadata
full_m@sam_data$location <- sample_names(full_m) %>% str_split("\\|") %>% map_chr(2)
full_m@sam_data$g_spp <- sample_names(full_m) %>% str_split("\\|") %>% map_chr(1)
# transform compositionally
full_m <- full_m %>% microbiome::transform('compositional')


# plot
p <- 
full_m %>% 
  plot_bar2(fill="Phylum",x="location") +
  facet_wrap(~g_spp) +
  scale_fill_viridis_d(option = 'turbo') +
  labs(y="Composition",x="") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face='bold.italic',size=12),
        legend.title = element_text(face='bold',size=16),
        legend.text = element_text(face='bold',size=12),
        axis.title = element_text(face='bold',size=16),
        axis.text.x = element_text(face='bold',size=14,hjust = 0,vjust=.5,angle=90))
saveRDS(p,"./output/figs/phylum_composition_for_each_species_by_location.RDS")

# Plot richness for each species/treatment/isolate ####

# glom at species level
full_spp <- full %>% 
  tax_glom("Species",NArm = FALSE)

full


# add richness to sample data
rich <- full_spp %>% estimate_richness()
full_spp@sam_data$richness <- rich$Observed
  

dat <- microbiome::meta(full_spp)
dat$Genus_species <- paste0(dat$genus," ",dat$species)
dat$Genus_species <- factor(dat$Genus_species,levels = sort(unique(dat$Genus_species)))


p <- 
dat %>% 
  ggplot(data=.,
         aes(x=isolate,y=richness,fill=location)) +
  geom_boxplot() +
  facet_wrap(~Genus_species,scales='free_x') +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = 'bold.italic',size=8),
        axis.text.x = element_text(face='bold',size=8,angle=30,hjust=1,vjust=1),
        legend.title = element_text(face='bold'),
        legend.text = element_text(face='bold')) +
  scale_fill_manual(values=pal$pal.earthtones) +
  labs(y="Species richness",x="",fill="Microbiome\nlocation")
saveRDS(p,"./output/figs/species_richness_by_isolate_boxplot.RDS")



