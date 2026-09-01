# SETUP ####

## packages ####
library(tidyverse)
library(phyloseq)
library(microbiome)
library(vegan)
library(corncob)
library(ranger)
library(vip)
library(Fidelity)
library(ggsci)
library(ggshroom)
library(ape)
library(glmmTMB)

## functions ####
source("./R/functions.R")
set.seed(666)

## get AMF phylogeny from TraitAM repository
if (!file.exists("./data/FinalTree_Oct2023.tre")) {
  download.file(
    "https://datadryad.org/downloads/file_stream/3985005",
    "./data/FinalTree_Oct2023.tre",
    method = "wget"
  )
}

## load data ####
full <- readRDS("./data/physeq_objects/endo_ecto_microbiomes_ps.RDS")
endo <- readRDS("./data/physeq_objects/endo_microbiome_ps.RDS")
ecto <- readRDS("./data/physeq_objects/ecto_microbiome_ps.RDS")

# extract phylogenetic tree from TraitAM file
x <- readLines("./data/FinalTree_Oct2023.tre")
tree_string <- x[grepl("^    TREE 1 =", x)]
# Everything after the "="
tree_string <- sub("^.*?=\\s*", "", tree_string)
# Remove BEAST/DendroPy annotations in square brackets
tree_string <- gsub("\\[&[^]]*\\]", "", tree_string)
# Read as Newick
amf_tree <- read.tree(text = tree_string)
length(amf_tree$tip.label)
plot(amf_tree, cex = 0.3)
# remove bootstrap values
amf_tree$node.label <- NULL

# check tree tip names
amf_tree$tip.label
# make them match phyloseq names
endo@sam_data$genus_species <- paste(
  endo@sam_data$genus,
  endo@sam_data$species,
  sep = "_"
)
ecto@sam_data$genus_species <- paste(
  ecto@sam_data$genus,
  ecto@sam_data$species,
  sep = "_"
)

study_species <- unique(endo@sam_data$genus_species)

tip_matches <- sapply(
  study_species,
  function(sp) {
    grep(
      paste0("^", sp, "_"),
      amf_tree$tip.label,
      value = TRUE
    )
  }
)

tip_matches
matched_tips <- unname(unlist(tip_matches))
amf_subtree <- ape::keep.tip(
  amf_tree,
  matched_tips
)
tip_map <- setNames(
  names(tip_matches),
  unname(unlist(tip_matches))
)

amf_subtree$tip.label <- tip_map[amf_subtree$tip.label]
setequal(
  amf_subtree$tip.label,
  study_species
)
saveRDS(
  amf_subtree,
  "./data/AMF_study_species_tree.rds"
)


phy <-
  full %>%
  tax_glom("Phylum", ) %>%
  psmelt()
phy$relabund <- phy$Abundance / sum(phy$Abundance)

## pie chart comparison w/ published metaanalysis phyla
pie_data <- phy %>%
  group_by(Phylum) %>%
  summarise(relabund = sum(relabund), .groups = "drop") %>%
  mutate(Phylum = if_else(relabund < 0.01, "Other", Phylum)) %>%
  group_by(Phylum) %>%
  summarise(relabund = sum(relabund), .groups = "drop")
ggplot(pie_data, aes(x = "", y = relabund, fill = Phylum)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  scale_fill_viridis_d(option = 'turbo') +
  theme(
    legend.title = element_text(face = 'bold', size = 16),
    legend.text = element_text(face = 'bold', size = 14)
  )
ggsave(
  "./output/figs/phylum_piechart_overall.png",
  dpi = 600,
  height = 6,
  width = 7
)

# ORDINATE #####

# ordinate and plot several ways
for (i in c("DCA", "CCA", "RDA", "NMDS", "MDS", "PCoA")) {
  x <- endo %>% microbiome::transform('pa')
  ord <- ordinate(
    x,
    method = i,
    distance = phyloseq::distance(x, method = "jaccard", binary = TRUE)
  )
  p <- plot_ordination(endo, ord, color = "trait_orn_height_mean") + ggtitle(i)
  print(p)
}

# PCoA is best
ord <- ordinate(
  ecto,
  method = "PCoA",
  distance = phyloseq::distance(x, method = "jaccard", binary = TRUE)
)
ord2 <- ordinate(
  endo,
  method = "PCoA",
  distance = phyloseq::distance(x, method = "jaccard", binary = TRUE)
)

x <- endo %>% microbiome::transform('pa')
ord <- ordinate(
  x,
  method = "PCoA",
  distance = phyloseq::distance(x, method = "jaccard", binary = TRUE)
)
plot_ordination(x, ord, color = "trait_orn_height_mean") +
  theme_bw() +
  scale_color_viridis_c(, end = .8) +
  labs(color = "Ornamentation\nheight (μm)\n", title = "Endospore") +
  theme(
    legend.title = element_text(face = 'bold', size = 12),
    legend.text = element_text(face = 'bold', size = 12),
    legend.position = "right",
    axis.title = element_text(face = 'bold', size = 14),
    plot.title = element_text(face = 'bold', size = 14)
  )
ggsave(
  "./output/figs/PCoA_plot_endospore.png",
  dpi = 600,
  height = 8,
  width = 8
)


x <- ecto %>% microbiome::transform('pa')
ord2 <- ordinate(
  x,
  method = "PCoA",
  distance = phyloseq::distance(x, method = "jaccard", binary = TRUE)
)
plot_ordination(x, ord2, color = "trait_orn_height_mean") +
  theme_bw() +
  scale_color_viridis_c(, end = .8) +
  labs(color = "Ornamentation\nheight (μm)\n", title = "Epispore") +
  theme(
    legend.title = element_text(face = 'bold', size = 12),
    legend.text = element_text(face = 'bold', size = 12),
    legend.position = "right",
    axis.title = element_text(face = 'bold', size = 14),
    plot.title = element_text(face = 'bold', size = 14)
  )
ggsave("./output/figs/PCoA_plot_epispore.png", dpi = 600, height = 8, width = 8)


env_data <- as(ecto@sam_data, 'data.frame') %>%
  dplyr::filter(sample_id %in% row.names(ord2$vectors)[1:nrow(ord2$values)]) %>%
  select(starts_with("trait_"))


# fit trait data to PCoA
env <- envfit(ord2$values, env_data, permutations = 999, na.rm = FALSE)


# export to file
sink("./output/envfit_results_table.txt")
print("Ectosporic")
env
sink(NULL)

env2_data <- as(endo@sam_data, 'data.frame') %>%
  dplyr::filter(sample_id %in% row.names(ord$vectors)[1:nrow(ord$values)]) %>%
  select(starts_with("trait_"))
env2 <- envfit(ord$values, env2_data, permutations = 999, na.rm = FALSE)
# export to file
sink("./output/envfit_results_table.txt", append = TRUE)
print("Endosporic")
env2
sink(NULL)


# PERMANOVA ####

## Prep phylogeny ####
phylo_dist <- cophenetic(amf_subtree)
phylo_pcoa <- ape::pcoa(phylo_dist)
phylo_pcoa$values
phylo_axes <- as.data.frame(phylo_pcoa$vectors[, 1:3])
colnames(phylo_axes) <- c(
  "phylo_PC1",
  "phylo_PC2",
  "phylo_PC3"
)
phylo_axes$genus_species <- rownames(phylo_axes)
endo_meta <- microbiome::meta(endo)
ecto_meta <- microbiome::meta(ecto)
endo_meta <- endo_meta %>%
  left_join(
    phylo_axes,
    by = "genus_species"
  )
ecto_meta <- ecto_meta %>%
  left_join(
    phylo_axes,
    by = "genus_species"
  )
endo_meta %>%
  select(genus_species, phylo_PC1, phylo_PC2, phylo_PC3) %>%
  distinct()

# check colinearity of phylogenetic decomposition with spore traits
cor_data <- endo_meta %>%
  select(
    genus_species,
    phylo_PC1,
    phylo_PC2,
    phylo_PC3,
    trait_orn_height_mean,
    trait_shape_median,
    trait_vol_mean,
    trait_investment_mean
  ) %>%
  distinct() %>%
  select(-genus_species) %>%
  cor()
write_csv(
  cor_data |> as.data.frame(),
  "./output/trait_phylogeny_correlation_matrix.csv"
)
cor_long <- as.data.frame(cor_data) %>%
  rownames_to_column("variable1") %>%
  pivot_longer(
    -variable1,
    names_to = "variable2",
    values_to = "correlation"
  )

ggplot(
  cor_long,
  aes(x = variable1, y = variable2, fill = correlation)
) +
  geom_tile() +
  geom_text(
    aes(label = sprintf("%.2f", correlation)),
    size = 3
  ) +
  scale_fill_gradient2(
    limits = c(-1, 1),
    midpoint = 0,
    name = "Pearson r"
  ) +
  coord_equal() +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
ggsave(
  "./output/figs/SI_trait_phylogeny_correlation.png",
  dpi = 300,
  height = 6,
  width = 6
)
# phylogeny and these traits are largely independent for this group of isolates!

## db-RDA/Variance partitioning ####
endo_ra <- transform_sample_counts(
  endo,
  function(x) x / sum(x)
)
ecto_ra <- transform_sample_counts(
  ecto,
  function(x) x / sum(x)
)
endo_mat <- as.matrix(otu_table(endo_ra))
if (taxa_are_rows(endo_ra)) {
  endo_mat <- t(endo_mat)
}
ecto_mat <- as.matrix(otu_table(ecto_ra))
if (taxa_are_rows(ecto_ra)) {
  ecto_mat <- t(ecto_mat)
}
isolate_id_endo <- as.character(endo@sam_data$isolate)

endo_mat_isolate <- rowsum(
  endo_mat,
  group = isolate_id_endo
)

rep_n_endo <- table(isolate_id_endo)

endo_mat_isolate <- sweep(
  endo_mat_isolate,
  1,
  rep_n_endo[rownames(endo_mat_isolate)],
  "/"
)

# Check
rowSums(endo_mat_isolate)

isolate_id_ecto <- as.character(ecto@sam_data$isolate)

ecto_mat_isolate <- rowsum(
  ecto_mat,
  group = isolate_id_ecto
)

rep_n_ecto <- table(isolate_id_ecto)

ecto_mat_isolate <- sweep(
  ecto_mat_isolate,
  1,
  rep_n_ecto[rownames(ecto_mat_isolate)],
  "/"
)

# Check
rowSums(ecto_mat_isolate)

endo_isolate_meta <- endo_meta %>%
  mutate(isolate = as.character(isolate)) %>%
  group_by(isolate) %>%
  summarise(
    genus_species = first(genus_species),
    genus = first(genus),
    trait_orn_height_mean = first(trait_orn_height_mean),
    trait_shape_median = first(trait_shape_median),
    trait_vol_mean = first(trait_vol_mean),
    trait_investment_mean = first(trait_investment_mean),
    phylo_PC1 = first(phylo_PC1),
    phylo_PC2 = first(phylo_PC2),
    phylo_PC3 = first(phylo_PC3),
    .groups = "drop"
  )
ecto_isolate_meta <- ecto_meta %>%
  mutate(isolate = as.character(isolate)) %>%
  group_by(isolate) %>%
  summarise(
    genus_species = first(genus_species),
    genus = first(genus),
    trait_orn_height_mean = first(trait_orn_height_mean),
    trait_shape_median = first(trait_shape_median),
    trait_vol_mean = first(trait_vol_mean),
    trait_investment_mean = first(trait_investment_mean),
    phylo_PC1 = first(phylo_PC1),
    phylo_PC2 = first(phylo_PC2),
    phylo_PC3 = first(phylo_PC3),
    .groups = "drop"
  )

endo_isolate_meta <- endo_isolate_meta[
  match(
    rownames(endo_mat_isolate),
    endo_isolate_meta$isolate
  ),
]

ecto_isolate_meta <- ecto_isolate_meta[
  match(
    rownames(ecto_mat_isolate),
    ecto_isolate_meta$isolate
  ),
]

# build dist matrices for dbRDA at isolate level
endo_dist_isolate <- vegdist(
  endo_mat_isolate,
  method = "bray"
)
ecto_dist_isolate <- vegdist(
  ecto_mat_isolate,
  method = "bray"
)

# initial adonis model
perm_endo_phylo <- adonis2(
  endo_dist_isolate ~
    phylo_PC1 +
    phylo_PC2 +
    phylo_PC3,
  data = endo_isolate_meta,
  by = "margin",
  permutations = 999
)

perm_endo_phylo

perm_ecto_phylo <- adonis2(
  ecto_dist_isolate ~
    phylo_PC1 +
    phylo_PC2 +
    phylo_PC3,
  data = ecto_isolate_meta,
  by = "margin",
  permutations = 9999
)

perm_ecto_phylo

# now overall model test, not partitioned phylogenetic axes:
perm_endo_phylo_global <- adonis2(
  endo_dist_isolate ~
    phylo_PC1 +
    phylo_PC2 +
    phylo_PC3,
  data = endo_isolate_meta,
  permutations = 9999
)
perm_ecto_phylo_global <- adonis2(
  ecto_dist_isolate ~
    phylo_PC1 +
    phylo_PC2 +
    phylo_PC3,
  data = endo_isolate_meta,
  permutations = 9999
)

perm_ecto_phylo_global
perm_endo_phylo_global

## PermANOVA ####
mat <- ecto %>%
  microbiome::transform('pa') %>%
  otu_table() %>%
  as('matrix')
perm_ecto <- adonis2(
  formula = mat ~ genus +
    trait_orn_height_mean +
    trait_shape_median +
    trait_vol_mean +
    trait_investment_mean,
  data = ecto@sam_data %>% as("data.frame"),
  by = "margin"
)
perm_ecto
adonis2(
  formula = mat ~ genus + family,
  data = ecto@sam_data %>% as("data.frame"),
  by = "margin"
)

mat <- endo %>%
  microbiome::transform('pa') %>%
  otu_table() %>%
  as('matrix')

perm_endo <- adonis2(
  formula = mat ~ genus +
    trait_orn_height_mean +
    trait_shape_median +
    trait_vol_mean +
    trait_investment_mean,
  data = endo@sam_data %>% as("data.frame"),
  by = "margin"
)
perm_endo
adonis2(
  formula = mat ~ genus + family,
  data = endo@sam_data %>% as("data.frame"),
  by = "margin"
)

mat <- full %>%
  microbiome::transform('pa') %>%
  otu_table() %>%
  as('matrix')
perm_full <- adonis2(
  formula = mat ~ genus * trait_orn_height_mean,
  data = full@sam_data %>% as("data.frame"),
  by = "margin"
)

perm_full

# save tables
a <- broom::tidy(perm_ecto) |> mutate(component = "Epispore")
b <- broom::tidy(perm_endo) |> mutate(component = "Endospore")
full_join(a, b) |>
  write_csv("./output/PermANOVA_results_tables_endo-ecto_df.csv")


# beta dispersion check
dist_endo <-
  endo %>%
  microbiome::transform('pa') %>%
  otu_table() %>%
  as('matrix') |>
  vegan::vegdist(method = 'jaccard', binary = TRUE)
betadisp_endo <- vegan::betadisper(dist_endo, endo@sam_data$genus)

dist_epi <-
  ecto %>%
  microbiome::transform('pa') %>%
  otu_table() %>%
  as('matrix') |>
  vegan::vegdist(method = 'jaccard', binary = TRUE)
betadisp_epi <- vegan::betadisper(dist_epi, ecto@sam_data$genus)

permutest(betadisp_endo, permutations = 999, pairwise = TRUE)
permutest(betadisp_epi, permutations = 999)

# dispersion values by genus
a <- data.frame(
  genus = betadisp_epi$group,
  distance = betadisp_epi$distances
) |>
  dplyr::group_by(genus) |>
  dplyr::summarise(
    mean_distance = mean(distance),
    sd_distance = sd(distance),
    n = dplyr::n()
  ) |>
  dplyr::mutate(component = "Epispore")

b <- data.frame(
  genus = betadisp_endo$group,
  distance = betadisp_endo$distances
) |>
  dplyr::group_by(genus) |>
  dplyr::summarise(
    mean_distance = mean(distance),
    sd_distance = sd(distance),
    n = dplyr::n()
  ) |>
  dplyr::mutate(component = "Endospore")
full_join(a, b)
write_csv(full_join(a, b), "./output/betadisp_by_host_genus.csv")
full_join(a, b) |>
  ggplot(aes(x = genus, y = mean_distance, fill = component)) +
  geom_col(position = 'dodge') +
  geom_errorbar(
    aes(ymin = mean_distance - sd_distance, ymax = mean_distance + sd_distance),
    position = 'dodge'
  ) +
  theme_bw() +
  labs(
    x = 'AMF genus',
    y = "Mean dispersion\n(distance to centroid)",
    fill = "Component"
  ) +
  scale_fill_manual(values = pal$pal.earthtones) +
  theme(
    axis.text.x = element_text(
      face = 'bold.italic',
      size = 12,
      angle = 60,
      hjust = 1
    ),
    axis.title = element_text(face = 'bold', size = 14),
    legend.title = element_text(face = 'bold', size = 14),
    legend.text = element_text(face = 'bold', size = 12),
    axis.text.y = element_text(face = 'bold', size = 12)
  )
ggsave(
  "./output/figs/betadisp_by_host_genus.png",
  dpi = 600,
  height = 8,
  width = 10
)


# Overlap test ####
# effect of ornamentation on epi/endo community overlap (genus level)
full_genus_melt <-
  full %>%
  tax_glom("Genus") %>%
  psmelt()

# sanity check
full_genus_melt %>%
  dplyr::filter(Abundance > 1) %>%
  group_by(isolate, trait_orn_height_mean, treatment, Genus) %>%
  summarize(N = n())

# genus-level overlap btwn epi/endo, using presence among ANY replicate
L <- list()
i <- 1
for (G in unique(full_genus_melt$Genus)) {
  (for (I in unique(full_genus_melt$isolate)) {
    C <- full_genus_melt %>%
      dplyr::filter(Abundance > 1 & Genus == G & isolate == I) %>%
      pluck("location") %>%
      unique() %>%
      length()

    L[[i]] <- data.frame(isolate = I, genus = G, count = C)
    i <- i + 1
  })
}
full_map <- reduce(L, full_join)


overlap_df <-
  full_map %>%
  full_join(data.frame(
    isolate = full@sam_data$isolate %>% as.character(),
    trait_orn_height_mean = full@sam_data$trait_orn_height_mean
  )) %>%
  unique.data.frame() %>%
  filter(count > 0) %>%
  mutate(
    overlap = case_when(
      count == 1 ~ FALSE,
      count == 2 ~ TRUE
    )
  ) %>%
  group_by(isolate, trait_orn_height_mean) %>%
  summarize(
    overlap = sum(overlap) %>% as.numeric(),
    .groups = "drop"
  )
overlap_df %>%
  summarise(
    mean = mean(overlap),
    variance = var(overlap),
    min = min(overlap),
    max = max(overlap)
  )
overlap_df <-
  full_map %>%
  filter(count > 0) %>%
  group_by(isolate) %>%
  summarise(
    shared = sum(count == 2),
    total = n(),
    proportion_shared = shared / total,
    .groups = "drop"
  ) %>%
  left_join(
    microbiome::meta(full) %>%
      select(isolate, trait_orn_height_mean) %>%
      distinct(),
    by = "isolate"
  )


overlap_df %>%
  ggplot(aes(x = trait_orn_height_mean, y = overlap)) +
  geom_point(size = 2) +
  theme_bw() +
  geom_smooth(method = 'lm', se = FALSE, color = 'black') +
  labs(y = "Scaled community overlap", x = "Ornamentation height (μm)") +
  theme(
    axis.title = element_text(face = 'bold', size = 14),
    axis.text = element_text(face = 'bold', size = 12)
  )


ggsave(
  "./output/figs/scaled_community_overlap_ornheight.png",
  dpi = 600,
  height = 6,
  width = 6
)

overlap_lm <- lm(
  overlap ~ trait_orn_height_mean,
  data = overlap_df
)
sink("./output/overlap_test.txt")
summary(overlap_lm)
sink(NULL)

overlap_df %>%
  summarise(
    mean = mean(overlap),
    variance = var(overlap),
    min = min(overlap),
    max = max(overlap)
  )

overlap_bb <- glmmTMB(
  cbind(shared, total - shared) ~ trait_orn_height_mean,
  data = overlap_df,
  family = betabinomial(link = "logit")
)

summary(overlap_bb)


full_genus_melt %>%
  dplyr::filter(
    Abundance > 1 & Genus == "Brevundimonas" & isolate == "MG101B-12"
  ) %>%
  pluck("location") %>%
  unique() %>%
  length()


# CORNCOB ####
# test for differential abundance by location (endo/ecto)
# watch for perfectly discriminatory taxa and remove them

ps_genus <- full %>% tax_glom("Genus", NArm = FALSE)

full@sam_data$genus
# account for AMF species or phylogeny here using genus as null
da_analysis <-
  differentialTest(
    formula = ~ location + genus, #abundance
    phi.formula = ~1, #dispersion
    formula_null = ~genus, # control for __ abundance
    phi.formula_null = ~1, # control for __ dispersion
    test = "Wald",
    boot = FALSE,
    B = 100,
    data = full,
    fdr_cutoff = 0.05,
    full_output = FALSE
  )

# saveRDS(da_analysis,"./output/da_analysis_location-g_spp.RDS")
da_analysis$significant_taxa
otu_to_taxonomy(da_analysis$significant_taxa, full)
plot(da_analysis)
ggsave("./output/figs/corncob_da.png", dpi = 600, height = 6, width = 12)

sink("./output/diffabund_taxa_full.txt")
da_analysis$significant_taxa
otu_to_taxonomy(da_analysis$significant_taxa, full)
da_analysis$significant_models
sink(NULL)


# RANDOM FOREST ####

# build data frame
dat <-
  data.frame(ectosporic = full@sam_data$location) %>%
  bind_cols(
    full %>%
      microbiome::transform('compositional') %>%
      otu_table() %>%
      as("matrix")
  ) %>%
  mutate(
    ectosporic = case_when(ectosporic == "Ectosporic" ~ TRUE, TRUE ~ FALSE)
  )

# fit model
rf_mod <-
  ranger::ranger(
    data = dat,
    formula = ectosporic ~ .,
    importance = 'permutation'
  )
# get most important variables (taxa)
vi(rf_mod)$Importance %>% plot
vi(rf_mod)$Importance |> summary()
vi_cutoff <- (vi(rf_mod)$Importance[vi(rf_mod)$Importance > 0] |>
  summary() |>
  pluck("Max.")) /
  2
# cutoff is half of max VI

important_taxa <- vi(rf_mod)[vi(rf_mod)$Importance >= vi_cutoff, ]
important_taxa$taxonomy <- otu_to_taxonomy(important_taxa$Variable, full)

# plot those between the locations
full@sam_data$location[full@sam_data$location == "Ectosporic"] <- "Episporic"

full %>%
  subset_taxa(taxa_names(full) %in% important_taxa$Variable) %>%
  merge_samples("location", fun = 'sum') %>%
  microbiome::transform('compositional') %>%
  plot_bar2(fill = "Family") +
  scale_fill_viridis_d(option = 'turbo') +
  labs(
    y = "Relative abundance",
    x = ""
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = 'bold.italic', size = 12),
    legend.title = element_text(face = 'bold', size = 16),
    legend.text = element_text(face = 'bold', size = 12),
    axis.title = element_text(face = 'bold', size = 16),
    axis.text.x = element_text(
      face = 'bold',
      size = 14,
      hjust = 1,
      vjust = .5,
      angle = 90
    )
  )
ggsave("./output/figs/important_taxa_relabund_by_location-family.png")

full %>%
  subset_taxa(taxa_names(full) %in% important_taxa$Variable) %>%
  merge_samples("location", fun = 'sum') %>%
  microbiome::transform('compositional') %>%
  plot_bar2(fill = "Genus") +
  scale_fill_viridis_d(option = 'turbo') +
  labs(
    y = "Relative abundance",
    x = "",
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = 'bold.italic', size = 12),
    legend.title = element_text(face = 'bold', size = 16),
    legend.text = element_text(face = 'bold', size = 12),
    axis.title = element_text(face = 'bold', size = 16),
    axis.text.x = element_text(
      face = 'bold',
      size = 14,
      hjust = 0,
      vjust = .5,
      angle = 90
    )
  )
ggsave("./output/figs/important_taxa_relabund_by_location-genus.png")

# FIDELITY ####
# run fidelity analysis
full_fidel <-
  full %>%
  Fidelity_physeq(groups = "location", ovp.plot = TRUE)
ggsave("./output/figs/fidelity_ovp_plot.png", dpi = 600, height = 5, width = 5)

csi <- full_fidel$community_index |> as.data.frame()
csi_aov <- aov(community_index ~ group, data = csi)
summary(csi_aov)
TukeyHSD(csi_aov) %>% plot()

sink("./output/community_specficity_aov.txt")
summary(csi_aov)
TukeyHSD(csi_aov)
sink(NULL)

# ectosporic communities are less host-specific
ggplot(csi, aes(x = group, y = community_index, fill = group)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = pal$pal.earthtones) +
  theme(
    axis.text = element_text(face = 'bold', size = 12),
    axis.title = element_text(face = 'bold', size = 14)
  ) +
  labs(x = "", y = "Community specificity index", fill = "Component") +
  theme(
    legend.title = element_text(face = 'bold', size = 12),
    legend.text = element_text(face = 'bold', size = 12),
    legend.position = "right",
    axis.text = element_text(face = 'bold', size = 14),
    plot.title = element_text(face = 'bold', size = 14)
  )
ggsave(
  "./output/figs/community_specficity.png",
  dpi = 600,
  height = 8,
  width = 8
)

# how host specificity changes with traits
spore_meta <- meta(full)
spore_meta_merge <- merge(csi, spore_meta)
spore_meta_merge$acaulospora_abund
csi_lm <-
  spore_meta_merge %>%
  dplyr::select(community_index, starts_with("trait_")) %>%
  lm(community_index ~ ., data = .)
summary(csi_lm)

spore_meta_merge %>%
  ggplot(aes(x = trait_color_most, y = community_index)) +
  geom_point() +
  geom_smooth(method = 'lm')

# We are interested in the host specificity of a few 'important' taxa
tsi <- full_fidel$taxon_index
mean_tsi <- mean(tsi$taxon_index)
tsi_important <-
  tsi %>%
  dplyr::filter(comm_name %in% important_taxa$Variable)
names_tsi_important <-
  important_taxa$taxonomy %>%
  str_extract("[^_]+_[^_]+$") %>%
  str_replace("_", " \\| ")
tsi_important$taxon_name <- names_tsi_important
tsi_important %>%
  ggplot(aes(x = taxon_name, y = taxon_index, color = most_loyal_to)) +
  geom_point(size = 5) +
  geom_hline(yintercept = mean_tsi, linetype = 2) +
  theme_bw() +
  scale_color_manual(values = pal$pal.earthtones) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = 'bold'),
    axis.text = element_text(face = 'bold', size = 12),
    axis.title = element_text(face = 'bold', size = 14),
    legend.title = element_text(face = 'bold', size = 14),
    legend.text = element_text(face = 'bold', size = 12)
  ) +
  labs(color = "Most 'loyal' to:", y = "Taxon index", x = "Taxonomy")
ggsave(
  "./output/figs/taxon_specificity_important_taxa.png",
  dpi = 600,
  height = 8,
  width = 12
)
