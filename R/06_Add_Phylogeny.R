# BUILD PHYLOEGENY ####

# SETUP ####

## Packages and functions ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(phangorn); packageVersion("phangorn")
library(ape); packageVersion("ape")
library(seqinr); packageVersion("seqinr")
library(DECIPHER); packageVersion("DECIPHER")
library(parallel); packageVersion("parallel")
library(ShortRead); packageVersion("ShortRead")
# library(fastreeR); packageVersion("fastreeR")

## number of threads for parallel processing
nthreads <- parallel::detectCores() - 1

## Data ####
ps <- readRDS("./data/physeq_objects/clean_physeq_object.RDS")
invam <- readRDS("./data/physeq_objects/clean_invam_bulk_soil_physeq_object.RDS") 
neon <- readRDS("./data/physeq_objects/clean_neon_soil_physeq_object.RDS")

## Get Seqeuences ####
# simplify ASV names
seqs <- rownames(tax_table(ps))
names(seqs) <- paste0("ASV_",1:length(seqs)) # This propagates to the tip labels of the tree

# create DNAStringSet object
seqs_StringSet <- DNAStringSet(seqs)
ShortRead::writeFasta(DNAStringSet(seqs),"./taxonomy/16S_ASVs_seqs_ecto-endo.fasta")

# Multiple sequence alignment  ####
decipher_alignment <- DECIPHER::AlignSeqs(seqs_StringSet, processors = parallel::detectCores() - 1,verbose = TRUE) # DECIPHER method
adj_decipher_alignment <- DECIPHER::AdjustAlignment(decipher_alignment,processors = parallel::detectCores() - 1)

saveRDS(adj_decipher_alignment,"./taxonomy/16S_dna_alignment_DECIPHER.RDS")
ShortRead::writeFasta(adj_decipher_alignment,"./taxonomy/16S_ASVs_seqs_aligned.fasta")


# Convert to various formats
align_character <- as.character(decipher_alignment)

# phang.align <- as.phyDat(align_character, type = "DNA")
# dnabin.align <- as.DNAbin(adj_decipher_alignment)

# distance - maximum likelihood ####
dm <- DECIPHER::DistanceMatrix(adj_decipher_alignment)

#save
saveRDS(dm,"./taxonomy/16S_ML_Distance.RDS")
dm <- readRDS("./Output/16S_ML_Distance.RDS")
dist_dm <- as.dist(dm)

# Build tree outside of R
# https://mafft.cbrc.jp/alignment/server/spool/_nj.250619021810594.html?done

# import tree
x <- ape::read.tree(file = "./taxonomy/_nj.250619021810594.on.nh")

# replace tip labels with ASV sequences for import into phyloseq
x$tip.label <- x$tip.label %>% sub("^[0-9]+_", "",.)
x$tip.label <- seqs[x$tip.label]

#save

# add tree to phyloseq object ####
ps2 <- phyloseq(tax_table(tax_table(ps)),
                otu_table(otu_table(ps)),
                sample_data(sample_data(ps)),
                phy_tree(x))

# sanity check that tip labels are in right order
taxa_names(ps2)[1] == ps2@phy_tree$tip.label[1]

# Save updated phyloseq object with tree
saveRDS(ps2, "./data/physeq_objects/clean_physeq_object_w_tree.RDS")

