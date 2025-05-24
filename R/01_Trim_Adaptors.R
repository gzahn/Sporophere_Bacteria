# SETUP ####

## Packages ####
library(tidyverse); packageVersion('tidyverse')
library(dada2); packageVersion("dada2")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")
library(ShortRead); packageVersion("ShortRead")
library(parallel); packageVersion("parallel")

## Functions ####
source("./R/functions.R")

## Data ####
metadata <- read_csv("./data/metadata/sporosphere_metadata.csv")

# subset to existing files only
metadata <- metadata[file.exists(metadata$fwd_fp_raw) & file.exists(metadata$rev_fp_raw),]


## primer sequences ####

# 16S SSU v3-4
FWD_341f <- "CCTACGGGAGGCAGCAG"
REV_806r <- "GGACTACHVGGGTWTCTAAT"

fwd_names <- metadata$sample_id

fwd_filtn_names <- file.path("./data/raw/filtN",paste(fwd_names,"filtN_fwd.fastq.gz",sep="_"))
rev_filtn_names <- file.path("./data/raw/filtN",paste(fwd_names,"filtN_rev.fastq.gz",sep="_"))

# RUN CUTADAPT ####

remove_primers(metadata=metadata, # metadata object for multi-seq-run samples; must contain "run" column and fwd/rev filepath columns
               amplicon.colname = "amplicon", # column name that contains the amplicon info for each sample
               amplicon = "16S_V3V4", # which amplicon from the run are you processing (ITS, SSU, LSU, etc)?
               sampleid.colname = "sample_id", # column name in metadata containing unique sample identifier
               fwd.fp.colname = "fwd_fp_raw", # name of column in metadata indicating fwd filepath to raw data
               rev.fp.colname = "rev_fp_raw",
               fwd_pattern="_R1_",
               rev_pattern="_R2_",
               fwd_primer=FWD_341f,
               rev_primer=REV_806r,
               multithread=parallel::detectCores()-1)

