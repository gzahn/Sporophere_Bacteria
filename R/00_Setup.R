# SETUP ####
# sra-toolkit must be installed in your $PATH to run this

## set up custom PATH if needed
fqdump <- system("which fasterq-dump",intern = TRUE)
if(length(fqdump)==0){
  stop("Make sure fasterq-dump is installed (sra-toolkit) and in your $PATH before running this script")
}

## load packages ####
library(tidyverse)

## load metadata ####
dat <- read_csv("./data/metadata/sporosphere_metadata.csv")
# remove failed samples for this script
dat <- dat[!is.na(dat$sra_run),]

# DOWNLOADS ####

## Download Silva taxonomy databases ####
download.file(url = "https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toGenus_trainset.fa.gz?download=1",
              destfile = "./taxonomy/silva_nr99_v138.2_toGenus_trainset.fa.gz"
              )
download.file(url = "https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toSpecies_trainset.fa.gz?download=1",
              destfile = "./taxonomy/silva_nr99_v138.2_toSpecies_trainset.fa.gz"
)
download.file(url = "https://zenodo.org/records/14169026/files/silva_v138.2_assignSpecies.fa.gz?download=1",
              destfile = "./taxonomy/silva_v138.2_assignSpecies.fa.gz"
)

## Download raw sequence data from SRA ####

# get SRA accessions
acc <- dat$sra_run

# get infile names
sra_fwd <- paste0(dat$sra_run,"_1.fastq")
sra_rev <- paste0(dat$sra_run,"_2.fastq")

# get outfile names
fwd_fp <- dat$fwd_fp_raw
rev_fp <- dat$rev_fp_raw

# run download and extraction with sra-toolkit
for(i in seq_along(acc)){
  # get sra_run fwd and rev files
  system2(fqdump, args = acc[i])
  # compress w/ gzip
  system2("gzip", args = c(sra_fwd[i]))
  system2("gzip", args = c(sra_rev[i]))
  # rename and move files to ./data/raw/
  file.rename(paste0(sra_fwd[i],".gz"),fwd_fp[i])
  file.rename(paste0(sra_rev[i],".gz"),rev_fp[i])
}
