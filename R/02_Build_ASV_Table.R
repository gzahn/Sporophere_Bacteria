# SETUP ####

# packages
library(tidyverse)
library(phyloseq)
library(dada2)
library(decontam)

# functions
source("./R/functions.R")

# metadata
meta <- read_csv("./data/metadata/sporosphere_metadata.csv")

# add cutadapt file paths
cutadapt_dir <- list.dirs(full.names = TRUE)[grep("cutadapt$",list.dirs(full.names = TRUE))]
filtN_dir <- list.dirs(full.names = TRUE)[grep("filtN$",list.dirs(full.names = TRUE))]
meta$cutadapt_fwd_paths <- paste0(cutadapt_dir,"/",meta$sample_id,"_cutadapt_fwd.fastq.gz")
meta$cutadapt_rev_paths <- paste0(cutadapt_dir,"/",meta$sample_id,"_cutadapt_rev.fastq.gz")


# subset metadata to samples clearly present in cutadapt
meta <- meta[file.exists(meta$cutadapt_fwd_paths),]


# list of sequencing runs
# just one run for this project...
meta$run_id <- "run1"
all_runs <- unique(meta$run_id)


# RUN ON ALL SSU DATA ####
for(seqrun in all_runs){
  
  # make sure these seq runs actually have samples from that amplicon present
  n.samples.in.run <- sum(meta[['run_id']] == seqrun & meta[['amplicon']] == "16S_V3V4")
  if(n.samples.in.run > 0){
    build_asv_table(metadata=meta, # metadata object for multi-seq-run samples; must contain "run" column and fwd/rev filepath columns
                    run.id.colname = "run_id", # name of column in metadata indicating which sequencing run a sample comes from
                    run.id = seqrun, # the run ID to perform function on, from the run.id.colname column. Enter as a character, e.g., "1"
                    amplicon.colname = "amplicon", # column name that contains the amplicon info for each sample
                    amplicon = "16S_V3V4", # which amplicon from the run are you processing (ITS, SSU, LSU, etc)?
                    sampleid.colname = "sample_id", # column name in metadata containing unique sample identifier
                    fwd.fp.colname = "cutadapt_fwd_paths", # name of column in metadata indicating fwd filepath to trimmed data (e.g., cutadapt)
                    rev.fp.colname = "cutadapt_rev_paths", # name of column in metadata indicating rev filepath to trimmed data (e.g., cutadapt)
                    fwd.pattern = "_R1_", # pattern in filepath column indicating fwd reads (to be safe)
                    rev.pattern = "_R2_", # pattern in filepath column indicating rev reads (to be safe),
                    maxEE = c(3,3), # max expected errors for filtration step of dada2 (for single-end cases like ITS, will default to maxEE=2)
                    trim.right=NA, # amount to trim off of 3' end after quality truncation
                    truncQ = 2, # special value denoting "end of good quality sequence"
                    rm.phix = TRUE, # remove phiX sequences?
                    compress = TRUE, # gzip compression of output?
                    multithread = (parallel::detectCores() -2), # how many cores to use? Set to FALSE on windows
                    single.end = FALSE, # use only forward reads and skip rev reads and merging (e.g., for ITS data)?
                    filtered.dir = "filtered", # name of output directory for all QC filtered reads. will be created if not extant. subdirectory of trimmed filepath
                    asv.table.dir = "./data/ASV_Tables", # path to directory where final ASV table will be saved
                    random.seed = 666, # random seed
                    control.col = "treatment", # column name where negative controls are indicated
                    control.indicator = "negative_ctl" # rows matching this in control.col will be considered negative controls
    )
  } else {break}
}

# move files to custom location
clean_files_from <- list.files("./data/raw/cutadapt/filtered",full.names = TRUE)
clean_files_to <- clean_files_from %>% str_replace("raw/cutadapt/filtered/","clean/")
file.rename(clean_files_from,clean_files_to)
