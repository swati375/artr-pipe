## R program that creates a granges object of ARTR-seq peaks from macs3 narrowpeak files, add region of peak, gene names, etc.
# rm(list=ls())

## ---------------------------------------------------------
## Locate repository and input files
## ---------------------------------------------------------

if (interactive()) {
  
  ## Running from RStudio
  script_dir <- dirname(
    rstudioapi::getSourceEditorContext()$path
  )
  
} else {
  
  ## Running with Rscript
  script_path <- commandArgs(trailingOnly = FALSE)[
    grep("^--file=", commandArgs(trailingOnly = FALSE))
  ]
  
  script_dir <- dirname(
    sub("^--file=", "", script_path)
  )
}

script_dir <- normalizePath(
  script_dir,
  mustWork = TRUE
)

## peak-analysis/ is the location of this script.
## Repository root is therefore one directory above it.

repo_dir <- normalizePath(
  file.path(script_dir, ".."),
  mustWork = TRUE
)

cat("Script directory:", script_dir, "\n")
cat("Repository root:", repo_dir, "\n")
## ---------------------------------------------------------
## Get peak directory and peak file
## ---------------------------------------------------------

library(optparse)

option_list <- list(
  make_option(
    c("-p", "--peak_dir"),
    type = "character",
    help = "Path to peak-file directory"
  ),
  make_option(
    c("-f", "--file"),
    type = "character",
    help = "MACS3 narrowPeak file name"
  ),
  make_option(
    c("-P", "--protein"),
    type = "character",
    help = "Protein name"
  )
)

opt <- parse_args(
  OptionParser(option_list = option_list)
)

print(opt)

prot <- opt$protein
cat("Protein:", prot, "\n")
dir_peak <- path.expand(opt$peak_dir)
cat("Peak directory:", dir_peak, "\n")
file <- opt$file
cat("Peak file:", file, "\n")

peaks_file<-file.path(dir_peak,file)

if (!file.exists(peaks_file)) {
  stop("Peak file not found: ", peaks_file)
}
###--------------------------------
##create sample wise directory and move peak file there
##-------------------------------------------------

dir_extract <- sub(
  "^peak\\.(.*)\\.narrowPeak\\.bed$",
  "\\1",
  file
)

# dir_extract <- gsub("_", ".", dir_extract)
sample_dir <- file.path(
  dir_peak,
  dir_extract
)

dir.create(
  sample_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

new_peaks_file <- file.path(
  sample_dir,
  file
)

if (!file.rename(peaks_file, new_peaks_file)) {
  stop(
    "Could not move peak file from:\n",
    peaks_file, "\nto:\n",
    new_peaks_file
  )
}

# ############### import crosslink sites from MACS3 output
suppressPackageStartupMessages(
  suppressWarnings({
    library(GenomicFeatures)
    library(GenomicRanges)
    library(rtracklayer)
    library(ChIPseeker)
    library(optparse)

  })
)
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(ChIPseeker)
name = sub("\\.narrowPeak\\.bed$", "", file)
print(name)
peak_df <- read.table(new_peaks_file, header = FALSE, sep = "\t")

# narrowPeak spec has 10 columns (strand is not standard).
# If you have strand in col 11, add it here:
colnames(peak_df)[1:10] <- c("chr","start","end","name","score","strand",
                             "signalValue","pValue","qValue","peak")

# Convert to GRanges with strand
peaks <- GRanges(
  seqnames = peak_df$chr,
  ranges   = IRanges(start = peak_df$start+1, end = peak_df$end),
  strand   = peak_df$strand,
  summit = peak_df$peak,
  qval = peak_df$qValue
)

## files for GRCh38 given in genome"_files. For other gemnomes download from UCSC browser/ gencode. The repeats should be downloaded from Repeat Masker website.
gff3_file <- file.path(repo_dir,"genome_files", "GRCh38", "gencode.v39.annotation.gff3" )
# --- build TxDb from GFF3 --------------------------------------------------
txdb <- makeTxDbFromGFF(gff3_file, format = "gff3")
repeats_file <- file.path(repo_dir, "genome_files", "GRCh38", "RepeatMask.csv")
repeats<-read.csv(repeats_file)
source(file.path(script_dir,'peak_analysis_functions.R'))
peaks<-addregion_annotation(peaks,txdb,repeats)
peaks<-peak_annotate(name,txdb,sample_dir,peaks,repeats)
l=peak_go(name,sample_dir)
##------------------------------------------------------
## call for creating sequence file if needed
##--------------------------------------------
# source(file.path(script_dir,'function_motif_seq_extract.R'))
# ## eg.
# peak_motif_seq_extract(sample_dir,file,c('Intron'),'intron_')
# peak_motif_seq_extract(sample_dir,file,c('3UTR','Intron'),'utr3_intron_')

##############check peak lengths in data using the hist_peaklength function in common functions in jupyter notebook- still to add ####################

peaks<-readRDS(paste(sample_dir,paste(name,".rds",sep=""),sep="/"))
peaks_genesym<-add_geneid_symbol_togranges(peaks,txdb)
saveRDS(peaks_genesym,paste(sample_dir,'/',name,'_geneid_sym.rds',sep=""))
##############################################################################
