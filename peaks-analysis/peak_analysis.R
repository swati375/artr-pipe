## R program that creates a granges object of ARTR-seq peaks from macs3 narrowpeak files, add region of peak, gene names, etc.

library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(ChIPseeker)

prot <- readline("What is your protein name?")
##Import the crosslink sites from macs3 output
dir="~/Desktop/ARTR-seq/scripts_server/edited/"
setwd(dir)
peak_dir_name<-readline("What is your peak-file directory name?") #eg. peak-files_min30"
file <- readline("What is your ARTR-seq processed peak file name?") #eg. "peak.WT_TIA.narrowPeak.bed"
setwd(dir_peak)

## ---------------------------------------------------------
## Locate repository and input files
## ---------------------------------------------------------

## peak-analysis/ is the location of this script.
## The repository root is therefore one directory above it.
script_dir <- normalizePath(
  dirname(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[
    grep("^--file=", commandArgs(trailingOnly = FALSE))
  ])),
  mustWork = TRUE
)

repo_dir <- normalizePath(
  file.path(script_dir, ".."),
  mustWork = TRUE
)

dir_peak <- file.path(repo_dir, peak_dir_name)
peaks_file <- file.path(dir_peak, file)

if (!file.exists(peaks_file)) {
  stop("Peak file not found: ", peaks_file)
}

############### import crosslink sites from MACS3 output
name = sub("\\.narrowPeak\\.bed$", "", file)
print(name)
peak_df <- read.table(peaks_file, header = FALSE, sep = "\t")

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

source(file.path(script_dir,'peak_analysis_functions.R'))
peaks<-addregion_annotation(peaks,txdb,repeats)
peaks<-peak_annotate(name,txdb,dir_peak,peaks,repeats)
l=peak_go(name,dir_peak)
source(file.path(script_dir,'function_motif_seq_extract.R'))
## eg.
peak_motif_seq_extract(dir_peak,file,c('Intron'),'intron_')
peak_motif_seq_extract(dir_peak,file,c('3UTR','Intron'),'utr3_intron_')

##############check peak lengths in data using the hist_peaklength function in common functions in jupyter notebook- still to add ####################

peaks<-readRDS(paste(dir_peak,paste(name,".rds",sep=""),sep="/"))
peaks_genesym<-add_geneid_symbol_togranges(peaks)
saveRDS(peaks_genesym,paste(dir_peak,'/',name,'_geneid_sym.rds',sep=""))
##############################################################################
