library(Guitar)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#######################################################################
################## Guitar Plot ########################################
guitarTxdb <- makeGuitarTxdb( txdb=txdb)
dir<-"~/Desktop/ARTR-seq/ezgi-data/2026-08-tia/peak-files-se-ext30"
# saveRDS(guitarTxdb,paste(dir,"guitartxdb.rds",sep='/'))

library(rtracklayer)

convert_narrowpeak <- function(infile, outfile){
  
  x <- read.table(
    infile,
    sep="\t",
    header=FALSE,
    stringsAsFactors=FALSE
  )
  
  colnames(x) <- c(
    "chr","start","end","name","score",
    "strand","signalValue",
    "pvalue","qvalue","peak"
  )
  
  gr <- GRanges(
    seqnames=x$chr,
    ranges=IRanges(
      start=x$start+1,
      end=x$end
    ),
    strand=x$strand
  )
  
  mcols(gr)$name <- x$name
  mcols(gr)$score <- x$signalValue
  
  common_chr <- intersect(
    seqlevels(gr),
    seqlevels(txdb)
  )
  
  gr <- keepSeqlevels(
    gr,
    common_chr,
    pruning.mode="coarse"
  )
  
  export(gr,outfile,format="BED")
}

convert_narrowpeak(
  paste(dir,"WT_water/peak.WT_water.narrowPeak.bed",sep='/'),
  paste(dir,"WT_water/WT_water.clean.bed",sep='/')
)

convert_narrowpeak(
  paste(dir,"NSUN2KO_water/peak.NSUN2KO_water.narrowPeak.bed",sep='/'),
  paste(dir,"NSUN2KO_water/NSUN2KO_water.clean.bed",sep='/')
)

convert_narrowpeak(
  paste(dir,"NSUN6KO_water/peak.NSUN6KO_water.narrowPeak.bed",sep='/'),
  paste(dir,"NSUN6KO_water/NSUN6KO_water.clean.bed",sep='/')
)

convert_narrowpeak(
  paste(dir,"DKO_water/peak.DKO_water.narrowPeak.bed",sep='/'),
  paste(dir,"DKO_water/DKO_water.clean.bed",sep='/'))


GuitarPlot(
  txTxdb = txdb,
    stBedFiles=list(
    WT=file.path(dir,"WT_TIA","WT.clean.bed"),
    NSUN2_KO=file.path(dir,"NSUN2KO","NSUN2KO.clean.bed"),
    NSUN6_KO=file.path(dir,"NSUN6KO","NSUN6KO.clean.bed"),
    DKO=file.path(dir,"DKO_TIA","DKO.clean.bed")
  ),
  pltTxType = c("mrna","ncrna"),
  enableCI=FALSE,adjust = 1
)


#######################################################################
