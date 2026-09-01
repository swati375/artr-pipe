################ extract sequences for motif analysis ##############################
## extracts sequences from hg38 genome

############### extract all peak sequences and select for type of region of interest to extract sequences
peak_motif_seq_extract<- function(dir,name,file,region_list,fname)
{

library(GenomicRanges)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)  # adjust genome if needed

peaks<-readRDS(paste(dir,paste(name,".rds",sep=""),sep="/"))
valid_chr <- seqnames(BSgenome.Hsapiens.UCSC.hg38)
peaks <- peaks[as.character(seqnames(peaks)) %in% valid_chr]
peaks_no_na <- peaks[!is.na(mcols(peaks)$distance_to_gene)]

seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, peaks_no_na)
names(seqs) <- paste0("peak_", seq_along(seqs))
names(seqs) <- make.unique(paste0(
  seqnames(peaks_no_na), ":",
  start(peaks_no_na), "-",
  end(peaks_no_na), "(",
  strand(peaks_no_na), ")"
))

## write to FASTA file
writeXStringSet(seqs, paste(dir,paste("peaks_",name,".fa",sep=""),sep="/"))

# 1. Filter peaks for region of interest
select_peaks <- peaks[mcols(peaks)$region %in% region_list]
##If a range in the GRanges object has the strand marked as "-", getSeq() automatically returns the reverse complement of that sequence.
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, select_peaks)

names(seqs) <- paste0("peak_", seq_along(seqs))
names(seqs) <- make.unique(paste0(
  seqnames(cds_3utr_peaks), ":",
  start(cds_3utr_peaks), "-",
  end(cds_3utr_peaks), "(",
  strand(cds_3utr_peaks), ")"
))

writeXStringSet(seqs, paste(dir,paste(fname,name,".fa",sep=""),sep="/"))
}


########################## extract sequences for selected peaks gives as a granges object
peak_seq_extract <- function(dir,peaks,fname)
{
    library(GenomicRanges)
    library(Biostrings)
    library(BSgenome.Hsapiens.UCSC.hg38)  # adjust genome if needed
    valid_chr <- seqnames(BSgenome.Hsapiens.UCSC.hg38)
  
    peaks <- peaks[as.character(seqnames(peaks)) %in% valid_chr]
    seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, peaks)
    names(seqs) <- paste0("peak_", seq_along(seqs))
    names(seqs) <- make.unique(paste0(
      seqnames(peaks), ":",
      start(peaks), "-",
      end(peaks), "(",
      strand(peaks), ")","_",peaks$genesymbol_unique
    ))
    print(names(seqs))
    print(paste0("no. of peaks=",length(peaks)))
    print(paste0("no. of seqs=",length(seqs)))
    writeXStringSet(seqs, paste(dir,paste(fname,".fa",sep=""),sep="/"))
}

####################################### extract +-10nt from peak summit
peak_seq_extract_summit <- function(dir,peaks,fname)
{
  library(GenomicRanges)
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38)  # adjust genome if needed
  
  valid_chr <- as.character(seqnames(BSgenome.Hsapiens.UCSC.hg38))
  
  peaks <- peaks[as.character(seqnames(peaks)) %in% valid_chr]
  new_start <- pmax(start(peaks), start(peaks) + peaks$summit - 10)
  new_end   <- pmin(end(peaks),   start(peaks) + peaks$summit + 10)
  peaks_centered<-GRanges(
    seqnames = seqnames(peaks),
    ranges   = IRanges(start = new_start, end = new_end),
    strand   = strand(peaks)
  )
  
  # optional: keep metadata
  mcols(peaks_centered) <- mcols(peaks)
  
  # trim just in case (chromosome boundaries)
  peaks_centered <- trim(peaks_centered)
  
  # extract sequences
  seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, peaks_centered)
  
  # naming
  names(seqs) <- make.unique(paste0(
    seqnames(peaks_centered), ":",
    start(peaks_centered), "-",
    end(peaks_centered), "(",
    strand(peaks_centered), ")"
  ))
  
  print(paste0("no. of peaks=",length(peaks)))
  print(paste0("no. of seqs=",length(seqs)))
  writeXStringSet(seqs, paste(dir,paste(fname,"_summit20.fa",sep=""),sep="/"))
}


############### for single nt positions, extract sequences +-len/2 from the position
singlent_peak_seqlen_extract <- function(dir,peaks,fname,len)
{
  library(GenomicRanges)
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38)  # adjust genome if needed
  
  valid_chr <- as.character(seqnames(BSgenome.Hsapiens.UCSC.hg38))
  peaks<-bs_gr
  len<-10
  peaks <- peaks[as.character(seqnames(peaks)) %in% valid_chr]
  new_start <- pmin(start(peaks) - (len/2), start(peaks) + (len/2))
  new_end   <- pmax(start(peaks) - (len/2), start(peaks) + (len/2))
  peaks_centered<-GRanges(
    seqnames = seqnames(peaks),
    ranges   = IRanges(start = new_start, end = new_end),
    strand   = strand(peaks)
  )
  
  # optional: keep metadata
  mcols(peaks_centered) <- mcols(peaks)
  
  # trim just in case (chromosome boundaries)
  peaks_centered <- trim(peaks_centered)
  
  # extract sequences
  seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, peaks_centered)
  
  # naming
  names(seqs) <- make.unique(paste0(
    seqnames(peaks_centered), ":",
    start(peaks_centered), "-",
    end(peaks_centered), "(",
    strand(peaks_centered), ")"
  ))
  
  print(paste0("no. of peaks=",length(peaks)))
  print(paste0("no. of seqs=",length(seqs)))
  writeXStringSet(seqs, paste(dir,paste(fname,"_summit_",len,".fa",sep=""),sep="/"))
}