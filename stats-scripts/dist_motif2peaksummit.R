rm(list=ls())
library(stringr)
library(GenomicFeatures)
library(openxlsx)

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
cat("Repository directory:", repo_dir, "\n")

dir_peakcount<-readline("Enter directory path for the consensus read counts? Dont put in quotes #eg. consensus_wt\n :")
name<-readline("Enter file name prefix for deg results. used before eg. deg_wt\n :")
df_gr<-readRDS(paste0(dir_peakcount,'/',name,'.rds'))

########################## add peak summit to allow distance from motif measurement

first_peak <- sub("\\|.*$", "", df_gr$peakcondition)
condition <- sub("\\.(fwd|rev)_peak_.*$", "", first_peak)
df_gr$peak_summit <- NA_integer_
peak_dir <- readline(
  prompt = "Enter directory containing peak folders (e.g. DKO_water, NSUN6KO_water): "
)

peak_dir <- normalizePath(peak_dir, mustWork = TRUE)

cat("Peak directory:", peak_dir, "\n\n")


for (cond in unique(condition)) {
  
  cat("Processing:", cond, "\n")
  
  # Rows belonging to this condition
  idx <- which(condition == cond)
  
  # Peak IDs for these rows
  peak_ids <- first_peak[idx]

  peak_file <- file.path(
    peak_dir,
    cond,
    paste0("peak.", cond, ".narrowPeak.bed")
  )
  
  # --------------------------------------------------------
  # Read MACS3 narrowPeak file
  #
  # Column 1  = chromosome
  # Column 2  = start (0-based BED coordinate)
  # Column 3  = end
  # Column 4  = peak name
  # Column 10 = summit offset
  # --------------------------------------------------------
  
  peaks <- read.delim(
    peak_file,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  
  # --------------------------------------------------------
  # Calculate GENOMIC summit coordinate
  # --------------------------------------------------------
  
  genomic_summit <- peaks[[2]] + peaks[[10]]
  
  # Create lookup table:
  # peak name -> genomic summit
  peak_summit_lookup <- setNames(
    genomic_summit,
    peaks[[4]]
  )
  
  # --------------------------------------------------------
  # Match the first peak in peakcondition
  # --------------------------------------------------------
  
  matched_summits <- peak_summit_lookup[peak_ids]
  # Check for peaks that could not be found
  missing <- is.na(matched_summits)
  
  if (any(missing)) {
    
    warning(
      sum(missing),
      " peak(s) could not be found in ",
      peak_file,
      " for condition ",
      cond
    )
    
    cat(
      "Missing peak IDs:\n",
      paste(unique(peak_ids[missing]), collapse = "\n"),
      "\n"
    )
  }
  
  # --------------------------------------------------------
  # Add genomic summit coordinates to df_gr
  # --------------------------------------------------------
  
  df_gr$peak_summit[idx] <- as.integer(matched_summits)
}

df_summit <- GRanges(
  seqnames = seqnames(df_gr),
  ranges = IRanges(
    start = df_gr$peak_summit,
    end = df_gr$peak_summit
  ),
  strand = strand(df_gr)
)

source(file.path(repo_dir,'peaks-analysis/function_motif_seq_extract.R'))
nsun6ko_lostpeak<-df_gr[df_gr$log2FC_NSUN6KO_vs_WT_water< -1.5,]
singlent_peak_seqlen_extract(dir_peakcount,nsun6ko_lostpeak,'nsun6ko_lost',800)
singlent_peak_seqlen_extract(dir_peakcount,nsun6ko_lostpeak,'nsun6ko_lost',400)
nsun6ko_gainpeak<-df_gr[df_gr$log2FC_NSUN6KO_vs_WT_water> 1.5,]
singlent_peak_seqlen_extract(dir_peakcount,nsun6ko_gainpeak,'nsun6ko_gain',800)
singlent_peak_seqlen_extract(dir_peakcount,nsun6ko_gainpeak,'nsun6ko_gain',400)

nsun2ko_lostpeak<-df_gr[df_gr$log2FC_NSUN2KO_vs_WT_water< -1.5,]
singlent_peak_seqlen_extract(dir_peakcount,nsun2ko_lostpeak,'nsun2ko_lost',800)
nsun2ko_gainpeak<-df_gr[df_gr$log2FC_NSUN2KO_vs_WT_water> 1.5,]
singlent_peak_seqlen_extract(dir_peakcount,nsun2ko_gainpeak,'nsun2ko_gain',800)

############################ motif-peak_summit distance ##########################3

library(Biostrings)
library(ggplot2)

fasta_file <- paste0(
  dir_peakcount,
  "/nsun6ko_lost_summit_800.fa"
)

motif_distance <- function(fasta_file, motif) {
  
  seqs <- readDNAStringSet(fasta_file)
  
  do.call(rbind, lapply(seq_along(seqs), function(i) {
    
    h <- matchPattern(
      DNAString(motif),
      seqs[[i]],
      fixed = FALSE
    )
    
    if (!length(h)) return(NULL)
    
    data.frame(
      sequence = i,
      distance = start(h) + (nchar(motif) - 1) / 2 -
        (nchar(seqs[[i]]) + 1) / 2,
      motif = motif
    )
  }))
}

motifs <- c(
  "CTCCA",
  "CTCTA",
  "CAGGG",
  "CTGGG",
  "CCGGG",
  "CGGGG",
  "GCATG"
)

colors <- c(
  "yellow",
  "red",
  "brown",
  "purple",
  "grey",
  "green",
  "blue"
)

motif_df <- do.call(
  rbind,
  lapply(motifs, function(m)
    motif_distance(fasta_file, m)
  )
)

############ if the motif occurs multiple times, all positions are taken.
ggplot(
  motif_df,
  aes(x = distance, colour = motif)
) +
  geom_density(linewidth = 1.2) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  scale_colour_manual(
    values = setNames(colors, motifs)
  ) +
  labs(
    x = "Distance of motif center to peak summit (nt)",
    y = "Density",
    colour = "Motif"
  ) +
  theme_classic()

# ## motif wise distance from peak summit
# library(tidyr)
# motif_table <- pivot_wider(
#   motif_df,
#   id_cols = sequence,
#   names_from = motif,
#   values_from = distance,
#   values_fn = \(x) paste(x, collapse = ", ")
# )


############ if the motif occurs multiple times, nearest motif taken only.
library(dplyr)
nearest_motif_df <- motif_df %>%
  group_by(sequence, motif) %>%
  slice_min(abs(distance), n = 1, with_ties = FALSE) %>%
  ungroup()

ggplot(
  nearest_motif_df,
  aes(x = distance, colour = motif)
) +
  geom_density(linewidth = 1.2) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  scale_colour_manual(
    values = setNames(colors, motifs)
  ) +
  labs(
    x = "Distance of nearest motif to peak summit (nt)",
    y = "Density",
    colour = "Motif"
  ) +
  theme_classic()

# library(tidyr)
# motif_table <- pivot_wider(
#   nearest_motif_df,
#   id_cols = sequence,
#   names_from = motif,
#   values_from = distance,
#   values_fn = \(x) paste(x, collapse = ", ")
# )
