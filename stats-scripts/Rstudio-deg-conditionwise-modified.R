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
cat("Repository root:", repo_dir, "\n")

#-------------------------------------------------------------------------------------

dir_peak<-readline("Enter directory path for the consensus read counts? Dont put in quotes #eg. consensus_wt\n :")

files <- list.files(paste0(dir_peak,'/count_reads'), full.names = FALSE)
file_names <- unique(
  tools::file_path_sans_ext(files)
)

file_names

########### load counts
count_list <- lapply(file_names, function(f) {
  fi<-paste0(dir_peak,'/count_reads/',f,'.txt')
  scan(fi)
})

# combine into matrix
counts <- as.data.frame(count_list)

# assign column names from prefix
colnames(counts) <- file_names

# assign peak IDs
rownames(counts) <- paste0("peak_", seq_len(nrow(counts)))

coldata <- data.frame(
  row.names = file_names,
  sample_prefix = file_names,
  condition= sub("^Input_", "", file_names),
  assay = ifelse(grepl("Input", file_names), "Input", "IP")
)

counts2<-counts+1
enrich_df <- data.frame(row.names = rownames(counts2))
for (ele in unique(coldata$condition)){
  print(ele)
  input_nm<-paste0('Input_',ele)
  enrich_df[[paste0(ele,"_enrich")]] <- counts2[[ele]] / counts2[[input_nm]]
}
dim(enrich_df)
colnames(enrich_df)
df<- counts2
peaks <- read.table(paste0(dir_peak,'/count_reads/',file_names[1],'.counts'))
df$chr<-peaks$V1
df$start<-peaks$V2
df$end<-peaks$V3
df$peakcondition<-peaks$V4
df$strand<-peaks$V5

############### Add region- 3UTR/coding/ 5UTR  and gene id/symbol etc.
## name them by conditions being tested. See example below


# IP conditions only
ip_conditions <- unique(coldata$condition[coldata$assay == "IP"])

cat("\nIP conditions present:\n")
print(ip_conditions)
cat("\nNumber of IP conditions:", length(ip_conditions), "\n")

############################################################
# Identify condition groups from the input names
############################################################

# Extract suffixes after the final underscore.
# Example: WT_water -> water; KO_stress -> stress
condition_suffix <- sub("^.*_", "", ip_conditions)

# Identify the available treatment/condition groups
treatment_groups <- unique(condition_suffix)

# Identify genotype names by removing the final suffix
genotypes <- unique(sub("_[^_]+$", "", ip_conditions))

############################################################
# 1. Compare treatment groups within each genotype
############################################################

for (genotype in genotypes) {
  
  genotype_conditions <- ip_conditions[
    grepl(paste0("^", genotype, "_"), ip_conditions)
  ]
  
  if (length(genotype_conditions) >= 2) {
    
    # Compare all available pairs for this genotype
    condition_pairs <- combn(genotype_conditions, 2, simplify = FALSE)
    
    for (pair in condition_pairs) {
      
      condition1 <- pair[1]
      condition2 <- pair[2]
      
      col1 <- paste0(condition1, "_enrich")
      col2 <- paste0(condition2, "_enrich")
      
      if (all(c(col1, col2) %in% colnames(enrich_df))) {
        
        comparison_name <- paste0(
          "log2FC_", condition2, "_vs_", condition1
        )
        
        df[[comparison_name]] <-
          log2(enrich_df[[col2]] / enrich_df[[col1]])
      }
    }
  }
}

############################################################
# 2. Compare genotypes within each treatment group
############################################################

for (treatment in treatment_groups) {
  
  conditions_in_group <- ip_conditions[
    condition_suffix == treatment
  ]
  
  # Use WT as reference if present; otherwise use the first
  # available genotype as reference
  reference <- paste0("WT_", treatment)
  
  if (!(reference %in% conditions_in_group)) {
    reference <- conditions_in_group[1]
  }
  
  reference_col <- paste0(reference, "_enrich")
  
  if (reference_col %in% colnames(enrich_df)) {
    
    for (condition in conditions_in_group) {
      
      if (condition != reference) {
        
        condition_col <- paste0(condition, "_enrich")
        
        if (condition_col %in% colnames(enrich_df)) {
          
          comparison_name <- paste0(
            "log2FC_", condition, "_vs_", reference
          )
          
          df[[comparison_name]] <-
            log2(enrich_df[[condition_col]] /
                   enrich_df[[reference_col]])
        }
      }
    }
  }
}

print(df)

##########################################

df_gr <- GRanges(
  seqnames = df$chr,
  ranges   = IRanges(start = df$start + 1, end = df$end),
  strand   = df$strand,
  peakcondition = df$peakcondition
)

# Add all log2FC columns automatically
log2FC_columns <- grep("^log2FC_", colnames(df), value = TRUE)

mcols(df_gr)[, log2FC_columns] <- df[, log2FC_columns, drop = FALSE]
head(df_gr)


gff3_file <- file.path(repo_dir,"genome_files", "GRCh38", "gencode.v39.annotation.gff3" )
# --- build TxDb from GFF3 --------------------------------------------------
txdb <- makeTxDbFromGFF(gff3_file, format = "gff3")
repeats_file <- file.path(repo_dir, "genome_files", "GRCh38", "RepeatMask.csv")
repeats<-read.csv(repeats_file)
source(file.path(repo_dir,'peaks-analysis/peak_analysis_functions.R'))
df_gr<-addregion_annotation(df_gr,txdb,repeats)
df_gr<-add_geneid_symbol_togranges(df_gr,txdb)
df_gr$gene_symbol[is.na(df_gr$gene_symbol)]<-''
df_gr$genesymbol_unique <- make.unique(as.character(df_gr$gene_symbol), sep = "_")
name<-readline("Enter file name prefix for deg results. eg. deg_wt\n :")
saveRDS(df_gr,paste0(dir_peak,'/',name,'.rds'))
write.xlsx(df_gr, paste0(dir_peak,'/',name,'.xlsx'))
