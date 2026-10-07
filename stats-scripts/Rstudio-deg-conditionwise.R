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
ip_conditions <- coldata$condition[coldata$assay == "IP"]
# remove duplicates
ip_conditions <- unique(ip_conditions)

cat("\nIP conditions present:\n")
print(ip_conditions)

cat("\nNumber of IP conditions:", length(ip_conditions), "\n")

############################################################
# 1. Stress vs water within the same genotype/condition
############################################################

water_conditions <- ip_conditions[grepl("_water$", ip_conditions)]
stress_conditions <- ip_conditions[grepl("_stress$", ip_conditions)]

for (water in water_conditions) {
  
  genotype <- sub("_water$", "", water)
  
  stress <- paste0(genotype, "_stress")
  
  if (stress %in% stress_conditions) {
    
    water_enrich <- paste0(water, "_enrich")
    stress_enrich <- paste0(stress, "_enrich")
    
    df[[paste0("log2FC_", genotype, "_stress_vs_water")]] <-
      log2(enrich_df[[stress_enrich]] / enrich_df[[water_enrich]])
  }
}

print(df)
############################################################
# 2. KO/genotype vs WT within water
############################################################

if ("WT_water" %in% ip_conditions) {

  for (condition in water_conditions) {

    if (condition != "WT_water") {

      genotype <- sub("_water$", "", condition)

      df[[paste0("log2FC_", genotype, "_vs_WT_water")]] <-
        log2(
          enrich_df[[paste0(condition, "_enrich")]] /
            enrich_df[["WT_water_enrich"]]
        )
    }
  }
}


############################################################
# 3. KO/genotype vs WT within stress
############################################################

if ("WT_stress" %in% ip_conditions) {

  for (condition in stress_conditions) {

    if (condition != "WT_stress") {

      genotype <- sub("_stress$", "", condition)

      df[[paste0("log2FC_", genotype, "_vs_WT_stress")]] <-
        log2(
          enrich_df[[paste0(condition, "_enrich")]] /
            enrich_df[["WT_stress_enrich"]]
        )
    }
  }
}

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

############################################# visualization plots ####################################
## move into a function file separately
library(openxlsx)
dir_peak<-readline("Enter directory path for the consensus read counts? Dont put in quotes #eg. consensus_wt\n :")
name<-readline("Enter file name prefix for deg results. used before eg. deg_wt\n :")
df<-read.xlsx(paste0(dir_peak,'/',name,'.xlsx'))

# ##scatter plot
# top_WT<-order(abs(df$log2FC_NSUN6KO_vs_WT_water),decreasing = TRUE)[1:100]
# top_NSUN6KO<-order(abs(df$log2FC_NSUN6KO_stress_vs_water),decreasing = TRUE)[1:100]
# 
# top_df_wt_nsun6ko<- dplyr::distinct(df[c(top_WT,top_NSUN6KO),])

source(file.path(repo_dir,'stats-scripts/plots_deg.R'))
## only wt or nsun6ko condition
plot_deg(df$log2FC_WT_stress_vs_water,seq_len(nrow(df)),"log2FC WT stress vs WT water",
         "Peaks","WT stress vs WT water",df$genesymbol_unique)
plot_deg(df$log2FC_NSUN6KO_stress_vs_water,seq_len(nrow(df)),"log2FC NSUN6KO stress vs water",
         "Peaks","NSUN6KO stress vs water",df$genesymbol_unique)

# nsun6_related<-df[abs(df$log2FC_NSUN6KO_vs_WT_water)> 3,]


## water samples deg
plot_deg(df$log2FC_NSUN6KO_vs_WT_water,seq_len(nrow(df)),"log2FC NSUN6KO vs WT water",
         "Peaks","NSUN6KO vs WT water",df$genesymbol_unique)
plot_deg(df$log2FC_NSUN2KO_vs_WT_water,seq_len(nrow(df)),"log2FC NSUN2KO vs WT water",
         "Peaks","NSUN2KO vs WT water",df$genesymbol_unique)
plot_deg(df$log2FC_DKO_vs_WT_water,seq_len(nrow(df)),"log2FC DKO vs WT water",
         "Peaks","DKO vs WT water",df$genesymbol_unique)

## stress samples deg
plot_deg(df$log2FC_NSUN6KO_vs_WT_stress,seq_len(nrow(df)),"log2FC NSUN6KO vs WT stress",
         "Peaks","NSUN6KO vs WT stress",df$genesymbol_unique)
plot_deg(df$log2FC_NSUN2KO_vs_WT_stress,seq_len(nrow(df)),"log2FC NSUN2KO vs WT stress",
         "Peaks","NSUN2KO vs WT stress",df$genesymbol_unique)
plot_deg(df$log2FC_DKO_vs_WT_stress,seq_len(nrow(df)),"log2FC DKO vs WT stress",
         "Peaks","DKO vs WT stress",df$genesymbol_unique)
## you can save the plots from export option or add ggsave command to save plots. ggsave saves last displayed plot so you have to run it after every plot. 

## only one medium (all water/ all stress, when common denominator used for normalization/ comparison eg. WT here for water or stress)


### for comparing withing conditions after normalizing also against WT of that condition (used eg. when comparing amongst all water/ all stress conditions)
## water
top_idx_nsun2 <- order(abs(df$log2FC_NSUN2KO_vs_WT_water),decreasing = TRUE)[1:100]
top_idx_nsun6 <- order(abs(df$log2FC_NSUN6KO_vs_WT_water),decreasing = TRUE)[1:100]
top_idx_dko<-order(abs(df$log2FC_DKO_vs_WT_water),decreasing = TRUE)[1:100]

top_df_nsun26<- dplyr::distinct(df[c(top_idx_nsun2,top_idx_nsun6),])
top_df_nsun6dko<-dplyr::distinct(df[c(top_idx_nsun6,top_idx_dko),])
top_df_nsun2dko<-dplyr::distinct(df[c(top_idx_nsun2,top_idx_dko),])

pdf(file = paste0(dir_peak,'/plots/nsun2ko_nsun6ko_topdeg.pdf'),
    width = 9, # The width of the plot in inches
    height = 8) # The height of the plot in inches

ggplot(df, aes(x = log2FC_NSUN2KO_vs_WT_water, y = log2FC_NSUN6KO_vs_WT_water)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_point(data = dplyr::distinct(top_df_nsun26, genesymbol_unique, .keep_all = TRUE),
             color = "red", size = 1.5) +
  geom_text_repel(
    data = dplyr::distinct(top_df_nsun26, genesymbol_unique, .keep_all = TRUE),
    aes(label = genesymbol_unique),
    size = 3,
    max.overlaps = 100
  ) +
  geom_hline(yintercept = 0, color = "red") +
  geom_vline(xintercept = 0, color = "red") +
  xlab("NSUN2 KO (log2FC)") +
  ylab("NSUN6 KO (log2FC)") +
  theme_bw()

dev.off()

pdf(file = paste0(dir_peak,'/plots/nsundko_nsun6ko_topdeg.pdf'),
    width = 9, # The width of the plot in inches
    height = 8) # The height of the plot in inches
ggplot(df, aes(x=log2FC_NSUN6KO_vs_WT_water, y=log2FC_DKO_vs_WT_water)) +
  geom_point(alpha=0.4, size=1) +
  geom_point(data=dplyr::distinct(top_df_nsun6dko, genesymbol_unique, .keep_all = TRUE), color="red", size=1.5) +
  geom_text_repel(
    data=dplyr::distinct(top_df_nsun6dko, genesymbol_unique, .keep_all = TRUE),
    aes(label=gene_symbol),
    size=3,
    max.overlaps=100
  ) +
  geom_hline(yintercept=0, color="red") +
  geom_vline(xintercept=0, color="red") +
  xlab("NSUN6 KO (log2FC)") +
  ylab("NSUN DKO (log2FC)") +
  theme_bw()
dev.off()

pdf(file = paste0(dir_peak,'/plots/nsundko_nsun2ko_topdeg.pdf'),
    width = 9, # The width of the plot in inches
    height = 8) # The height of the plot in inches
ggplot(df, aes(x=log2FC_NSUN2KO_vs_WT_water, y=log2FC_DKO_vs_WT_water)) +
  geom_point(alpha=0.4, size=1) +
  geom_point(data=dplyr::distinct(top_df_nsun2dko, genesymbol_unique, .keep_all = TRUE), color="red", size=1.5) +
  geom_text_repel(
    data=dplyr::distinct(top_df_nsun2dko, genesymbol_unique, .keep_all = TRUE),
    aes(label=gene_symbol),
    size=3,
    max.overlaps=100
  ) +
  geom_hline(yintercept=0, color="red") +
  geom_vline(xintercept=0, color="red") +
  xlab("NSUN2 KO (log2FC)") +
  ylab("NSUN DKO (log2FC)") +
  theme_bw()
dev.off()

top_df_nsun<- dplyr::distinct(df[c(top_idx_nsun2,top_idx_nsun6,top_idx_dko),])
write.xlsx(top_df_nsun, paste0(dir_peak,'/deg-allcondition-top100each.xlsx'))


##stress
top_idx_nsun2 <- order(abs(df$log2FC_NSUN2KO_vs_WT_stress),decreasing = TRUE)[1:100]
top_idx_nsun6 <- order(abs(df$log2FC_NSUN6KO_vs_WT_stress),decreasing = TRUE)[1:100]
top_idx_dko<-order(abs(df$log2FC_DKO_vs_WT_stress),decreasing = TRUE)[1:100]

top_df_nsun26<- dplyr::distinct(df[c(top_idx_nsun2,top_idx_nsun6),])
top_df_nsun6dko<-dplyr::distinct(df[c(top_idx_nsun6,top_idx_dko),])
top_df_nsun2dko<-dplyr::distinct(df[c(top_idx_nsun2,top_idx_dko),])

pdf(file = paste0(dir_peak,'/nsun2ko_nsun6ko_topdeg.pdf'),
    width = 9, # The width of the plot in inches
    height = 8) # The height of the plot in inches

ggplot(df, aes(x = log2FC_NSUN2KO_vs_WT_stress, y = log2FC_NSUN6KO_vs_WT_stress)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_point(data = dplyr::distinct(top_df_nsun26, genesymbol_unique, .keep_all = TRUE),
             color = "red", size = 1.5) +
  geom_text_repel(
    data = dplyr::distinct(top_df_nsun26, genesymbol_unique, .keep_all = TRUE),
    aes(label = genesymbol_unique),
    size = 3,
    max.overlaps = 100
  ) +
  geom_hline(yintercept = 0, color = "red") +
  geom_vline(xintercept = 0, color = "red") +
  xlab("NSUN2 KO (log2FC)") +
  ylab("NSUN6 KO (log2FC)") +
  theme_bw()

dev.off()

pdf(file = paste0(dir_peak,'/nsundko_nsun6ko_topdeg.pdf'),
    width = 9, # The width of the plot in inches
    height = 8) # The height of the plot in inches
ggplot(df, aes(x=log2FC_NSUN6KO_vs_WT_stress, y=log2FC_DKO_vs_WT_stress)) +
  geom_point(alpha=0.4, size=1) +
  geom_point(data=dplyr::distinct(top_df_nsun6dko, genesymbol_unique, .keep_all = TRUE), color="red", size=1.5) +
  geom_text_repel(
    data=dplyr::distinct(top_df_nsun6dko, genesymbol_unique, .keep_all = TRUE),
    aes(label=gene_symbol),
    size=3,
    max.overlaps=100
  ) +
  geom_hline(yintercept=0, color="red") +
  geom_vline(xintercept=0, color="red") +
  xlab("NSUN6 KO (log2FC)") +
  ylab("NSUN DKO (log2FC)") +
  theme_bw()
dev.off()

pdf(file = paste0(dir_peak,'/nsundko_nsun2ko_topdeg.pdf'),
    width = 9, # The width of the plot in inches
    height = 8) # The height of the plot in inches
ggplot(df, aes(x=log2FC_NSUN2KO_vs_WT_stress, y=log2FC_DKO_vs_WT_stress)) +
  geom_point(alpha=0.4, size=1) +
  geom_point(data=dplyr::distinct(top_df_nsun2dko, genesymbol_unique, .keep_all = TRUE), color="red", size=1.5) +
  geom_text_repel(
    data=dplyr::distinct(top_df_nsun2dko, genesymbol_unique, .keep_all = TRUE),
    aes(label=gene_symbol),
    size=3,
    max.overlaps=100
  ) +
  geom_hline(yintercept=0, color="red") +
  geom_vline(xintercept=0, color="red") +
  xlab("NSUN2 KO (log2FC)") +
  ylab("NSUN DKO (log2FC)") +
  theme_bw()
dev.off()

top_df_nsun<- dplyr::distinct(df[c(top_idx_nsun2,top_idx_nsun6,top_idx_dko),])
write.xlsx(top_df_nsun, paste0(dir_peak,'/deg-allcondition-top100each.xlsx'))


# ################################################
# ## Genes to highlight in GREEN
# ################################################
# 
# genes_green <- c(
#   "G3BP1",
#   "TIA1",
#   "RUNXL1","USO1"
# )
# 
# green_idx <- which(
#   df$gene_symbol %in% genes_green
# )
# 
# ## Plot green points
# points(
#   df$log2FC_WT_stress_vs_water[green_idx],
#   green_idx,
#   col = "green",
#   pch = 16,
#   cex = 1.2
# )
# 
# ## Label green genes
# text(
#   df$log2FC_WT_stress_vs_water[green_idx],
#   green_idx,
#   labels = df$gene_symbol[green_idx],
#   col = "green",
#   pos = 4,
#   cex = 0.7
# )
# 
# abline(v = 0, lty = 2)
