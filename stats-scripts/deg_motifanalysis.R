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
############################################################# motif analysis ##############################################
## extract sequences for motif search
## NSUN6 motif: CTCCA/CTCTA- loops of hairpin- type II (3UTR)- hek293//// UCCA- hek293T
## NSUN2 motif: G-rich triplet at 5' end of stem loop structures- type I CNGGG

##nsun6ko vs wt water
dir_peak<-readline("Enter directory path for the consensus read counts? Dont put in quotes #eg. consensus_wt\n :")
name<-readline("Enter file name prefix for deg results. used before eg. deg_wt\n :")
df_gr<-readRDS(paste0(dir_peak,'/',name,'.rds'))

source(file.path(repo_dir,'peaks-analysis/function_motif_seq_extract.R'))
nsun6ko_lostpeak<-df_gr[df_gr$log2FC_NSUN6KO_vs_WT_water< -1.5,]
peak_seq_extract(dir_peak,nsun6ko_lostpeak,'nsun6ko_lost_water')
nsun6ko_gainpeak<-df_gr[df_gr$log2FC_NSUN6KO_vs_WT_water> 1.5,]
peak_seq_extract(dir_peak,nsun6ko_gainpeak,'nsun6ko_gain_water')
nsun6ko_midcloudpeak<-df_gr[abs(df_gr$log2FC_NSUN6KO_vs_WT_water)<= 1.5,]
peak_seq_extract(dir_peak,nsun6ko_midcloudpeak,'nsun6ko_midcloud_water')

nsun2ko_lostpeak<-df_gr[df_gr$log2FC_NSUN2KO_vs_WT_water< -1.5,]
peak_seq_extract(dir_peak,nsun2ko_lostpeak,'nsun2ko_lost_water')
nsun2ko_gainpeak<-df_gr[df_gr$log2FC_NSUN2KO_vs_WT_water> 1.5,]
peak_seq_extract(dir_peak,nsun2ko_gainpeak,'nsun2ko_gain_water')
nsun2ko_midcloudpeak<-df_gr[abs(df_gr$log2FC_NSUN2KO_vs_WT_water)<= 1.5,]
peak_seq_extract(dir_peak,nsun2ko_midcloudpeak,'nsun2ko_midcloud_water')

dko_lostpeak<-df_gr[df_gr$log2FC_DKO_vs_WT_water< -1.5,]
peak_seq_extract(dir_peak,dko_lostpeak,'dko_lost_water')
dko_gainpeak<-df_gr[df_gr$log2FC_DKO_vs_WT_water> 1.5,]
peak_seq_extract(dir_peak,dko_gainpeak,'dko_gain_water')
dko_midcloudpeak<-df_gr[abs(df_gr$log2FC_DKO_vs_WT_water)<= 1.5,]
peak_seq_extract(dir_peak,dko_midcloudpeak,'dko_midcloud_water')

##nsun6ko vs wt stress
dir_peak<-readline("Enter directory path for the consensus read counts? Dont put in quotes #eg. consensus_wt\n :")
name<-readline("Enter file name prefix for deg results. used before eg. deg_wt\n :")
df_gr<-readRDS(paste0(dir_peak,'/',name,'.rds'))

source(file.path(repo_dir,'peaks-analysis/function_motif_seq_extract.R'))
nsun6ko_lostpeak<-df_gr[df_gr$log2FC_NSUN6KO_vs_WT_stress< -1.5,]
peak_seq_extract(dir_peak,nsun6ko_lostpeak,'nsun6ko_lost_stress')
singlent_peak_seqlen_extract(dir_peak,nsun6ko_lostpeak,'nsun6ko_lost',400)
nsun6ko_gainpeak<-df_gr[df_gr$log2FC_NSUN6KO_vs_WT_stress> 1.5,]
peak_seq_extract(dir_peak,nsun6ko_gainpeak,'nsun6ko_gain_stress')
nsun6ko_midcloudpeak<-df_gr[abs(df_gr$log2FC_NSUN6KO_vs_WT_stress)<= 1.5,]
peak_seq_extract(dir_peak,nsun6ko_midcloudpeak,'nsun6ko_midcloud_stress')

nsun2ko_lostpeak<-df_gr[df_gr$log2FC_NSUN2KO_vs_WT_stress< -1.5,]
peak_seq_extract(dir_peak,nsun2ko_lostpeak,'nsun2ko_lost_stress')
nsun2ko_gainpeak<-df_gr[df_gr$log2FC_NSUN2KO_vs_WT_stress> 1.5,]
peak_seq_extract(dir_peak,nsun2ko_gainpeak,'nsun2ko_gain_stress')
nsun2ko_midcloudpeak<-df_gr[abs(df_gr$log2FC_NSUN2KO_vs_WT_stress)<= 1.5,]
peak_seq_extract(dir_peak,nsun2ko_midcloudpeak,'nsun2ko_midcloud_stress')

dko_lostpeak<-df_gr[df_gr$log2FC_DKO_vs_WT_stress< -1.5,]
peak_seq_extract(dir_peak,dko_lostpeak,'dko_lost_stress')
dko_gainpeak<-df_gr[df_gr$log2FC_DKO_vs_WT_stress> 1.5,]
peak_seq_extract(dir_peak,dko_gainpeak,'dko_gain_stress')
dko_midcloudpeak<-df_gr[abs(df_gr$log2FC_DKO_vs_WT_stress)<= 1.5,]
peak_seq_extract(dir_peak,dko_midcloudpeak,'dko_midcloud_stress')


#### finding motif :
#use environment environment_motif
#run "~/Desktop/ARTR-seq/scripts_server/edited/stats-scripts/getmotif_fromfastaseq.sh"
source(file.path(repo_dir,'stats-scripts/plots_deg.R'))
p<-plot_piechart(as.data.frame(table(nsun6ko_lostpeak$region)),'nsun6ko_lost_stress')
p
p<-plot_piechart(as.data.frame(table(nsun6ko_gainpeak$region)),'nsun6ko_gain_stress')
p
p<-plot_piechart(as.data.frame(table(nsun6ko_midcloudpeak$region)),'nsun6ko_midcloud_stress')
p


########################
# select only 3utr sequences
dfgr_water<-readRDS("~/Desktop/ARTR-seq/ezgi-data/2026-08-tia/consensus_water/deg_water.rds")
nsun6_lost_water<-dfgr_water[dfgr_water$log2FC_NSUN6KO_vs_WT_water< -1.5,]
source(file.path(repo_dir,'stats-scripts/plots_deg.R'))
p<-plot_piechart(as.data.frame(table(nsun6_lost_water$region)),'nsun6ko_lost_water')
p
nsun6_gain_water<-dfgr_water[dfgr_water$log2FC_NSUN6KO_vs_WT_water> 1.5,]
p<-plot_piechart(as.data.frame(table(nsun6_gain_water$region)),'nsun6ko_gain_water')
p
nsun6_midcloud_water<-dfgr_water[abs(dfgr_water$log2FC_NSUN6KO_vs_WT_water)< 1.5,]
p<-plot_piechart(as.data.frame(table(nsun6_midcloud_water$region)),'nsun6ko_midcloud_water')
p
nsun6_lost_water_3utr<-nsun6_lost_water[nsun6_lost_water$region=='3UTR']
source(file.path(repo_dir,'peaks-analysis/function_motif_seq_extract.R'))
peak_seq_extract(dir_peak,nsun6_lost_water_3utr,'nsun6ko_lost_water_3utr')
singlent_peak_seqlen_extract(dir_peak,nsun6_lost_water_3utr,'nsun6ko_lost_3utr',400)

nsun6_lost_water_intron<-nsun6_lost_water[nsun6_lost_water$region=='Intron']
source(file.path(repo_dir,'peaks-analysis/function_motif_seq_extract.R'))
peak_seq_extract(dir_peak,nsun6_lost_water_intron,'nsun6ko_lost_water_intron')
singlent_peak_seqlen_extract(dir_peak,nsun6_lost_water_intron,'nsun6ko_lost_intron',400)
########################3
#Peaks/ gene with peaks lost only caz of KO and not stress

dfgr_water<-readRDS("~/Desktop/ARTR-seq/ezgi-data/2026-08-tia/consensus_water/deg_water.rds")
dfgr_stress<-readRDS("~/Desktop/ARTR-seq/ezgi-data/2026-08-tia/consensus_stress/deg_stress.rds")
# nsun6_lost_water<-dfgr_water[dfgr_water$log2FC_NSUN6KO_vs_WT_water< -1.5,]
# nsun6_gain_water<-dfgr_water[dfgr_water$log2FC_NSUN6KO_vs_WT_water> 1.5,]
# nsun6_mid_water<-dfgr_water[abs(dfgr_water$log2FC_NSUN6KO_vs_WT_water)< 1.5,]
# nsun6_lost_stress<-dfgr_stress[dfgr_stress$log2FC_NSUN6KO_vs_WT_stress< -1.5,]
# nsun6_gain_stress<-dfgr_stress[dfgr_stress$log2FC_NSUN6KO_vs_WT_stress> 1.5,]
# nsun6_mid_stress<-dfgr_stress[abs(dfgr_stress$log2FC_NSUN6KO_vs_WT_stress)< 1.5,]

# nsun6lost_stress_not_inwater <- nsun6_lost_stress[
#   !is.na(nsun6_lost_stress$gene_id) &
#     !nsun6_lost_stress$gene_id %in% nsun6_lost_water$gene_id
# ]
# ## extract sequence for these and check motif
# source(file.path(repo_dir,'peaks-analysis/function_motif_seq_extract.R'))
# peak_seq_extract(dir_peak,nsun6ko_lostpeak,'nsun6ko_lost_stress_notinwtstress')


## go analysis
library(clusterProfiler)
library(org.Hs.eg.db)
# Extract unique gene symbols
genes_unique <- unique(na.omit(nsun6_lost_water$gene_symbol))

# Convert gene symbols to Entrez IDs
gene_ids <- bitr(
  genes_unique,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# GO enrichment
ego <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

dotplot(ego, showCategory = 20)
