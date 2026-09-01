add_geneid_symbol_togranges <- function (obj,txdb)
{
  library(GenomicFeatures)
  # Gene assignement more that Finding nearest gene
  gene_ranges <- genes(txdb)
  gene_idx <- nearest(obj, gene_ranges, select="arbitrary")
  # Safely subset: only non-NA indices
  gene_id <- rep(NA, length(obj))
  gene_id[!is.na(gene_idx)] <- names(gene_ranges)[gene_idx[!is.na(gene_idx)]]
  mcols(obj)$gene_id <- gene_id
  
  ##convert geneid to gene symbol
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  
  gene_ids_clean <- sub("\\..*$", "", obj$gene_id)
  # Map gene symbols to Entrez IDs for enrichment
  
  gene_symbols <- mapIds(
    org.Hs.eg.db,
    keys = gene_ids_clean,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  obj$gene_symbol <- gene_symbols
  return(obj)
}

addregion_annotation <- function(par_gr,txdb,repeats) {
  
  library(GenomicFeatures)
  library(GenomicRanges)
  library(rtracklayer)
  library(ChIPseeker)
  library(ggplot2)

  colnames(repeats)[colnames(repeats) == 'genoName'] <- 'csv_chr_col'
  colnames(repeats)[colnames(repeats) == 'genoStart'] <- 'csv_start_col'
  colnames(repeats)[colnames(repeats) == 'genoEnd'] <- 'csv_end_col'
  csv_class_col <- "repClass"  # name of the column in your repeats dataframe
  
  # extract feature GRanges (unlist GRangesList -> GRanges)
  utr5_gr    <- unlist(fiveUTRsByTranscript(txdb), use.names = FALSE)
  utr3_gr    <- unlist(threeUTRsByTranscript(txdb), use.names = FALSE)
  cds_gr     <- unlist(cdsBy(txdb, by = "tx"), use.names = FALSE)
  exons_gr   <- unlist(exonsBy(txdb, by = "tx"), use.names = FALSE)
  introns_gr <- unlist(intronsByTranscript(txdb, use.names = FALSE), use.names = FALSE)
  
  # create noncoding exons = exons that do not overlap any CDS
  hit_exon_cds <- findOverlaps(exons_gr, cds_gr, ignore.strand = TRUE)
  noncoding_exons <- exons_gr[-unique(queryHits(hit_exon_cds))]
  # ensure seqname style matches (e.g., "chr1" vs "1"). We'll coerce feature seqlevels to the same style as peaks.
  # If peaks use "chr" prefix, convert txdb features to UCSC style:
  seqlevelsStyle(utr5_gr) <- seqlevelsStyle(par_gr)
  seqlevelsStyle(utr3_gr) <- seqlevelsStyle(par_gr)
  seqlevelsStyle(cds_gr)  <- seqlevelsStyle(par_gr)
  seqlevelsStyle(exons_gr)<- seqlevelsStyle(par_gr)
  seqlevelsStyle(introns_gr) <- seqlevelsStyle(par_gr)
  seqlevelsStyle(noncoding_exons) <- seqlevelsStyle(par_gr)
  
  repeat_gr <- GRanges(
    seqnames = repeats[["csv_chr_col"]],
    ranges   = IRanges(start = as.integer(repeats[["csv_start_col"]]),
                       end   = as.integer(repeats[["csv_end_col"]])),
    strand   = if("strand" %in% colnames(repeats)) repeats[["strand"]] else "*"
  )
  # attach class/family if available
  if (csv_class_col %in% colnames(repeats)) {
    mcols(repeat_gr)$repClass <- repeats[[csv_class_col]]
  }
  seqlevelsStyle(repeat_gr) <- seqlevelsStyle(par_gr)
  
  
  # --- fast vectorized annotation -------------------------------------------
  n <- length(par_gr)
  region <- rep("Intergenic", n)   # default
  
  # helper to assign only peaks still "Intergenic"
  assign_region <- function(feature_gr, label) {
    if (length(feature_gr) == 0) return(invisible(NULL))
    ol <- findOverlaps(par_gr, feature_gr, ignore.strand = TRUE)
    if (length(ol) == 0) return(invisible(NULL))
    idx <- unique(queryHits(ol))
    to_change <- idx[ region[idx] == "Intergenic" ]
    if (length(to_change)) region[to_change] <<- label
    invisible(NULL)
  }
  
  # priority order: UTR5 > UTR3 > CDS > Noncoding_exon > Intron > Repeat > Intergenic
  assign_region(utr5_gr, "5UTR")
  assign_region(utr3_gr, "3UTR")
  assign_region(cds_gr,  "CDS")
  assign_region(noncoding_exons, "Noncoding_exon")
  assign_region(introns_gr, "Intron")
  assign_region(repeat_gr, "Repeat")
  
  # add region metadata to peaks
  mcols(par_gr)$region <- region
  print(table(region))
  
  return(par_gr)
}


peak_annotate <- function (name,txdb,dir,peaks,repeats)
{
  # optionally add nearest gene / gene id using TxDb
  # has cds, 3UTR,5UTR, intron,exon
  # so nearest will return the same gene where it occurs if its in these regions of the gene.
  gene_ranges <- genes(txdb)
  
  # Peak center
  peak_center <- GRanges(seqnames = seqnames(peaks),
                         ranges = IRanges(
                           start = start(peaks) + peaks$summit,
                           width = 1),
                         strand = strand(peaks))
  
  # Gene assignement more that Finding nearest gene
  nearest_gene_idx <- nearest(peak_center, gene_ranges, select="arbitrary")
  # Safely subset: only non-NA indices
  nearest_genes <- rep(NA, length(peak_center))
  nearest_genes[!is.na(nearest_gene_idx)] <- names(gene_ranges)[nearest_gene_idx[!is.na(nearest_gene_idx)]]
  
  # Initialize distances
  distances <- rep(NA, length(peak_center))
  
  for (i in seq_along(peak_center)) {
    if (is.na(nearest_gene_idx[i])) {
      distances[i] <- NA
      next
    }
    
    pc <- peak_center[i]
    gene <- gene_ranges[nearest_gene_idx[i]]
    
    if (countOverlaps(pc, gene) > 0) {
      distances[i] <- 0
    } else {
      # Gene strand aware
      if (as.character(strand(gene)) == "+") {
        boundary <- ifelse(start(pc) < start(gene), start(gene), end(gene))
      } else {
        boundary <- ifelse(start(pc) < start(gene), end(gene), start(gene))
      }
      distances[i] <- start(pc) - boundary
    }
  }
  
  # Add back to peaks
  mcols(peaks)$nearest_gene <- nearest_genes
  mcols(peaks)$distance_to_gene <- distances
  
  # Make a data frame with region counts
  region_counts <- as.data.frame(table(peaks$region))
  colnames(region_counts) <- c("Region", "Count")
  
  # Calculate percentage labels
  region_counts$Perc <- round(region_counts$Count / sum(region_counts$Count) * 100, 1)
  region_counts$Label <- paste0(region_counts$Region, " (", region_counts$Perc, "%)")
  
  # Pie chart with ggplot2
  ggplot(region_counts, aes(x = "", y = Count, fill = Region)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    theme_void() +
    labs(title = paste0("Peak annotation distribution for ",name)) +
    geom_text(aes(label = Perc), position = position_stack(vjust = 0.5), size = 4)
  
  ggsave(file.path(dir,paste0(name,"_regions.pdf")),width=5,height=5,dpi=300)
  
  ############################## plot with extended intergenic and intron labels
  # promoters ±2kb around TSS
  tss_gr <- promoters(genes(txdb), upstream = 2000, downstream = 2000)
  
  # TTS ±2kb: get 1bp range at transcript end, then extend ±2kb
  tts_ends <- resize(genes(txdb), width=1, fix="end")
  tts_gr <- promoters(tts_ends, upstream=2000, downstream=2000)
  
  # Match seqlevels style to peaks
  seqlevelsStyle(tss_gr) <- seqlevelsStyle(peaks)
  seqlevelsStyle(tts_gr) <- seqlevelsStyle(peaks)
  
  ### Assign subcategories for Intergenic and Intron peaks
  # Initialize sub_category vector with NA
  sub_category <- rep(NA_character_, length(peaks))
  
  # Filter indices of Intergenic and Intron peaks
  focus_idx <- which(peaks$region %in% c("Intergenic", "Intron"))
  
  # Subset these peaks for overlap queries
  focus_peaks <- peaks[focus_idx]
  
  colnames(repeats)[colnames(repeats) == 'genoName'] <- 'csv_chr_col'
  colnames(repeats)[colnames(repeats) == 'genoStart'] <- 'csv_start_col'
  colnames(repeats)[colnames(repeats) == 'genoEnd'] <- 'csv_end_col'
  csv_class_col <- "repClass"  # name of the column in your repeats dataframe
  repeat_gr <- GRanges(
    seqnames = repeats[["csv_chr_col"]],
    ranges   = IRanges(start = as.integer(repeats[["csv_start_col"]]),
                       end   = as.integer(repeats[["csv_end_col"]])),
    strand   = if("strand" %in% colnames(repeats)) repeats[["strand"]] else "*"
  )
  # attach class/family if available
  if (csv_class_col %in% colnames(repeats)) {
    mcols(repeat_gr)$repClass <- repeats[[csv_class_col]]
  }
  seqlevelsStyle(repeat_gr) <- seqlevelsStyle(peaks)
  
  # Helper to assign subcategories based on overlap
  assign_subcat <- function(feature_gr, label) {
    ol <- findOverlaps(focus_peaks, feature_gr, ignore.strand=TRUE)
    hits <- unique(queryHits(ol))
    if(length(hits)) {
      sub_category[focus_idx[hits]] <<- label
    }
  }
  
  # Assign repeat classes from RepeatMasker (LINE, SINE, DNA, LTR)
  LINE_gr <- repeat_gr[grep("^LINE", mcols(repeat_gr)$repClass, ignore.case=TRUE)]
  SINE_gr <- repeat_gr[grep("^SINE", mcols(repeat_gr)$repClass, ignore.case=TRUE)]
  DNA_gr  <- repeat_gr[grep("^DNA",  mcols(repeat_gr)$repClass, ignore.case=TRUE)]
  LTR_gr  <- repeat_gr[grep("^LTR",  mcols(repeat_gr)$repClass, ignore.case=TRUE)]
  
  assign_subcat(LINE_gr, "LINE")
  assign_subcat(SINE_gr, "SINE")
  assign_subcat(DNA_gr,  "DNA")
  assign_subcat(LTR_gr,  "LTR")
  
  # Assign promoter and TTS regions
  assign_subcat(tss_gr, "Promoter")
  assign_subcat(tts_gr, "TTS")
  
  # For any remaining NA in sub_category at these indices, assign original broad label
  sub_category[focus_idx[is.na(sub_category[focus_idx])]] <- peaks$region[focus_idx[is.na(sub_category[focus_idx])]]
  
  # Add to peaks metadata
  mcols(peaks)$sub_category <- sub_category
  table(peaks$sub_category)
  table(peaks$region)
  
  ## Make pie chart of subcategories only for Intergenic and Intron peaks
  library(ggplot2)
  
  # Filter only Intergenic/Intron peaks for plotting
  plot_peaks <- peaks[mcols(peaks)$region %in% c("Intergenic", "Intron")]
  
  # Make table of subcategories (should be no NAs now)
  subcat_table <- table(mcols(plot_peaks)$sub_category)
  if (length(subcat_table)==0){
    # print('yes')
    saveRDS(peaks,paste(dir,'/',name,'.rds',sep=""))
    return(peaks)
  }
  subcat_df <- as.data.frame(subcat_table)
  colnames(subcat_df) <- c("Subcategory", "Count")
  
  # Calculate percentage labels
  subcat_df$Perc <- round(subcat_df$Count / sum(subcat_df$Count) * 100, 1)
  
  # Plot pie chart
  ggplot(subcat_df, aes(x = "", y = Count, fill = Subcategory)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    theme_void() +
    labs(title = "Subcategories of Intergenic and Intronic Peaks") +
    geom_text(aes(label = Perc), position = position_stack(vjust = 0.5), size = 3)
  
  ggsave(paste(dir,'/',name,"_intergenic_intron.pdf",sep=""),width=5,height=5,dpi=300)
  
  ## save the object
  saveRDS(peaks,paste(dir,'/',name,'.rds',sep=""))
  return(peaks)
}

peak_go <- function (name,dir)
{
  library(GenomicFeatures)
  library(GenomicRanges)
  library(rtracklayer)
  library(ChIPseeker)
  library(ggplot2)
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  
  peaks<-readRDS(paste(dir,paste(name,".rds",sep=""),sep="/"))
  peaks_near_genes <- peaks[!is.na(mcols(peaks)$distance_to_gene) &
                              mcols(peaks)$distance_to_gene < 500]
  # Get unique gene IDs from your nearest_gene column
  gene_symbols <- unique(na.omit(mcols(peaks_near_genes)$nearest_gene))
  # print(gene_symbols[1])
  gene_ids_clean <- sub("\\..*$", "", gene_symbols)
  # Map gene symbols to Entrez IDs for enrichment
  entrez_ids <- mapIds(org.Hs.eg.db,
                       keys = gene_ids_clean,
                       column = "ENTREZID",
                       keytype = "ENSEMBL",
                       multiVals = "first")
  entrez_ids <- na.omit(unique(entrez_ids))
  ego_mf <- enrichGO(gene = entrez_ids,
                     OrgDb   = org.Hs.eg.db,
                     keyType = "ENTREZID",
                     ont     = "MF",   # can also use "MF" or "CC" or "BP"
                     pAdjustMethod= "BH",
                     qvalueCutoff = 0.05,
                     readable     = TRUE)
  
  ego_bp <- enrichGO(gene         = entrez_ids,
                     OrgDb        = org.Hs.eg.db,
                     keyType      = "ENTREZID",
                     ont          = "BP",   # can also use "MF" or "CC" or "BP"
                     pAdjustMethod= "BH",
                     qvalueCutoff = 0.05,
                     readable     = TRUE)
  
  ego_cc <- enrichGO(gene         = entrez_ids,
                     OrgDb        = org.Hs.eg.db,
                     keyType      = "ENTREZID",
                     ont          = "CC",   # can also use "MF" or "CC" or "BP"
                     pAdjustMethod= "BH",
                     qvalueCutoff = 0.05,
                     readable     = TRUE)
  
  # head(ego)
  
  ekegg <- enrichKEGG(gene         = entrez_ids,
                      organism     = "hsa",  # human
                      pvalueCutoff = 0.05)
  
  # Convert Entrez IDs back to gene symbols in results
  ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  
  library(ReactomePA)
  
  ereact <- enrichPathway(gene         = entrez_ids,
                          organism     = "human",
                          pvalueCutoff = 0.05,
                          readable     = TRUE)
  
  library(enrichplot)
  
  # Dotplot for GO
  dotplot(ego_bp, showCategory=20) + ggtitle("GO Biological Process")
  ggsave(paste(dir,'/',name,"_enrichgo_bp_dotplot.pdf",sep=""),width=5,height=7,dpi=300)
  dotplot(ego_mf, showCategory=20) + ggtitle("GO Mol. Function")
  ggsave(paste(dir,'/',name,"_enrichgo_mf_dotplot.pdf",sep=""),width=5,height=7,dpi=300)
  dotplot(ego_cc, showCategory=20) + ggtitle("GO Cell. compartment")
  ggsave(paste(dir,'/',name,"_enrichgo_cc_dotplot.pdf",sep=""),width=5,height=7,dpi=300)
  
  # Barplot for KEGG
  barplot(ekegg, showCategory=15) + ggtitle("KEGG Pathways")
  ggsave(paste(dir,'/',name,"_enrichkegg_barplot.pdf",sep=""),width=5,height=5,dpi=300)

  return(list(ego_bp,ego_mf,ego_cc,ekegg,ereact))
}
