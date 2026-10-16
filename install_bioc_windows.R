# Run this script in R after activating the Conda environment.
# This installs Bioconductor packages from R/Bioconductor rather than trying
# to install Linux-specific Conda builds on Windows.

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Use the Bioconductor release appropriate for your installed R version.
# Do not force an old release unless you specifically need it.
bioc_pkgs <- c(
  "AnnotationDbi", "AnnotationFilter", "apeglm", "BindingSiteFinder",
  "Biobase", "BiocFileCache", "BiocGenerics", "BiocIO",
  "BiocParallel", "BiocStyle", "biomaRt", "Biostrings", "biovizBase",
  "BSgenome", "BSgenome.Hsapiens.UCSC.hg19", "BSgenome.Hsapiens.UCSC.hg38",
  "clusterProfiler", "ComplexHeatmap", "DelayedArray", "DESeq2", "DOSE",
  "EnrichedHeatmap", "enrichplot", "ensembldb", "fgsea", "GenomeInfoDb",
  "GenomeInfoDbData", "GenomicAlignments", "GenomicFeatures", "GenomicRanges",
  "ggtree", "GO.db", "GOSemSim", "graph", "graphite", "Guitar", "Gviz",
  "HDO.db", "IRanges", "KEGGREST", "matrixGenerics", "msa", "org.Hs.eg.db",
  "plyranges", "pqsfinder", "ProtGenerics", "qvalue", "Reactome.db",
  "ReactomePA", "Rhtslib", "Rsamtools", "rtracklayer", "S4Arrays",
  "S4Vectors", "SparseArray", "SummarizedExperiment", "treeio",
  "TxDb.Hsapiens.UCSC.hg38.knownGene", "VariantAnnotation", "wavClusteR",
  "XVector", "zlibbioc"
)

BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)
