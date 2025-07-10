#!/usr/bin/env Rscript

# Bioconductor packages installation script
cat("Installing Bioconductor packages...\n")

# Install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Bioconductor packages
bioc_packages <- c(
  "Biostrings",
  "GenomicRanges",
  "IRanges",
  "GenomeInfoDb",
  "biomaRt",
  "AnnotationDbi",
  "DESeq2",
  "edgeR",
  "limma",
  "GOstats",
  "GO.db",
  "KEGGREST",
  "reactome.db",
  "org.Hs.eg.db",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "BSgenome.Hsapiens.UCSC.hg38",
  "VariantAnnotation",
  "Rsamtools",
  "rtracklayer",
  "ChIPseeker",
  "clusterProfiler"
)

# Install packages with error handling
for (pkg in bioc_packages) {
  tryCatch({
    if (!require(pkg, character.only = TRUE)) {
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
      cat(paste("Successfully installed:", pkg, "\n"))
    } else {
      cat(paste("Already installed:", pkg, "\n"))
    }
  }, error = function(e) {
    cat(paste("Error installing", pkg, ":", e$message, "\n"))
  })
}

cat("Bioconductor packages installation completed.\n")