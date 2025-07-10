#!/usr/bin/env Rscript

# Bioconductor packages installation script
cat("Installing Bioconductor packages...\n")

# Install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Core Bioconductor packages (excluding large genome packages)
bioc_packages <- c(
  # Core infrastructure
  "Biostrings",
  "GenomicRanges", 
  "IRanges",
  "GenomeInfoDb",
  "biomaRt",
  "AnnotationDbi",
  
  # Statistical analysis
  "DESeq2",
  "edgeR", 
  "limma",
  
  # Annotation databases (smaller ones)
  "GO.db",
  "KEGGREST",
  "org.Hs.eg.db",
  
  # File I/O and manipulation
  "VariantAnnotation",
  "Rsamtools",
  "rtracklayer",
  
  # Analysis tools
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