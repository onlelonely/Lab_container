# CI-friendly R packages (fast installation, no compilation)
# Generated for reliable CI builds

cat("Installing CI-friendly R packages...\n")

# Essential packages from CRAN (binary packages, fast install)
basic_packages <- c(
  "ggplot2",     # Visualization
  "dplyr",       # Data manipulation  
  "readr",       # Data reading
  "tidyr",       # Data tidying
  "stringr"      # String manipulation
)

# Install packages with error handling
for (pkg in basic_packages) {
  tryCatch({
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org/", 
                      dependencies = FALSE, quiet = TRUE)
      cat(paste("✓ Installed", pkg, "\n"))
    } else {
      cat(paste("✓", pkg, "already installed\n"))
    }
  }, error = function(e) {
    cat(paste("⚠ Failed to install", pkg, ":", e$message, "\n"))
  })
}

# Packages excluded from CI (too slow/complex):
# - Bioconductor packages (BiocManager, Biostrings, GenomicRanges)
# - Packages requiring compilation
# - Large statistical packages
cat("CI R package installation complete\n")