#!/usr/bin/env Rscript

# Core R packages installation script
cat("Installing core R packages...\n")

# Base packages (avoiding tidyverse duplicates)
core_packages <- c(
  "tidyverse",      # Includes ggplot2, dplyr, readr, stringr, lubridate, etc.
  "data.table",     # High-performance data manipulation
  "devtools",       # Development tools
  "roxygen2",       # Documentation
  "testthat",       # Testing framework
  "lintr",          # Code linting
  "rmarkdown",      # R Markdown documents
  "knitr"           # Dynamic report generation
)

# Install packages with error handling
for (pkg in core_packages) {
  tryCatch({
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org/")
      cat(paste("Successfully installed:", pkg, "\n"))
    } else {
      cat(paste("Already installed:", pkg, "\n"))
    }
  }, error = function(e) {
    cat(paste("Error installing", pkg, ":", e$message, "\n"))
  })
}

cat("Core R packages installation completed.\n")