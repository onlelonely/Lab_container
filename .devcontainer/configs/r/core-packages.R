#!/usr/bin/env Rscript

# Core R packages installation script
cat("Installing core R packages...\n")

# Base packages
core_packages <- c(
  "tidyverse",
  "data.table",
  "ggplot2",
  "dplyr",
  "readr",
  "stringr",
  "lubridate",
  "devtools",
  "roxygen2",
  "testthat",
  "lintr",
  "rmarkdown",
  "knitr"
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