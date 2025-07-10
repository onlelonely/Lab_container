#!/usr/bin/env Rscript

# Extended R packages installation script
cat("Installing extended R packages...\n")

# Extended packages
extended_packages <- c(
  "plotly",
  "DT",
  "shiny",
  "shinydashboard",
  "shinyWidgets",
  "htmlwidgets",
  "leaflet",
  "RColorBrewer",
  "viridis",
  "ggsci",
  "ggthemes",
  "patchwork",
  "cowplot",
  "scales",
  "reshape2",
  "tidyr",
  "broom",
  "modelr",
  "caret",
  "randomForest",
  "glmnet",
  "survival",
  "survminer",
  "pheatmap",
  "ComplexHeatmap",
  "circlize",
  "VennDiagram",
  "corrplot",
  "psych",
  "Hmisc",
  "doParallel",
  "foreach",
  "future",
  "furrr"
)

# Install packages with error handling
for (pkg in extended_packages) {
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

cat("Extended R packages installation completed.\n")