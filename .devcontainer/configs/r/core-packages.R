#!/usr/bin/env Rscript

# Core R packages installation script
cat("Installing core R packages...\n")

# Set options for non-interactive installation
options(repos = c(CRAN = "https://cloud.r-project.org/"))
options(timeout = 300)  # 5 minute timeout

# Base packages (minimal set to reduce failure risk)
core_packages <- c(
  "data.table",     # High-performance data manipulation
  "devtools",       # Development tools
  "rmarkdown",      # R Markdown documents
  "knitr"           # Dynamic report generation
)

# Install packages with robust error handling
success_count <- 0
total_packages <- length(core_packages)

for (pkg in core_packages) {
  cat(paste("Installing package:", pkg, "\n"))
  
  tryCatch({
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE, 
                      repos = "https://cloud.r-project.org/",
                      quiet = FALSE)
      
      # Verify installation
      if (require(pkg, character.only = TRUE, quietly = TRUE)) {
        cat(paste("✓ Successfully installed:", pkg, "\n"))
        success_count <- success_count + 1
      } else {
        cat(paste("✗ Installation verification failed for:", pkg, "\n"))
      }
    } else {
      cat(paste("✓ Already installed:", pkg, "\n"))
      success_count <- success_count + 1
    }
  }, error = function(e) {
    cat(paste("✗ Error installing", pkg, ":", e$message, "\n"))
  })
}

cat(paste("R packages installation completed:", success_count, "of", total_packages, "packages successful\n"))

# Install tidyverse separately as it's large and can fail
cat("Installing tidyverse (may take longer)...\n")
tryCatch({
  if (!require("tidyverse", character.only = TRUE, quietly = TRUE)) {
    install.packages("tidyverse", dependencies = TRUE, 
                    repos = "https://cloud.r-project.org/")
    cat("✓ Tidyverse installed successfully\n")
  } else {
    cat("✓ Tidyverse already installed\n")
  }
}, error = function(e) {
  cat(paste("✗ Tidyverse installation failed:", e$message, "\n"))
  cat("This is not critical - continuing with core functionality\n")
})

cat("Core R packages installation process completed.\n")