# Auto-generated R packages installation script for profile: Minimal
# Generated on: 西元2025年07月12日 (週六) 08時04分16秒 CST

# 設定 CRAN mirror
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# 安裝套件
packages <- c("base", "utils")

for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
    }
}

cat("✅ R packages installation completed for profile: Minimal\n")
