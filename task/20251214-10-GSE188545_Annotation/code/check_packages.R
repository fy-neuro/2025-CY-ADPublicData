# Check and install required packages for annotation pipeline
cat("Checking required packages...\n")

required_packages <- c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "patchwork",
  "here",
  "viridis",
  "plyr",
  "stringr",
  "purrr",
  "clusterProfiler",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "ggrepel",
  "cols4all",
  "ggpubr"
)

# Function to check and install packages
check_and_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste("Installing", pkg, "...\n"))
    tryCatch({
      install.packages(pkg, repos = "https://cloud.r-project.org")
      cat(paste("  Successfully installed", pkg, "\n"))
    }, error = function(e) {
      cat(paste("  Failed to install", pkg, ":", e$message, "\n"))
    })
  } else {
    cat(paste("  ", pkg, "already installed\n"))
  }
}

# Check each package
for (pkg in required_packages) {
  check_and_install(pkg)
}

# Special case for Bioconductor packages
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  cat("Installing org.Hs.eg.db from Bioconductor...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("org.Hs.eg.db")
}

if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  cat("Installing clusterProfiler from Bioconductor...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("clusterProfiler")
}

cat("\nAll packages checked.\n")
cat("To load all packages, run:\n")
cat("lapply(required_packages, library, character.only = TRUE)\n")