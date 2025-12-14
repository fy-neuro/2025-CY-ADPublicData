if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_path <- rstudioapi::getSourceEditorContext()$path
  script_dir <- dirname(script_path)
  cat("Script path (RStudio):", script_path, "\n")
  cat("Script directory:", script_dir, "\n")
}

# Method 2: Using here package for project root
project_root <- here::here()
cat("Project root (here):", project_root, "\n")
