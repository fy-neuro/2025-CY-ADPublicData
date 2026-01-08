# Load required libraries
library(Seurat)   # For single-cell RNA-seq analysis
library(loupeR)   # For creating Loupe Browser files


# Function to prepare Seurat object for loupeR
# Handles Seurat v5 layers and cluster column detection
# 准备 Seurat 对象用于 loupeR 导出
# 处理 Seurat v5 多层数据并自动检测聚类列
prepare_seurat_for_loupe <- function(sobj) {
  cat("  Preparing Seurat object for loupeR...\n")
  
  # Check if it's Seurat v5 with multiple layers
  assay_name <- DefaultAssay(sobj)
  assay_obj <- sobj[[assay_name]]
  
  # Check for multiple layers (Seurat v5 feature)
  if (inherits(assay_obj, "Assay5") && length(Layers(assay_obj, search = "counts")) > 1) {
    cat("  Object has multiple layers. Joining layers...\n")
    tryCatch({
      sobj <- JoinLayers(sobj, assay = assay_name)
      cat("  Layers joined successfully.\n")
    }, error = function(e) {
      cat("  Warning: Failed to join layers:", e$message, "\n")
    })
  }
  
  # Ensure there's a cluster column for loupeR
  # loupeR looks for 'active_cluster' by default, but we'll specify our own
  meta_cols <- colnames(sobj@meta.data)
  
  # Prefer seurat_clusters, then celltype, then any cluster column
  preferred_cluster <- NULL
  if ("seurat_clusters" %in% meta_cols) {
    preferred_cluster <- "seurat_clusters"
  } else if ("celltype" %in% meta_cols) {
    preferred_cluster <- "celltype"
  } else {
    # Look for any column with 'cluster' in name
    cluster_cols <- grep("cluster", meta_cols, ignore.case = TRUE, value = TRUE)
    if (length(cluster_cols) > 0) {
      preferred_cluster <- cluster_cols[1]
    }
  }
  
  cat("  Preferred cluster column:", ifelse(is.null(preferred_cluster), "None found", preferred_cluster), "\n")
  
  return(list(sobj = sobj, cluster_col = preferred_cluster))
}


# Main function to create Loupe Browser (.cloupe) file from Seurat object
# 从 Seurat 对象创建 Loupe Browser (.cloupe) 文件的主函数
# Parameters:
#   sobj: Seurat object OR path to .rds file containing Seurat object
#   output_dir: Directory where .cloupe file will be saved (default: "data/loupe_files")
#   output_name: Base name for output file (without extension) (default: inferred from sobj)
#   force: Overwrite existing file (default: TRUE)
#   verbose: Print progress messages (default: TRUE)
# Returns:
#   logical TRUE if successful, FALSE otherwise
create_loupe_file <- function(sobj,
                              output_dir = "data/loupe_files",
                              output_name = NULL,
                              force = TRUE,
                              verbose = TRUE) {
  
  # Start timing
  start_time <- Sys.time()
  
  # Helper function for verbose output
  vcat <- function(...) {
    if (verbose) cat(...)
  }
  
  vcat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
  vcat("Creating Loupe Browser file from Seurat object\n")
  
  # Check if sobj is a path (character string) or already a Seurat object
  # 检查输入是文件路径还是 Seurat 对象
  if (is.character(sobj) && length(sobj) == 1) {
    # Assume it's a file path
    input_path <- sobj
    vcat("Input path:", input_path, "\n")
    
    # Check if file exists
    if (!file.exists(input_path)) {
      vcat("ERROR: File not found.\n")
      return(FALSE)
    }
    
    # Load the Seurat object
    vcat("Loading Seurat object...\n")
    tryCatch({
      sobj <- readRDS(input_path)
      vcat("Object loaded successfully.\n")
    }, error = function(e) {
      vcat("ERROR loading RDS file:", e$message, "\n")
      return(FALSE)
    })
    
    # Set default output name based on file name if not provided
    if (is.null(output_name)) {
      output_name <- tools::file_path_sans_ext(basename(input_path))
      vcat("Using output name from file:", output_name, "\n")
    }
  } else {
    # sobj is already a Seurat object
    if (!inherits(sobj, "Seurat")) {
      vcat("ERROR: Input is not a Seurat object.\n")
      return(FALSE)
    }
    vcat("Input is a Seurat object.\n")
    
    # Set default output name if not provided
    if (is.null(output_name)) {
      output_name <- deparse(substitute(sobj))
      vcat("Using output name from object variable:", output_name, "\n")
    }
  }
  
  # Print basic info
  vcat("Number of cells:", ncol(sobj), "\n")
  vcat("Number of features:", nrow(sobj), "\n")
  vcat("Assays available:", paste(names(sobj@assays), collapse = ", "), "\n")
  
  # Prepare object for loupeR
  prep_result <- prepare_seurat_for_loupe(sobj)
  sobj <- prep_result$sobj
  cluster_col <- prep_result$cluster_col
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    vcat("Created output directory:", output_dir, "\n")
  }
  
  # Try to create loupe file using create_loupe_from_seurat (preferred)
  # 尝试使用 create_loupe_from_seurat 创建 loupe 文件（首选）
  vcat("Creating loupe file...\n")
  
  if (!is.null(cluster_col)) {
    vcat("Using cluster column:", cluster_col, "\n")
    
    # Try with explicit metadata_cols parameter
    tryCatch({
      result <- create_loupe_from_seurat(
        obj = sobj,
        output_dir = output_dir,
        output_name = output_name,
        metadata_cols = cluster_col,  # Explicitly specify cluster column
        force = force
      )
      vcat("SUCCESS: Loupe file created using cluster column", cluster_col, "\n")
      
      # Print summary
      end_time <- Sys.time()
      vcat("Time elapsed:", format(end_time - start_time, digits = 3), "\n")
      return(TRUE)
    }, error = function(e) {
      vcat("Failed with explicit cluster column:", e$message, "\n")
      vcat("Trying without specifying cluster column...\n")
    })
  }
  
  # Try without specifying cluster column
  tryCatch({
    result <- create_loupe_from_seurat(
      obj = sobj,
      output_dir = output_dir,
      output_name = output_name,
      force = force
    )
    vcat("SUCCESS: Loupe file created with default parameters.\n")
    
    end_time <- Sys.time()
    vcat("Time elapsed:", format(end_time - start_time, digits = 3), "\n")
    return(TRUE)
  }, error = function(e) {
    vcat("Failed with default parameters:", e$message, "\n")
    
    # Last resort: manual create_loupe
    # 最后手段：手动调用 create_loupe
    vcat("Trying manual create_loupe as last resort...\n")
    tryCatch({
      # Get default assay
      default_assay <- DefaultAssay(sobj)
      
      # Get count matrix
      count_mat <- GetAssayData(sobj, assay = default_assay, layer = "counts")
      if (is.null(count_mat) || all(dim(count_mat) == 0)) {
        count_mat <- GetAssayData(sobj, assay = default_assay, layer = "data")
        vcat("  Using data layer instead of counts\n")
      }
      
      vcat("  Count matrix dimensions:", dim(count_mat), "\n")
      
      # Get clusters
      if (!is.null(cluster_col)) {
        clusters <- sobj@meta.data[[cluster_col]]
      } else {
        clusters <- Idents(sobj)
      }
      
      # Convert to list format expected by create_loupe
      clusters_list <- list(cluster = clusters)
      
      # Get projections
      projections_list <- list()
      reduc_names <- names(sobj@reductions)
      
      for (reduc_name in reduc_names) {
        emb <- sobj@reductions[[reduc_name]]@cell.embeddings
        if (ncol(emb) >= 2) {
          projections_list[[reduc_name]] <- emb[, 1:2]
        }
      }
      
      if (length(projections_list) == 0) {
        vcat("  Warning: No projections found. Creating dummy projection.\n")
        set.seed(123)
        dummy_proj <- matrix(rnorm(ncol(sobj) * 2), ncol = 2)
        rownames(dummy_proj) <- colnames(sobj)
        projections_list[["dummy"]] <- dummy_proj
      }
      
      # Get feature IDs
      feature_ids <- rownames(sobj)
      
      # Create loupe file
      result <- create_loupe(
        count_mat = count_mat,
        clusters = clusters_list,
        projections = projections_list,
        output_dir = output_dir,
        output_name = output_name,
        feature_ids = feature_ids,
        force = force
      )
      
      vcat("SUCCESS: Loupe file created manually.\n")
      
      end_time <- Sys.time()
      vcat("Time elapsed:", format(end_time - start_time, digits = 3), "\n")
      return(TRUE)
      
    }, error = function(e2) {
      vcat("Manual create_loupe also failed:", e2$message, "\n")
      return(FALSE)
    })
  })
}

# Usage examples 使用示例
# 
# 1. From a Seurat object already loaded:
#    create_loupe_file(sobj = my_seurat_object,
#                      output_dir = "data/loupe_files",
#                      output_name = "my_dataset")
#
# 2. From a .rds file path:
#    create_loupe_file(sobj = "path/to/annotated_seurat.rds",
#                      output_dir = "results/loupe",
#                      output_name = "annotated_dataset")
#
# 3. Minimal usage with defaults:
#    create_loupe_file("data/processed/sobj.rds")
#
# Note: Requires loupeR and Seurat packages installed.
# 注意：需要先安装 loupeR 和 Seurat 包。