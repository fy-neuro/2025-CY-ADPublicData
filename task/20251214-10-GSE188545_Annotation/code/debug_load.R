cat("Debug loading RDS...\n")
file_path <- "GSE188545_sobj_annotated_fast.rds"
cat("File exists:", file.exists(file_path), "\n")
cat("File size:", file.size(file_path), "bytes\n")
cat("Attempting to read...\n")
tryCatch({
  sobj <- readRDS(file_path)
  cat("Success!\n")
  cat("Class:", class(sobj), "\n")
  cat("Cells:", ncol(sobj), "\n")
}, error = function(e) {
  cat("Error:", e$message, "\n")
  # Try reading raw bytes
  cat("Attempting to read raw bytes...\n")
  con <- file(file_path, "rb")
  header <- readBin(con, "raw", n = 10)
  close(con)
  cat("First 10 bytes:", paste(header, collapse = " "), "\n")
})