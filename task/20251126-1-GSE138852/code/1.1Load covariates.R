# ...existing code...

# Load the covariates file as a data.frame
cat("Loading covariates from data/GSE138852_covariates.csv.gz...\n")
covariates <- read.csv(here("data", "GSE138852_covariates.csv.gz"), 
                       row.names = 1, 
                       header = TRUE, 
                       check.names = FALSE)

# Display basic information about the covariates
cat("Covariates summary:\n")
cat("Number of samples:", nrow(covariates), "\n")
cat("Number of variables:", ncol(covariates), "\n")
cat("Column names:", paste(colnames(covariates), collapse = ", "), "\n")

# Show first few rows
cat("First few rows:\n")
print(head(covariates))

# ...existing code...