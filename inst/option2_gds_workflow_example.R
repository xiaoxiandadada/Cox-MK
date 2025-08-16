# CoxMK Option 2: Step-by-step workflow with GDS file reuse
# This example demonstrates how to:
# 1. Generate knockoffs and save to GDS file
# 2. Fit null model separately
# 3. Load knockoffs from GDS file for analysis

library(CoxMK)

# Load example data
extdata_path <- system.file("extdata", package = "CoxMK")
plink_prefix <- file.path(extdata_path, "sample")
pheno_data <- prepare_phenotype(file.path(extdata_path, "tte_phenotype.txt"))
covar_data <- load_covariates(file.path(extdata_path, "covariates.txt"))

cat("=== CoxMK Option 2: GDS Workflow Example ===\n\n")

# Step 1: Generate knockoffs and save to GDS file
cat("Step 1: Generating knockoffs and saving to GDS file...\n")
plink_data <- load_plink_data(plink_prefix)
knockoffs <- create_knockoffs(
  X = plink_data$genotypes,
  pos = plink_data$positions,
  M = 3
)
gds_file <- knockoffs$gds_file  # GDS file path
cat("✓ Knockoffs saved to:", basename(gds_file), "\n\n")

# Step 2: Fit null model separately (can be done independently)
cat("Step 2: Fitting null Cox model...\n")
null_model <- fit_null_cox_model(
  time = pheno_data$time,
  status = pheno_data$status, 
  covariates = covar_data
)
cat("✓ Null model fitted successfully\n\n")

# Step 3: Load knockoffs from GDS file and run analysis
cat("Step 3: Running analysis with pre-generated knockoffs...\n")
result <- cox_knockoff_analysis(
  plink_prefix = plink_prefix,
  time = pheno_data$time,
  status = pheno_data$status,
  covariates = covar_data,
  null_model = null_model,
  gds_file = gds_file,  # Use pre-generated GDS file
  M = 3,  # This should match the M used in step 1
  fdr = 0.05
)

# Display results
cat("\n=== Final Results ===\n")
cat("Total variables tested:", result$summary$n_variables, "\n")
cat("Variables selected:", result$summary$n_selected, "\n")
cat("Selection rate:", round(result$summary$selection_rate * 100, 2), "%\n")
cat("FDR level:", result$fdr, "\n")
cat("GDS file used:", basename(result$gds_file), "\n")

cat("\n=== Benefits of Option 2 Workflow ===\n")
cat("✓ Knockoffs can be generated once and reused multiple times\n")
cat("✓ Different null models can be tested with same knockoffs\n")
cat("✓ Different FDR levels can be tested without regenerating knockoffs\n")
cat("✓ More efficient for parameter exploration and sensitivity analysis\n")
cat("✓ GDS files can be shared across different analysis sessions\n")

cat("\nWorkflow completed successfully!\n")
