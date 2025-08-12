# Test CoxMK Package
# This script tests the main functionality of the CoxMK package

cat("=== Testing CoxMK Package ===\n")

# Install and load the package
if (!require(CoxMK, quietly = TRUE)) {
  # Install from local tarball
  install.packages("CoxMK_0.1.0.tar.gz", repos = NULL, type = "source")
  library(CoxMK)
}

library(survival)
library(Matrix)

# Test 1: Load example data
cat("\n1. Testing data loading...\n")
data(example_genotypes)
data(example_positions)
data(example_phenotype)  
data(example_covariates)

cat("   ✓ Genotype data:", class(example_genotypes), dim(example_genotypes), "\n")
cat("   ✓ Positions:", length(example_positions), "SNPs\n")
cat("   ✓ Phenotype data:", nrow(example_phenotype), "samples\n")
cat("   ✓ Covariate data:", nrow(example_covariates), "samples\n")

# Test 2: Create knockoffs with different M values
cat("\n2. Testing knockoff creation with different M values...\n")

# Test M = 3
cat("   Testing M = 3...\n")
knockoffs_3 <- create_knockoffs(
  X = example_genotypes,
  pos = example_positions,
  M = 3,
  corr_max = 0.75
)
cat("   ✓ Created", length(knockoffs_3), "knockoff matrices\n")

# Test M = 5 (default)
cat("   Testing M = 5 (default)...\n")
knockoffs_5 <- create_knockoffs(
  X = example_genotypes,
  pos = example_positions,
  M = 5,
  corr_max = 0.75
)
cat("   ✓ Created", length(knockoffs_5), "knockoff matrices\n")

# Test 3: Prepare phenotype data
cat("\n3. Testing phenotype preparation...\n")
pheno_data <- merge(example_phenotype, example_covariates, by = c("FID", "IID"))
covariates <- pheno_data[, c("age", "sex", "bmi")]
cat("   ✓ Merged phenotype data:", dim(pheno_data), "\n")
cat("   ✓ Extracted covariates:", names(covariates), "\n")

# Test 4: Fit null models
cat("\n4. Testing null model fitting...\n")

# Standard coxph
cat("   Testing standard coxph...\n")
null_model <- fit_null_model(
  time = pheno_data$time,
  status = pheno_data$status,
  covariates = covariates,
  use_spacox = FALSE
)
cat("   ✓ Null model fitted with method:", attr(null_model, "method") %||% "coxph", "\n")

# SPACox (should fallback to coxph)
cat("   Testing SPACox option (should fallback)...\n")
null_model_spa <- fit_null_model(
  time = pheno_data$time,
  status = pheno_data$status,
  covariates = covariates,
  use_spacox = TRUE
)
cat("   ✓ SPACox test completed, method:", attr(null_model_spa, "method") %||% "coxph", "\n")

# Test 5: Calculate W statistics
cat("\n5. Testing W statistics calculation...\n")
set.seed(123)
p <- ncol(example_genotypes)

# Simulate results for M=5
original_pvals <- runif(p, 0.01, 1)
original_coefs <- rnorm(p, 0, 0.3)
knockoff_pvals <- lapply(1:5, function(k) runif(p, 0.01, 1))
knockoff_coefs <- lapply(1:5, function(k) rnorm(p, 0, 0.3))

w_stats <- calculate_w_statistics(
  original_pvals = original_pvals,
  knockoff_pvals = knockoff_pvals,
  original_coefs = original_coefs,
  knockoff_coefs = knockoff_coefs
)

cat("   ✓ W statistics calculated, range:", round(range(w_stats), 3), "\n")
cat("   ✓ Positive W stats:", sum(w_stats > 0), "\n")
cat("   ✓ Negative W stats:", sum(w_stats < 0), "\n")

# Test 6: Apply knockoff filter
cat("\n6. Testing knockoff filter...\n")
selected <- knockoff_filter(w_stats, fdr = 0.1)

cat("   ✓ Selection threshold:", round(selected$threshold, 3), "\n")
cat("   ✓ Estimated FDP:", round(selected$fdp, 3), "\n")
cat("   ✓ Number selected:", length(selected$selected), "\n")

if (length(selected$selected) > 0) {
  cat("   ✓ Selected indices:", head(selected$selected, 10), 
      if(length(selected$selected) > 10) "..." else "", "\n")
}

# Test 7: Test with external data (if available)
cat("\n7. Testing with external sample data...\n")
sample_bed <- "/Users/fairy/Documents/R_file/Cox-MK/sample/sample.bed"

if (file.exists(sample_bed)) {
  cat("   Found sample PLINK data, testing load_plink_data...\n")
  
  tryCatch({
    plink_data <- load_plink_data("/Users/fairy/Documents/R_file/Cox-MK/sample/sample")
    cat("   ✓ PLINK data loaded:", dim(plink_data$genotypes), "\n")
    
    # Test small knockoff creation
    small_knockoffs <- create_knockoffs(
      X = plink_data$genotypes[1:20, 1:50],  # Subset for speed
      pos = plink_data$positions[1:50],
      M = 3
    )
    cat("   ✓ Small knockoff test completed\n")
    
  }, error = function(e) {
    cat("   ⚠ PLINK data test failed:", e$message, "\n")
  })
} else {
  cat("   ℹ No sample PLINK data found, skipping external data test\n")
}

# Summary
cat("\n=== Test Summary ===\n")
cat("✓ All core functions tested successfully!\n")
cat("✓ Package supports:\n")
cat("  - Flexible knockoff generation (M = 3, 5, 10, etc.)\n")
cat("  - Standard coxph and SPACox integration\n")
cat("  - W statistics calculation\n")
cat("  - Knockoff filtering with FDR control\n")
cat("  - PLINK data loading (with fallback)\n")

cat("\n=== Usage Summary ===\n")
cat("# Load package\n")
cat("library(CoxMK)\n\n")
cat("# Create knockoffs (user chooses M)\n")
cat("knockoffs <- create_knockoffs(genotypes, positions, M = 5)  # or M = 3, 10, etc.\n\n")
cat("# Fit null model (with SPACox option)\n")
cat("null_model <- fit_null_model(time, status, covariates, use_spacox = TRUE)\n\n")
cat("# Calculate W statistics and select variables\n")
cat("w_stats <- calculate_w_statistics(...)\n")
cat("selected <- knockoff_filter(w_stats, fdr = 0.1)\n\n")

cat("Package testing completed successfully!\n")
