# CoxKnockoff Package Example
# This script demonstrates the complete workflow using real PLINK data

# Load required packages
library(CoxMK)
library(survival)
library(Matrix)

cat("=== CoxMK Package Example with Real PLINK Data ===\n\n")

# 1. Load real PLINK data from extdata
cat("1. Loading real PLINK data...\n")
extdata_path <- system.file("extdata", package = "CoxMK")

# Load PLINK files
plink_data <- load_plink_data(file.path(extdata_path, "sample"))
cat("   - PLINK data loaded successfully\n")
cat("   - Genotypes:", dim(plink_data$genotypes), "(samples x SNPs)\n")
cat("   - Positions:", length(plink_data$positions), "SNPs\n")

# Load phenotype data
phenotype_data <- read.table(file.path(extdata_path, "tte_phenotype.txt"), 
                            header = TRUE, stringsAsFactors = FALSE)
cat("   - Phenotypes:", nrow(phenotype_data), "samples\n")

# Load covariate data  
covariate_data <- read.table(file.path(extdata_path, "covariates.txt"), 
                            header = TRUE, stringsAsFactors = FALSE)
cat("   - Covariates:", nrow(covariate_data), "samples\n\n")

# 2. Create knockoff variables using real data
cat("2. Creating knockoff variables...\n")
knockoffs <- create_knockoffs(
  X = plink_data$genotypes,
  pos = plink_data$positions,
  M = 5,  # Default 5 knockoffs, users can change this (e.g., M = 3, 10, etc.)
  corr_max = 0.75,
  maxBP.neighbor = 1e5
)

cat("   - Created", length(knockoffs$knockoffs), "knockoff matrices (M=5 by default)\n")
cat("   - Users can specify different M values (e.g., M=3, M=10)\n")
cat("   - Each matrix dimensions:", dim(knockoffs$knockoffs[[1]]), "\n")
if (!is.null(knockoffs$gds_file)) {
  cat("   - Knockoffs saved to GDS file:", basename(knockoffs$gds_file), "\n")
}
cat("\n")

# 3. Prepare phenotype data
cat("3. Preparing phenotype data...\n")
pheno_data <- merge(phenotype_data, covariate_data, by = c("FID", "IID"))
covariates <- pheno_data[, c("age", "sex", "bmi")]

cat("   - Merged data dimensions:", dim(pheno_data), "\n")
cat("   - Covariates used:", names(covariates), "\n")
cat("   - Time to event range:", range(pheno_data$time), "days\n")
cat("   - Events observed:", sum(pheno_data$status), "out of", nrow(pheno_data), "samples\n\n")

# 4. Fit Cox model and perform association analysis
cat("4. Fitting Cox model and performing association analysis...\n")
cat("   - Using Cox regression with survival analysis\n")
cat("   - For large-scale studies, consider SPACox: https://github.com/WenjianBI/SPACox\n")

# Perform Cox regression analysis
cox_results <- fit_cox_spa(
  X = plink_data$genotypes,
  time = pheno_data$time,
  status = pheno_data$status,
  covariates = covariates,
  use_spa = FALSE  # Set to TRUE if SPAtest package is available
)

cat("   - Cox analysis completed using", cox_results$method, "\n")
cat("   - Number of events:", sum(pheno_data$status), "\n")
cat("   - Median survival time:", median(pheno_data$time), "days\n")
cat("   - P-value range:", range(cox_results$p_values), "\n\n")

# 5. Use real Cox regression results instead of simulation
cat("5. Using real Cox regression results...\n")
p <- ncol(plink_data$genotypes)

# Use actual Cox regression results
original_pvals <- cox_results$p_values
original_coefs <- cox_results$test_stats  # Using test statistics as coefficients

# Perform knockoff Cox analysis  
cat("   - Performing knockoff analysis on", length(knockoffs$knockoffs), "knockoff matrices...\n")
knockoff_results <- lapply(knockoffs$knockoffs, function(ko_matrix) {
  fit_cox_spa(
    X = ko_matrix,
    time = pheno_data$time,
    status = pheno_data$status,
    covariates = covariates,
    use_spa = FALSE
  )
})

# Extract p-values and coefficients
knockoff_pvals <- lapply(knockoff_results, function(x) x$p_values)
knockoff_coefs <- lapply(knockoff_results, function(x) x$test_stats)

cat("   - Analysis completed for", p, "SNPs\n")
cat("   - Original p-value range:", range(original_pvals), "\n")
cat("   - Significant SNPs (p < 0.05):", sum(original_pvals < 0.05), "\n\n")

# 6. Calculate W statistics
cat("6. Calculating W statistics...\n")
w_stats <- calculate_w_statistics(
  t_orig = cox_results$test_stats,
  t_knock = lapply(knockoff_results, function(x) x$test_stats),
  method = "difference"
)

cat("   - W statistics range:", round(range(w_stats), 3), "\n")
cat("   - Positive W stats:", sum(w_stats > 0), "\n")
cat("   - Negative W stats:", sum(w_stats < 0), "\n")
cat("   - Significant original SNPs (p < 0.05):", sum(original_pvals < 0.05), "\n\n")

# 7. Apply knockoff filter
cat("7. Applying knockoff filter...\n")
selected <- knockoff_filter(w_stats, fdr = 0.1)

cat("   - Selection threshold:", round(attr(selected, "threshold"), 3), "\n")
cat("   - Target FDR:", attr(selected, "fdr"), "\n")
cat("   - Number selected:", length(selected), "\n")

if (length(selected) > 0) {
  cat("   - Selected SNPs:", selected, "\n")
  # Calculate overlap with significant original SNPs for evaluation
  significant_original <- which(original_pvals < 0.05)
  cat("   - Overlap with significant originals:", sum(selected %in% significant_original), "\n")
}
cat("\n")

# 8. Visualization
cat("8. Creating visualization...\n")
if (requireNamespace("graphics", quietly = TRUE)) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists("plots")) dir.create("plots")
  
  png("plots/w_statistics_real_data.png", width = 800, height = 600)
  par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
  
  # Plot 1: W statistics
  significant_original <- which(original_pvals < 0.05)
  colors <- rep("black", length(w_stats))
  colors[selected] <- "red"
  colors[significant_original] <- "blue"
  colors[intersect(selected, significant_original)] <- "purple"
  
  plot(seq_along(w_stats), w_stats, 
       xlab = "SNP Index", ylab = "W Statistic",
       main = paste("W Statistics for Variable Selection (", p, "real SNPs)"),
       pch = 16, col = colors, cex = 1.2)
  
  if (is.finite(attr(selected, "threshold"))) {
    abline(h = attr(selected, "threshold"), col = "green", lty = 2, lwd = 2)
  }
  abline(h = 0, col = "gray", lty = 1)
  
  legend("topright", 
         c("True Positive", "False Positive", "True Negative", "False Negative", "Threshold"), 
         col = c("purple", "red", "black", "blue", "green"), 
         pch = c(16, 16, 16, 16, NA), 
         lty = c(NA, NA, NA, NA, 2),
         cex = 0.8)
  
  # Plot 2: P-value comparison  
  first_knockoff_pvals <- knockoff_results[[1]]$p_values
  plot(-log10(original_pvals), -log10(first_knockoff_pvals),
       xlab = "-log10(Original P-values)", ylab = "-log10(Knockoff P-values)",
       main = "Original vs Knockoff P-values (First Knockoff)",
       pch = 16, col = colors, cex = 1.2)
  abline(0, 1, col = "gray", lty = 2)
  
  dev.off()
  cat("   - Plot saved to: plots/w_statistics_real_data.png\n")
}

cat("\n=== Example completed successfully! ===\n")
cat("Results summary:\n")
cat("- Used real PLINK data from extdata\n")
cat("- Created", length(knockoffs$knockoffs), "knockoff matrices\n")
cat("- Analyzed", p, "genetic variants\n") 
cat("- Selected", length(selected), "significant variants\n")
cat("- Target FDR control at level", attr(selected, "fdr"), "\n")

# Demonstrate GDS file workflow
cat("\n=== GDS File Workflow Demonstration ===\n")
if (!is.null(knockoffs$gds_file) && file.exists(knockoffs$gds_file)) {
  cat("Testing GDS file loading...\n")
  
  # Load knockoffs from the saved GDS file
  loaded_knockoffs <- load_knockoff_gds(knockoffs$gds_file)
  cat("- Successfully loaded", length(loaded_knockoffs), "knockoff matrices from GDS\n")
  cat("- GDS file path:", knockoffs$gds_file, "\n")
  
  # Verify data integrity
  is_equal <- isTRUE(all.equal(knockoffs$knockoffs[[1]], loaded_knockoffs$knockoffs[[1]]))
  if (is_equal) {
    cat("- Data integrity verified: loaded knockoffs match original\n")
  } else {
    cat("- Warning: loaded knockoffs differ from original\n")
  }
} else {
  cat("No GDS file was created in this run.\n")
}

cat("\nFor more details, see the package documentation:\n")
cat("help(package = 'CoxMK')\n")
