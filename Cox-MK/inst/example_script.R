# CoxKnockoff Package Example
# This script demonstrates the complete workflow

# Load required packages
library(CoxMK)
library(survival)
library(Matrix)

cat("=== CoxMK Package Example ===\n\n")

# 1. Load example data
cat("1. Loading example data...\n")
data(example_genotypes)
data(example_positions)
data(example_phenotype)
data(example_covariates)

cat("   - Genotypes:", dim(example_genotypes), "(samples x SNPs)\n")
cat("   - Positions:", length(example_positions), "SNPs\n")
cat("   - Phenotypes:", nrow(example_phenotype), "samples\n")
cat("   - Covariates:", nrow(example_covariates), "samples\n\n")

# 2. Create knockoff variables (default M=5, user can specify different values)
cat("2. Creating knockoff variables...\n")
knockoffs <- create_knockoffs(
  X = example_genotypes,
  pos = example_positions,
  M = 5,  # Default 5 knockoffs, users can change this (e.g., M = 3, 10, etc.)
  corr_max = 0.75,
  maxBP.neighbor = 1e5
)

cat("   - Created", length(knockoffs), "knockoff matrices (M=5 by default)\n")
cat("   - Users can specify different M values (e.g., M=3, M=10)\n")
cat("   - Each matrix dimensions:", dim(knockoffs[[1]]), "\n\n")

# 3. Prepare phenotype data
cat("3. Preparing phenotype data...\n")
pheno_data <- merge(example_phenotype, example_covariates, by = c("FID", "IID"))
covariates <- pheno_data[, c("age", "sex", "bmi")]

cat("   - Merged data dimensions:", dim(pheno_data), "\n")
cat("   - Covariates used:", names(covariates), "\n\n")

# 4. Fit null model (with SPACox option)
cat("4. Fitting null Cox model...\n")
cat("   - Standard coxph method (default)\n")
cat("   - For large-scale studies, consider SPACox: https://github.com/WenjianBI/SPACox\n")

null_model <- fit_null_model(
  time = pheno_data$time,
  status = pheno_data$status,
  covariates = covariates,
  use_spacox = FALSE  # Set to TRUE if SPACox is installed
)

cat("   - Null model fitted successfully using", attr(null_model, "method") %||% "coxph", "\n")
cat("   - Number of events:", sum(pheno_data$status), "\n")
cat("   - Median survival time:", median(pheno_data$time), "days\n\n")

# 5. Simulate knockoff screening results (for demonstration)
cat("5. Performing knockoff screening (simulated for demo)...\n")
set.seed(123)
p <- ncol(example_genotypes)

# Simulate some significant associations
significant_snps <- sample(1:p, size = 5)
original_pvals <- runif(p, 0.1, 1)
original_pvals[significant_snps] <- runif(5, 0.001, 0.05)  # Make some significant
original_coefs <- rnorm(p, 0, 0.3)
original_coefs[significant_snps] <- rnorm(5, 0, 0.8)  # Larger effects

# Knockoff results (less significant on average)
knockoff_pvals <- lapply(1:5, function(k) {  # Changed to 5 knockoffs
  kp <- runif(p, 0.1, 1)
  kp[significant_snps] <- runif(5, 0.05, 0.3)  # Less significant
  kp
})

knockoff_coefs <- lapply(1:5, function(k) {  # Changed to 5 knockoffs
  kc <- rnorm(p, 0, 0.3)
  kc[significant_snps] <- rnorm(5, 0, 0.4)  # Smaller effects
  kc
})

cat("   - Simulated p-values for", p, "SNPs\n")
cat("   - True significant SNPs:", significant_snps, "\n\n")

# 6. Calculate W statistics
cat("6. Calculating W statistics...\n")
w_stats <- calculate_w_statistics(
  original_pvals = original_pvals,
  knockoff_pvals = knockoff_pvals,
  original_coefs = original_coefs,
  knockoff_coefs = knockoff_coefs
)

cat("   - W statistics range:", round(range(w_stats), 3), "\n")
cat("   - Positive W stats:", sum(w_stats > 0), "\n")
cat("   - Negative W stats:", sum(w_stats < 0), "\n\n")

# 7. Apply knockoff filter
cat("7. Applying knockoff filter...\n")
selected <- knockoff_filter(w_stats, fdr = 0.1)

cat("   - Selection threshold:", round(selected$threshold, 3), "\n")
cat("   - Estimated FDP:", round(selected$fdp, 3), "\n")
cat("   - Number selected:", length(selected$selected), "\n")

if (length(selected$selected) > 0) {
  cat("   - Selected SNPs:", selected$selected, "\n")
  cat("   - True positives:", sum(selected$selected %in% significant_snps), "\n")
  cat("   - False positives:", sum(!(selected$selected %in% significant_snps)), "\n")
}
cat("\n")

# 8. Visualization
cat("8. Creating visualization...\n")
if (requireNamespace("graphics", quietly = TRUE)) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists("plots")) dir.create("plots")
  
  png("plots/w_statistics_example.png", width = 800, height = 600)
  par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
  
  # Plot 1: W statistics
  colors <- rep("black", length(w_stats))
  colors[selected$selected] <- "red"
  colors[significant_snps] <- "blue"
  colors[intersect(selected$selected, significant_snps)] <- "purple"
  
  plot(seq_along(w_stats), w_stats, 
       xlab = "SNP Index", ylab = "W Statistic",
       main = "W Statistics for Variable Selection",
       pch = 16, col = colors, cex = 1.2)
  
  if (is.finite(selected$threshold)) {
    abline(h = selected$threshold, col = "green", lty = 2, lwd = 2)
  }
  abline(h = 0, col = "gray", lty = 1)
  
  legend("topright", 
         c("True Positive", "False Positive", "True Negative", "False Negative", "Threshold"), 
         col = c("purple", "red", "black", "blue", "green"), 
         pch = c(16, 16, 16, 16, NA), 
         lty = c(NA, NA, NA, NA, 2),
         cex = 0.8)
  
  # Plot 2: P-value comparison
  plot(-log10(original_pvals), -log10(sapply(knockoff_pvals, function(x) x[1])),
       xlab = "-log10(Original P-values)", ylab = "-log10(Knockoff P-values)",
       main = "Original vs Knockoff P-values (First Knockoff)",
       pch = 16, col = colors, cex = 1.2)
  abline(0, 1, col = "gray", lty = 2)
  
  dev.off()
  cat("   - Plot saved to: plots/w_statistics_example.png\n")
}

cat("\n=== Example completed successfully! ===\n")
cat("Results summary:\n")
cat("- Created", length(knockoffs), "knockoff matrices\n")
cat("- Analyzed", p, "genetic variants\n") 
cat("- Selected", length(selected$selected), "significant variants\n")
cat("- Estimated FDR control at level", selected$fdp, "\n")

# Optionally run real analysis with the provided sample data
cat("\n=== Optional: Analysis with Real Sample Data ===\n")
sample_data_path <- system.file("extdata", package = "CoxKnockoff")

if (file.exists(file.path(sample_data_path, "sample.bed"))) {
  cat("Sample PLINK data found! You can run a real analysis with:\n\n")
  cat("# Load real sample data\n")
  cat("plink_data <- load_plink_data(file.path(sample_data_path, 'sample'))\n")
  cat("pheno <- read.table(file.path(sample_data_path, 'tte_phenotype.txt'), header = TRUE)\n")
  cat("covar <- read.table(file.path(sample_data_path, 'covariates.txt'), header = TRUE)\n")
  cat("\n# Run full analysis\n")
  cat("real_knockoffs <- create_knockoffs(plink_data$genotypes, plink_data$positions, M = 5)\n")
  cat("# ... continue with analysis\n")
} else {
  cat("No sample PLINK data found in package.\n")
  cat("Copy your PLINK files to run a real analysis.\n")
}

cat("\nFor more details, see the package vignette:\n")
cat("vignette('introduction', package = 'CoxKnockoff')\n")
