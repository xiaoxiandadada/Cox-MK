# CoxMK Package Complete Workflow Example
library(CoxMK)
library(survival)
library(Matrix)

cat("=== CoxMK Package Complete Workflow Example ===\n\n")

# 1. Load PLINK data
cat("1. Loading PLINK data...\n")
extdata_path <- system.file("extdata", package = "CoxMK")

plink_data <- load_plink_data(file.path(extdata_path, "sample"))
phenotype_data <- read.table(file.path(extdata_path, "tte_phenotype.txt"), 
                            header = TRUE, stringsAsFactors = FALSE)
covariate_data <- read.table(file.path(extdata_path, "covariates.txt"), 
                            header = TRUE, stringsAsFactors = FALSE)

cat("   - PLINK data loaded successfully\n")
cat("   - Genotypes:", dim(plink_data$genotypes), "(samples x SNPs)\n")
cat("   - Phenotypes:", nrow(phenotype_data), "samples\n")
cat("   - Covariates:", nrow(covariate_data), "samples\n\n")

# 2. Create knockoffs
cat("2. Creating knockoff variables...\n")
knockoffs <- create_knockoffs(
  X = plink_data$genotypes,
  pos = plink_data$positions,
  chr_info = plink_data$chr_info,
  sample_ids = plink_data$sample_ids,
  M = 5,
  corr_max = 0.75,
  maxBP.neighbor = 1e5
)

cat("   - Created", length(knockoffs$knockoffs), "knockoff matrices\n")
if (!is.null(knockoffs$gds_file)) {
  cat("   - Knockoffs saved to GDS file:", basename(knockoffs$gds_file), "\n")
}
cat("\n")

# 3. Prepare phenotype data
cat("3. Preparing phenotype data for Cox regression...\n")

pheno_data <- prepare_phenotype(
  file.path(extdata_path, "tte_phenotype.txt"),
  time_col = "time",
  status_col = "status"
)

covariates <- load_covariates(file.path(extdata_path, "covariates.txt"))

cat("   - Time to event range:", range(pheno_data$time), "days\n")
cat("   - Events observed:", sum(pheno_data$status), "out of", nrow(pheno_data), "samples\n\n")

# 4. Fit null Cox model
cat("4. Fitting null Cox model...\n")

null_model <- fit_null_cox_model(
  time = pheno_data$time,
  status = pheno_data$status,
  covariates = covariates[, setdiff(names(covariates), "IID"), drop = FALSE]
)

cat("   - Null Cox model fitted successfully\n")
coef_summary <- coef(null_model)
for (i in seq_along(coef_summary)) {
  cat("     -", names(coef_summary)[i], ":", round(coef_summary[i], 4), "\n")
}

null_model_file <- file.path(getwd(), "inst", "extdata", "null_model.rds")
if (!dir.exists(dirname(null_model_file))) {
  dir.create(dirname(null_model_file), recursive = TRUE)
}
saveRDS(null_model, null_model_file)
cat("   - Null model saved to:", basename(null_model_file), "\n\n")

# 5. Perform association testing
cat("5. Performing association testing...\n")

cat("   - Testing original data...\n")
original_results <- perform_association_testing(plink_data$genotypes, null_model)

cat("   - Testing", length(knockoffs$knockoffs), "knockoff matrices...\n")
knockoff_results <- vector("list", length(knockoffs$knockoffs))

for (k in seq_along(knockoffs$knockoffs)) {
  cat("     - Testing knockoff matrix", k, "of", length(knockoffs$knockoffs), "...\n")
  knockoff_results[[k]] <- perform_association_testing(knockoffs$knockoffs[[k]], null_model)
}

original_stats <- original_results$test_stats
knockoff_stats <- lapply(knockoff_results, function(x) x$test_stats)
cat("   - Association testing completed\n\n")

# 6. Apply knockoff filter
cat("6. Applying knockoff filter for variable selection...\n")

if (is.list(knockoff_stats)) {
  knockoff_matrix <- do.call(cbind, knockoff_stats)
  knockoff_stats_final <- rowMeans(knockoff_matrix)
} else {
  knockoff_stats_final <- knockoff_stats
}

# Use calculate_w_statistics instead of mk_statistic directly
W_stats <- calculate_w_statistics(
  original_stats,
  knockoff_stats_final,
  method = "median"
)

selected_results <- knockoff_filter(W_stats, fdr = 0.05)

if (is.logical(selected_results)) {
  selected_indices <- which(selected_results)
  n_selected <- sum(selected_results)
} else if (is.numeric(selected_results)) {
  selected_indices <- selected_results
  n_selected <- length(selected_results)
} else {
  selected_indices <- integer(0)
  n_selected <- 0
}

cat("   - Variables selected:", n_selected, "out of", ncol(plink_data$genotypes), "\n")
cat("   - Selection rate:", round(100 * n_selected / ncol(plink_data$genotypes), 2), "%\n\n")

# 7. Save results
cat("7. Saving analysis results...\n")

# Create summary following main.R format
summary_stats <- list(
  n_variables = length(W_stats),
  n_selected = n_selected,
  selection_rate = n_selected / length(W_stats)
)

# Prepare complete results
analysis_results <- list(
  selected_vars = selected_indices,
  W_stats = W_stats,
  original_results = original_results,
  knockoff_results = knockoff_results,
  null_model = null_model,
  summary = summary_stats,
  data_info = list(
    n_samples = nrow(plink_data$genotypes),
    n_snps = ncol(plink_data$genotypes),
    n_events = sum(pheno_data$status),
    event_rate = sum(pheno_data$status) / nrow(pheno_data),
    fdr_level = 0.05
  )
)

# Save results to RDS file
results_file <- file.path(getwd(), "inst", "extdata", "analysis_results.rds")
saveRDS(analysis_results, results_file)
cat("   - Analysis results saved to:", basename(results_file), "\n")

# Save selected SNP information to text file
if (n_selected > 0) {
  selected_snps <- data.frame(
    index = selected_indices,
    position = plink_data$positions[selected_indices],
    W_statistic = W_stats[selected_indices],
    p_value = original_results$p_values[selected_indices]
  )
  
  selected_file <- file.path(getwd(), "inst", "extdata", "selected_snps.txt")
  write.table(selected_snps, selected_file, row.names = FALSE, quote = FALSE, sep = "\t")
  cat("   - Selected SNPs saved to:", basename(selected_file), "\n")
}
cat("\n")

# 8. Analysis Summary
cat("8. Analysis Summary\n")
cat("================================================================================\n")
cat("Data processed:\n")
cat("- Samples:", analysis_results$data_info$n_samples, "\n")
cat("- SNPs:", analysis_results$data_info$n_snps, "\n")
cat("- Events:", analysis_results$data_info$n_events, 
    "(", round(100 * analysis_results$data_info$event_rate, 1), "%)\n")
cat("- FDR level:", analysis_results$data_info$fdr_level, "\n\n")

cat("Results:\n")
cat("- Variables selected:", analysis_results$summary$n_selected, "\n")
cat("- Selection rate:", round(100 * analysis_results$summary$selection_rate, 2), "%\n")

if (!is.null(knockoffs$gds_file)) {
  cat("- Knockoff data saved to:", basename(knockoffs$gds_file), "\n")
}

cat("\n=== CoxMK Workflow Completed Successfully ===\n")
cat("Files saved to inst/extdata/:\n")
cat("✓ null_model.rds - Fitted Cox model\n")
cat("✓ analysis_results.rds - Complete analysis results\n")
if (n_selected > 0) {
  cat("✓ selected_snps.txt - Selected SNP information\n")
}
if (!is.null(knockoffs$gds_file)) {
  cat("✓", basename(knockoffs$gds_file), "- Knockoff variables\n")
}
cat("\n")
