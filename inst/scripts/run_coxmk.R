#!/usr/bin/env Rscript

# CoxMK: Cox Model-X Knockoffs for Survival Analysis
# Command-line interface script
# 
# This script provides a command-line interface for the CoxMK package,
# allowing users to run Cox Model-X knockoff analysis without writing R code.

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
})

# Define command line options
option_list <- list(
  make_option(c("--plink_prefix"), type="character", default=NULL,
              help="Path prefix to PLINK files (.bed, .bim, .fam)", metavar="PATH"),
  
  make_option(c("--phenotype"), type="character", default=NULL,
              help="Path to phenotype file (IID, time, status columns)", metavar="FILE"),
  
  make_option(c("--covariates"), type="character", default=NULL,
              help="Path to covariates file (IID + covariate columns)", metavar="FILE"),
  
  make_option(c("-M", "--M"), type="integer", default=5,
              help="Number of knockoffs to generate [default: %default]", metavar="INTEGER"),
  
  make_option(c("--fdr"), type="double", default=0.05,
              help="False discovery rate [default: %default]", metavar="DOUBLE"),
  
  make_option(c("--method"), type="character", default="median",
              help="W statistic method: 'median' or 'difference' [default: %default]", metavar="STRING"),
  
  make_option(c("--output_dir"), type="character", default="coxmk_results",
              help="Output directory [default: %default]", metavar="DIR"),
  
  make_option(c("--output_prefix"), type="character", default="analysis",
              help="Output file prefix [default: %default]", metavar="STRING"),
  
  make_option(c("--use_spacox"), action="store_true", default=FALSE,
              help="Use SPACox for null model fitting (if available)"),
  
  make_option(c("--save_knockoffs"), action="store_true", default=FALSE,
              help="Save knockoff matrices to file"),
  
  make_option(c("--save_gds"), action="store_true", default=FALSE,
              help="Save results in GDS format"),
  
  make_option(c("--verbose"), action="store_true", default=FALSE,
              help="Enable verbose output"),
  
  make_option(c("--test"), action="store_true", default=FALSE,
              help="Run with example data for testing"),
  
  make_option(c("--help_extended"), action="store_true", default=FALSE,
              help="Show extended help with examples")
)

# Parse command line arguments
opt_parser <- OptionParser(
  option_list=option_list,
  description="CoxMK: Cox Model-X Knockoffs for Survival Analysis",
  epilogue="For more information, visit: https://github.com/xiaoxiandadada/Cox-MK"
)

opt <- parse_args(opt_parser)

# Extended help function
show_extended_help <- function() {
  cat("\n=== CoxMK Extended Help ===\n\n")
  
  cat("DESCRIPTION:\n")
  cat("  CoxMK implements Cox regression with Model-X knockoffs for survival analysis,\n")
  cat("  providing finite-sample FDR control in high-dimensional genetic studies.\n\n")
  
  cat("BASIC USAGE:\n")
  cat("  Rscript run_coxmk.R --plink_prefix data/sample --phenotype pheno.txt --covariates cov.txt\n\n")
  
  cat("EXAMPLES:\n")
  cat("  # Basic analysis with default parameters\n")
  cat("  Rscript run_coxmk.R --plink_prefix data/study --phenotype survival.txt --covariates baseline.txt\n\n")
  
  cat("  # Custom parameters\n")
  cat("  Rscript run_coxmk.R --plink_prefix data/study --phenotype survival.txt \\\n")
  cat("    --covariates baseline.txt --M 10 --fdr 0.01 --method sdp --output_dir results/\n\n")
  
  cat("  # Use SPACox for large datasets\n")
  cat("  Rscript run_coxmk.R --plink_prefix data/large_study --phenotype survival.txt \\\n")
  cat("    --covariates baseline.txt --use_spacox --save_gds\n\n")
  
  cat("  # Test with example data\n")
  cat("  Rscript run_coxmk.R --test\n\n")
  
  cat("INPUT FILE FORMATS:\n")
  cat("  Phenotype file (tab-separated):\n")
  cat("    IID      time    status\n")
  cat("    sample1  10.5    1\n")
  cat("    sample2  8.2     0\n\n")
  
  cat("  Covariates file (tab-separated):\n")
  cat("    IID      age  sex  pc1    pc2\n")
  cat("    sample1  45   1    0.12  -0.05\n")
  cat("    sample2  52   0    0.08   0.03\n\n")
  
  cat("OUTPUT FILES:\n")
  cat("  - [prefix]_selected_snps.txt: Selected SNP information\n")
  cat("  - [prefix]_w_statistics.txt: W statistics for all SNPs\n")
  cat("  - [prefix]_summary.txt: Analysis summary\n")
  cat("  - [prefix]_knockoffs.rds: Knockoff matrices (if --save_knockoffs)\n\n")
  
  cat("For more details, see: https://github.com/xiaoxiandadada/Cox-MK\n\n")
}

# Show extended help if requested
if (opt$help_extended) {
  show_extended_help()
  quit(status=0)
}

# Test mode with example data
if (opt$test) {
  cat("=== Running CoxMK Test Mode ===\n")
  
  # Try to load CoxMK package
  if (!requireNamespace("CoxMK", quietly = TRUE)) {
    cat("Error: CoxMK package not found. Please install it first:\n")
    cat("  devtools::install_github('xiaoxiandadada/Cox-MK')\n")
    quit(status=1)
  }
  
  library(CoxMK)
  
  # Use example data
  extdata_path <- system.file("extdata", package = "CoxMK")
  
  if (!file.exists(file.path(extdata_path, "sample.bed"))) {
    cat("Error: Example data not found in package installation.\n")
    quit(status=1)
  }
  
  cat("Using example data from:", extdata_path, "\n")
  
  opt$plink_prefix <- file.path(extdata_path, "sample")
  opt$phenotype <- file.path(extdata_path, "tte_phenotype.txt")
  opt$covariates <- file.path(extdata_path, "covariates.txt")
  opt$output_dir <- "coxmk_test_results"
  opt$M <- 3  # Smaller for testing
  
  cat("Test parameters:\n")
  cat("  PLINK prefix:", opt$plink_prefix, "\n")
  cat("  Phenotype:", opt$phenotype, "\n")
  cat("  Covariates:", opt$covariates, "\n")
  cat("  M =", opt$M, ", FDR =", opt$fdr, "\n")
}

# Validate required arguments (unless in test mode)
if (!opt$test) {
  if (is.null(opt$plink_prefix) || is.null(opt$phenotype)) {
    cat("Error: --plink_prefix and --phenotype are required.\n")
    cat("Use --help for usage information or --test to run with example data.\n")
    quit(status=1)
  }
}

# Validate input files exist
validate_files <- function() {
  required_files <- c(
    paste0(opt$plink_prefix, ".bed"),
    paste0(opt$plink_prefix, ".bim"), 
    paste0(opt$plink_prefix, ".fam"),
    opt$phenotype
  )
  
  if (!is.null(opt$covariates)) {
    required_files <- c(required_files, opt$covariates)
  }
  
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0) {
    cat("Error: The following required files are missing:\n")
    for (f in missing_files) {
      cat("  -", f, "\n")
    }
    quit(status=1)
  }
}

validate_files()

# Load CoxMK package
if (!requireNamespace("CoxMK", quietly = TRUE)) {
  cat("Error: CoxMK package not found. Please install it first:\n")
  cat("  devtools::install_github('xiaoxiandadada/Cox-MK')\n")
  quit(status=1)
}

library(CoxMK)

# Create output directory
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
  if (opt$verbose) cat("Created output directory:", opt$output_dir, "\n")
}

# Print analysis parameters
cat("\n=== CoxMK Analysis Parameters ===\n")
cat("PLINK prefix:", opt$plink_prefix, "\n")
cat("Phenotype file:", opt$phenotype, "\n")
cat("Covariates file:", ifelse(is.null(opt$covariates), "None", opt$covariates), "\n")
cat("Number of knockoffs (M):", opt$M, "\n")
cat("FDR level:", opt$fdr, "\n")
cat("W statistic method:", opt$method, "\n")
cat("Output directory:", opt$output_dir, "\n")
cat("Use SPACox:", opt$use_spacox, "\n")
cat("Save knockoffs:", opt$save_knockoffs, "\n")
cat("Save GDS:", opt$save_gds, "\n")
cat("================================\n\n")

# Start timing
start_time <- Sys.time()

# Run the analysis
tryCatch({
  cat("Starting Cox Model-X knockoff analysis...\n")
  
  # Load phenotype data
  if (opt$verbose) cat("Loading phenotype data...\n")
  pheno_data <- read.table(opt$phenotype, header = TRUE, stringsAsFactors = FALSE)
  
  # Check required columns
  required_cols <- c("IID", "time", "status")
  missing_cols <- required_cols[!required_cols %in% colnames(pheno_data)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns in phenotype file: ", paste(missing_cols, collapse = ", "))
  }
  
  # Load covariates if provided
  covariates <- NULL
  if (!is.null(opt$covariates)) {
    if (opt$verbose) cat("Loading covariates...\n")
    cov_data <- read.table(opt$covariates, header = TRUE, stringsAsFactors = FALSE)
    # Remove IID column and use as matrix
    covariates <- as.matrix(cov_data[, !colnames(cov_data) %in% "IID", drop = FALSE])
  }
  
  # Run main analysis
  result <- cox_knockoff_analysis(
    plink_prefix = opt$plink_prefix,
    time = pheno_data$time,
    status = pheno_data$status,
    covariates = covariates,
    sample_ids = pheno_data$IID,
    M = opt$M,
    fdr = opt$fdr,
    method = opt$method,
    output_dir = opt$output_dir
  )
  
  # Save results
  output_prefix <- file.path(opt$output_dir, opt$output_prefix)
  
  # Save selected SNPs
  if (length(result$selected_vars) > 0) {
    selected_file <- paste0(output_prefix, "_selected_snps.txt")
    write.table(data.frame(
      SNP_Index = result$selected_vars,
      stringsAsFactors = FALSE
    ), file = selected_file, 
    row.names = FALSE, quote = FALSE, sep = "\t")
    cat("Selected SNPs saved to:", selected_file, "\n")
  } else {
    cat("No SNPs selected at FDR =", opt$fdr, "\n")
  }
  
  # Save W statistics
  w_stats_file <- paste0(output_prefix, "_w_statistics.txt")
  write.table(data.frame(
    SNP_Index = seq_along(result$W_stats),
    W_statistic = result$W_stats,
    stringsAsFactors = FALSE
  ), file = w_stats_file,
  row.names = FALSE, quote = FALSE, sep = "\t")
  cat("W statistics saved to:", w_stats_file, "\n")
  
  # Save summary
  summary_file <- paste0(output_prefix, "_summary.txt")
  summary_text <- c(
    "=== CoxMK Analysis Summary ===",
    paste("Analysis completed:", Sys.time()),
    paste("Input PLINK prefix:", opt$plink_prefix),
    paste("Phenotype file:", opt$phenotype),
    paste("Covariates file:", ifelse(is.null(opt$covariates), "None", opt$covariates)),
    paste("Number of knockoffs (M):", opt$M),
    paste("FDR level:", opt$fdr),
    paste("Knockoff method:", opt$method),
    "",
    "=== Results ===",
    paste("Number of SNPs selected:", length(result$selected_vars)),
    paste("Total runtime:", round(difftime(Sys.time(), start_time, units = "mins"), 2), "minutes"),
    ""
  )
  
  if (length(result$selected_vars) > 0) {
    summary_text <- c(summary_text,
      "Selected SNP indices:",
      paste(result$selected_vars, collapse = ", ")
    )
  }
  
  writeLines(summary_text, summary_file)
  cat("Analysis summary saved to:", summary_file, "\n")
  
  # Save knockoffs if requested
  if (opt$save_knockoffs && !is.null(result$knockoffs)) {
    knockoffs_file <- paste0(output_prefix, "_knockoffs.rds")
    saveRDS(result$knockoffs, knockoffs_file)
    cat("Knockoff matrices saved to:", knockoffs_file, "\n")
  }
  
  # Print final summary
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "mins")
  
  cat("\n=== Analysis Complete ===\n")
  cat("Selected SNPs:", length(result$selected_vars), "\n")
  cat("Total runtime:", round(runtime, 2), "minutes\n")
  cat("Results saved in:", opt$output_dir, "\n")
  
  if (opt$test) {
    cat("\n=== Test Mode Complete ===\n")
    cat("CoxMK package is working correctly!\n")
    cat("You can now use it with your own data.\n")
  }
  
}, error = function(e) {
  cat("Error during analysis:\n")
  cat(conditionMessage(e), "\n")
  quit(status=1)
})

cat("\nCoxMK analysis finished successfully!\n")
