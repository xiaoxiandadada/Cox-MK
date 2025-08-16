#' @title CoxMK: Cox Regression with Multiple Knockoffs
#' @name CoxMK
#' @description
#' Main interface functions for Cox regression analysis with Multiple knockoffs.
#' This package provides a complete workflow for survival analysis with 
#' variable selection using the multiple knockoffs methodology.
#' 
#' The workflow follows four main steps:
#' 1. **Generate Knockoffs**: Create knockoff variables using \code{\link{create_knockoffs}}
#' 2. **Fit Null Model**: Fit null Cox model using \code{\link{fit_null_cox_model}}
#' 3. **Perform Testing**: Conduct association testing using \code{\link{perform_association_testing}}
#' 4. **Apply Filter**: Select variables using \code{\link{knockoff_filter}}
#' 
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{cox_knockoff_analysis}} - Complete knockoff analysis workflow
#'   \item \code{\link{create_knockoffs}} - Step 1: Generate knockoff variables
#'   \item \code{\link{fit_null_cox_model}} - Step 2: Fit null Cox model for testing
#'   \item \code{\link{perform_association_testing}} - Step 3: Perform association testing
#'   \item \code{\link{knockoff_filter}} - Step 4: Apply knockoff filter for variable selection
#' }
#' 
#' @importFrom stats as.formula complete.cases median na.omit pchisq var
#' @importFrom Matrix Matrix
#' @importFrom survival coxph Surv
#' @importFrom irlba irlba
NULL

#' Complete Cox Knockoff Analysis Workflow
#'
#' Performs a complete Multiple knockoff analysis following the four-step workflow:
#' 1. Generate knockoff variables from PLINK data and save to GDS format
#' 2. Fit null Cox model using SPACox for efficient large-scale analysis
#' 3. Perform SPA testing using original and knockoff variables
#' 4. Apply knockoff filter for variable selection with FDR control
#'
#' @param plink_prefix Path prefix for PLINK files (without extension)
#' @param time Survival times
#' @param status Event indicators (1=event, 0=censored)
#' @param covariates Optional covariate matrix/data.frame
#' @param sample_ids Sample IDs (optional, will be generated from .fam file)
#' @param null_model Pre-fitted null Cox model or path to RDS file with fitted model (optional)
#' @param gds_file Path to pre-generated GDS file with knockoff data (optional, if provided, knockoffs will be loaded instead of generated)
#' @param M Number of knockoff copies to generate (default: 5)
#' @param fdr Target false discovery rate (default: 0.05)
#' @param method Statistical method for W statistics ("median", "difference")
#' @param output_dir Directory to save GDS files (default: extdata folder)
#' @return List containing:
#'   \item{selected_vars}{Indices of selected variables}
#'   \item{W_stats}{W statistics for all variables}
#'   \item{threshold}{Knockoff threshold used}
#'   \item{gds_file}{Path to GDS file used}
#'   \item{null_model}{Fitted null Cox model}
#'   \item{test_results}{SPA test results}
#' @export
#' @examples
#' \dontrun{
#' # Standard workflow with PLINK data
#' extdata_path <- system.file('extdata', package = 'CoxMK')
#' plink_prefix <- file.path(extdata_path, 'sample')
#' pheno_data <- prepare_phenotype(file.path(extdata_path, 'tte_phenotype.txt'))
#' covar_data <- load_covariates(file.path(extdata_path, 'covariates.txt'))
#' 
#' # Option 1: Complete analysis in one step
#' result <- cox_knockoff_analysis(
#'   plink_prefix = plink_prefix,
#'   time = pheno_data$time,
#'   status = pheno_data$status,
#'   covariates = covar_data,
#'   M = 3,
#'   fdr = 0.1
#' )
#' 
#' # Option 2: Step-by-step workflow with GDS file reuse
#' # Step 2a: Generate knockoffs and save to GDS file
#' knockoffs <- create_knockoffs(
#'   X = load_plink_data(plink_prefix)$genotypes,
#'   pos = load_plink_data(plink_prefix)$positions,
#'   M = 3
#' )
#' gds_file <- knockoffs$gds_file  # GDS file path
#' 
#' # Step 2b: Fit null model separately  
#' null_model <- fit_null_cox_model(
#'   time = pheno_data$time,
#'   status = pheno_data$status, 
#'   covariates = covar_data
#' )
#' 
#' # Step 2c: Load knockoffs from GDS file and run analysis
#' result <- cox_knockoff_analysis(
#'   plink_prefix = plink_prefix,
#'   time = pheno_data$time,
#'   status = pheno_data$status,
#'   covariates = covar_data,
#'   null_model = null_model,
#'   gds_file = gds_file,  # Use pre-generated GDS file
#'   M = 3,
#'   fdr = 0.05
#' )
#' # View selected variables
#' print(result$selected_vars)
#' print(result$summary)
#' }
cox_knockoff_analysis <- function(plink_prefix, time, status,
                                 covariates = NULL, sample_ids = NULL,
                                 null_model = NULL, gds_file = NULL, M = 5, fdr = 0.05, method = "median", 
                                 output_dir = NULL) {

  # Validate input parameters
  if (!file.exists(paste0(plink_prefix, ".bed")) || 
      !file.exists(paste0(plink_prefix, ".bim")) ||
      !file.exists(paste0(plink_prefix, ".fam"))) {
    stop("PLINK files (.bed, .bim, .fam) not found with prefix: ", plink_prefix)
  }
  
  if (length(time) != length(status)) {
    stop("time and status must have the same length")
  }
  
  # Setup output directory for GDS files
  if (is.null(output_dir)) {
    output_dir <- system.file("extdata", package = "CoxMK")
    if (!dir.exists(output_dir)) {
      output_dir <- file.path(getwd(), "extdata")
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  cat("=== Cox Multiple Knockoff Analysis Workflow ===\n")
  
  # Step 1: Generate knockoff variables from PLINK data
  cat("\n1. GENERATING KNOCKOFF VARIABLES FROM PLINK DATA\n")
  cat("   Loading PLINK data from:", plink_prefix, "\n")
  
  # Load PLINK data
  plink_data <- load_plink_data(plink_prefix)
  X_original <- plink_data$genotypes
  pos <- plink_data$positions
  if (is.null(sample_ids)) {
    sample_ids <- plink_data$sample_ids
  }
  
  # Extract chromosome information
  bim_file <- paste0(plink_prefix, ".bim")
  bim_data <- read.table(bim_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(bim_data) <- c("chr", "snp", "cM", "pos", "A1", "A2")
  chromosomes <- unique(bim_data$chr)
  
  # Step 1a: Check if GDS file is provided, if so load knockoffs from it
  if (!is.null(gds_file) && file.exists(gds_file)) {
    cat("   Loading pre-generated knockoffs from GDS file:", gds_file, "\n")
    knockoff_data <- load_knockoff_gds(gds_file)
    X_knockoffs <- knockoff_data$knockoffs
    gds_output <- gds_file
    
    # Validate that loaded data matches current PLINK data
    if (ncol(X_original) != ncol(knockoff_data$original)) {
      stop("Number of SNPs in GDS file (", ncol(knockoff_data$original), 
           ") does not match PLINK data (", ncol(X_original), ")")
    }
    if (nrow(X_original) != nrow(knockoff_data$original)) {
      stop("Number of samples in GDS file (", nrow(knockoff_data$original), 
           ") does not match PLINK data (", nrow(X_original), ")")
    }
    
    M <- length(X_knockoffs)
    cat("   - Loaded", M, "knockoff matrices from GDS file\n")
    cat("   - Data dimensions:", nrow(X_original), "samples x", ncol(X_original), "SNPs\n")
    
  } else {
    # Step 1b: Generate new knockoffs
    cat("   Creating", M, "knockoff copies...\n")
    cat("   Data dimensions:", nrow(X_original), "samples x", ncol(X_original), "SNPs\n")
    cat("   Chromosomes found:", paste(chromosomes, collapse = ", "), "\n")
    
    # Generate knockoffs
    knockoff_result <- create_knockoffs(X = X_original, pos = pos, M = M)
    X_knockoffs <- knockoff_result$knockoffs
    gds_output <- knockoff_result$gds_file
  }
  
  cat("   - Knockoff generation/loading complete!\n")

  # Validate sample sizes
  n_samples <- length(time)
  if (nrow(X_original) != n_samples) {
    stop("Number of genotype samples (", nrow(X_original), 
         ") does not match phenotype samples (", n_samples, ")")
  }

  # Step 2: Fit Null Model
  cat("\n2. Fit Null Model\n")

  # Handle null_model parameter (can be object, file path, or NULL)
  if (is.null(null_model)) {
    # Fit new model
    if (requireNamespace("SPACox", quietly = TRUE)) {
      cat("Fitting null Cox model using SPACox with", ifelse(is.null(covariates), "no", ncol(covariates)), "covariates...\n")
      null_model <- fit_null_cox_model(time = time, status = status, covariates = covariates)
      cat("- SPACox null model fitted successfully!\n")
    } else {
      cat("SPACox not available, will use traditional Cox regression...\n")
      cat("- Traditional Cox regression will be used for testing\n")
    }
  } else if (is.character(null_model) && length(null_model) == 1) {
    # Load model from file
    cat("   Loading pre-fitted null model from:", null_model, "\n")
    if (!file.exists(null_model)) {
      stop("Model file not found: ", null_model)
    }

    model_info <- readRDS(null_model)
    if (is.list(model_info) && "model" %in% names(model_info)) {
      # Model saved with metadata (from fit_model.R script)
      null_model <- model_info$model
      cat("- Model loaded successfully (Type:", model_info$model_type, ")\n")
      cat("- Original samples:", model_info$n_samples, ", Events:", model_info$n_events, "\n")
      if (length(model_info$covariate_names) > 0) {
        cat("- Covariates:", paste(model_info$covariate_names, collapse = ", "), "\n")
      }
    } else {
      # Model saved directly
      null_model <- model_info
      cat("- Model loaded successfully\n")
    }
  } else {
    # Use provided model object
    cat("Using provided null model object\n")
  }

  # Step 3: SPA testing and association analysis
  cat("\n3. SPA TESTING AND ASSOCIATION ANALYSIS\n")

  cat("   Testing original variables...\n")
  orig_results <- perform_association_testing(X_original, null_model)

  cat("   Testing knockoff variables...\n")
  M <- length(X_knockoffs)
  knockoff_results <- vector("list", M)
  
  for (k in seq_len(M)) {
    cat("     Knockoff copy", k, "/", M, "\n")
    knockoff_results[[k]] <- perform_association_testing(X_knockoffs[[k]], null_model)
  }
  
  test_results <- list(
    original = orig_results,
    knockoffs = knockoff_results
  )
  
  cat("   - Association testing complete!\n")
  
  # Step 4: Apply knockoff filter
  cat("\n4. APPLYING KNOCKOFF FILTER\n")
  cat("   Computing W statistics using method:", method, "\n")
  
  # Extract test statistics
  t_orig <- orig_results$test_stats
  
  # For multiple knockoffs, combine test statistics (use median as example)
  t_knock_matrix <- do.call(cbind, lapply(knockoff_results, function(x) x$test_stats))
  t_knock_combined <- apply(t_knock_matrix, 1, median)
  
  # Calculate W statistics using the same method as calculate_w_statistics
  W_stats <- calculate_w_statistics(t_orig, t_knock_combined, method = method)
  
  cat("   Applying knockoff filter with FDR =", fdr, "\n")
  selected_vars <- knockoff_filter(W_stats, fdr = fdr)
  threshold <- attr(selected_vars, "threshold")

  filter_results <- list(
    selected_vars = selected_vars,
    W_stats = W_stats,
    threshold = threshold
  )
  
  cat("   - Variable selection complete!\n")
  
  # Print analysis summary
  cat("\n=== ANALYSIS SUMMARY ===\n")
  cat("   Total variables tested:", length(filter_results$W_stats), "\n")
  cat("   Variables selected:", length(filter_results$selected_vars), "\n")
  cat("   Selection proportion:", round(length(filter_results$selected_vars) / 
                                        length(filter_results$W_stats) * 100, 2), "%\n")
  cat("   Threshold used:", round(filter_results$threshold, 4), "\n")
  if (!is.null(gds_output)) {
    cat("   Knockoff data saved to:", gds_output, "\n")
  }
  
  # Return results
  results <- list(
    selected_vars = filter_results$selected_vars,
    W_stats = filter_results$W_stats,
    threshold = filter_results$threshold,
    gds_file = gds_output,
    null_model = null_model,
    test_results = test_results,
    method = method,
    fdr = fdr,
    summary = list(
      n_variables = length(filter_results$W_stats),
      n_selected = length(filter_results$selected_vars),
      selection_rate = length(filter_results$selected_vars) / length(filter_results$W_stats)
    )
  )
  
  return(results)
}
