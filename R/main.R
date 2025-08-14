#' @title CoxMK: Cox Regression with Model-X Knockoffs
#' @name CoxMK
#' @description
#' Main interface functions for Cox regression analysis with Model-X knockoffs.
#' This package provides a complete workflow for survival analysis with 
#' variable selection using the knockoff methodology.
#' 
#' The workflow follows three main steps:
#' 1. **Generate Knockoffs**: Create knockoff variables using \code{\link{create_knockoffs}}
#' 2. **SPA Testing**: Perform association testing using \code{\link{fit_cox_spa}}
#' 3. **Apply Filter**: Select variables using \code{\link{knockoff_filter}}
#' 
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{cox_knockoff_analysis}} - Complete knockoff analysis workflow
#'   \item \code{\link{create_knockoffs}} - Step 1: Generate knockoff variables
#'   \item \code{\link{fit_cox_spa}} - Step 2: SPA testing and association analysis
#'   \item \code{\link{knockoff_filter}} - Step 3: Apply knockoff filter for variable selection
#' }
#' 
#' @importFrom stats as.formula complete.cases median na.omit pchisq var
#' @importFrom Matrix Matrix
#' @importFrom survival coxph Surv
#' @importFrom irlba irlba
NULL

#' Complete Cox Knockoff Analysis Workflow
#'
#' Performs a complete Model-X knockoff analysis following the three-step workflow:
#' 1. Generate knockoff variables and save to GDS format
#' 2. Perform SPA testing using GDS data
#' 3. Apply knockoff filter for variable selection
#'
#' @param X Genotype matrix (samples x SNPs) - only needed if gds_file is NULL
#' @param pos SNP positions vector - only needed if gds_file is NULL
#' @param time Survival times
#' @param status Event indicators (1=event, 0=censored)
#' @param gds_file Path to existing GDS file with knockoffs (optional)
#' @param sample_ids Sample IDs (optional, will be generated if NULL)
#' @param covariates Optional covariate matrix/data.frame
#' @param M Number of knockoff copies to generate (default: 5)
#' @param fdr Target false discovery rate (default: 0.1)
#' @param method Statistical method for W statistics ("difference", "mk_median")
#' @param use_spa Whether to use SPA test when available (default: TRUE)
#' @param save_gds Whether to save knockoffs to GDS format (default: TRUE)
#' @param output_dir Directory to save GDS files (default: extdata folder)
#' @return List containing:
#'   \item{selected_vars}{Indices of selected variables}
#'   \item{W_stats}{W statistics for all variables}
#'   \item{threshold}{Knockoff threshold used}
#'   \item{gds_file}{Path to GDS file used}
#'   \item{test_results}{SPA test results}
#' @export
#' @examples
#' \dontrun{
#' # Method 1: Standard workflow with PLINK data
#' extdata_path <- system.file('extdata', package = 'CoxMK')
#' plink_data <- load_plink_data(file.path(extdata_path, 'sample'))
#' pheno_data <- prepare_phenotype(file.path(extdata_path, 'tte_phenotype.txt'))
#' covar_data <- load_covariates(file.path(extdata_path, 'covariates.txt'))
#' 
#' # Generate knockoffs and save to GDS format
#' result1 <- cox_knockoff_analysis(
#'   X = plink_data$genotypes,
#'   pos = plink_data$positions,
#'   time = pheno_data$time,
#'   status = pheno_data$status,
#'   covariates = covar_data,
#'   M = 3,
#'   fdr = 0.1,
#'   save_gds = TRUE
#' )
#' 
#' # Method 2: Reuse existing GDS file
#' gds_file <- result1$gds_file
#' result2 <- cox_knockoff_analysis(
#'   gds_file = gds_file,
#'   time = pheno_data$time,
#'   status = pheno_data$status,
#'   covariates = covar_data,
#'   fdr = 0.05  # Different FDR threshold
#' )
#' 
#' # Method 3: Quick test with example data (no GDS)
#' data(example_genotypes)
#' data(example_positions)
#' data(example_phenotype)
#' 
#' result3 <- cox_knockoff_analysis(
#'   X = example_genotypes,
#'   pos = example_positions,
#'   time = example_phenotype$time,
#'   status = example_phenotype$status,
#'   save_gds = FALSE  # Skip GDS for quick testing
#' )
#' 
#' # View selected variables
#' print(result3$selected_vars)
#' }
cox_knockoff_analysis <- function(X = NULL, pos = NULL, time, status, 
                                 gds_file = NULL, sample_ids = NULL,
                                 covariates = NULL, M = 5, fdr = 0.1,
                                 method = "difference", use_spa = TRUE,
                                 save_gds = TRUE, output_dir = NULL) {
  
  # Validate input parameters
  if (is.null(gds_file) && (is.null(X) || is.null(pos))) {
    stop("Either gds_file must be provided, or both X and pos must be provided")
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
  
  cat("=== Cox Model-X Knockoff Analysis Workflow ===\n")
  
  # Step 1: Generate or load knockoff variables
  if (!is.null(gds_file) && file.exists(gds_file)) {
    cat("\n1. LOADING KNOCKOFFS FROM GDS FILE\n")
    cat("   Loading knockoffs from:", gds_file, "\n")
    knockoff_data <- load_knockoff_gds(gds_file)
    X_original <- knockoff_data$original
    X_knockoffs <- knockoff_data$knockoffs
    pos <- knockoff_data$positions
    sample_ids <- knockoff_data$sample_ids
    gds_output <- gds_file
    cat("   - Knockoff data loaded successfully!\n")
  } else {
    cat("\n1. GENERATING NEW KNOCKOFF VARIABLES\n")
    cat("   Creating", M, "knockoff copies...\n")
    if (is.null(sample_ids)) {
      sample_ids <- paste0("sample_", seq_len(nrow(X)))
    }
    
    result <- create_knockoffs(X = X, pos = pos, M = M)
    # Save to GDS format if requested
    if (save_gds) {
      # Generate descriptive filename based on genomic region
      # Try to get chromosome info from existing data (could be enhanced)
      chr <- "chr1"  # Default to chr1, could be passed as parameter in future
      
      start_pos <- min(pos)
      end_pos <- max(pos)
      
      gds_filename <- paste0(chr, "_", start_pos, "_", end_pos, "_knockoff.gds")
      gds_output <- file.path(output_dir, gds_filename)
      
      cat("   Saving knockoffs to GDS format:", gds_filename, "\n")
      
      # Create GDS file with knockoff data
      if (requireNamespace("gdsfmt", quietly = TRUE)) {
        gdsfile <- gdsfmt::createfn.gds(gds_output)
        
        # Add sample information
        gdsfmt::add.gdsn(gdsfile, "sample.id", sample_ids)
        gdsfmt::add.gdsn(gdsfile, "positions", pos)
        
        # Convert and add original genotype data
        X_dense <- safe_as_matrix(X, sparse = FALSE)
        gdsfmt::add.gdsn(gdsfile, "original", X_dense)
        
        # Convert and add knockoff data  
        knockoffs_dense <- array(dim = c(nrow(X), ncol(X), M))
        for (k in seq_len(M)) {
          knockoffs_dense[,,k] <- safe_as_matrix(result$knockoffs[[k]], sparse = FALSE)
        }
        gdsfmt::add.gdsn(gdsfile, "knockoffs", knockoffs_dense)
        
        gdsfmt::closefn.gds(gdsfile)
        cat("   * GDS file written successfully\n")
      } else {
        warning("gdsfmt package not available, knockoffs not saved to GDS format")
        gds_output <- NULL
      }
    } else {
      gds_output <- NULL
    }
    
    X_original <- X
    X_knockoffs <- result$knockoffs
    cat("   - Knockoff generation complete!\n")
  }
  
  # Validate sample consistency
  n_samples <- length(time)
  if (nrow(X_original) != n_samples) {
    stop("Number of genotype samples (", nrow(X_original), 
         ") does not match phenotype samples (", n_samples, ")")
  }
  
  # Step 2: Perform SPA testing
  cat("\n2. SPA TESTING AND ASSOCIATION ANALYSIS\n")
  
  # Test original variables
  cat("   Testing original variables...\n")
  orig_results <- fit_cox_spa(X_original, time, status, covariates, use_spa = use_spa)
  
  # Test knockoff variables
  cat("   Testing knockoff variables...\n")
  M <- length(X_knockoffs)
  knockoff_results <- vector("list", M)
  
  for (k in seq_len(M)) {
    cat("     Knockoff copy", k, "/", M, "\n")
    knockoff_results[[k]] <- fit_cox_spa(X_knockoffs[[k]], time, status, covariates, use_spa = use_spa)
  }
  
  # Combine test results
  test_results <- list(
    original = orig_results,
    knockoffs = knockoff_results
  )
  
  cat("   - Association testing complete!\n")
  
  # Step 3: Apply knockoff filter
  cat("\n3. APPLYING KNOCKOFF FILTER\n")
  cat("   Computing W statistics using method:", method, "\n")
  
  # Calculate W statistics
  t_orig <- orig_results$test_stats
  t_knock_list <- lapply(knockoff_results, function(x) x$test_stats)
  
  W_stats <- calculate_w_statistics(t_orig, t_knock_list, method = method)
  
  cat("   Applying knockoff filter with FDR =", fdr, "\n")
  selected_vars <- knockoff_filter(W_stats, fdr = fdr)
  threshold <- attr(selected_vars, "threshold")
  
  # Create filter results structure  
  filter_results <- list(
    selected_vars = selected_vars,
    W_stats = W_stats,
    threshold = threshold
  )
  
  cat("   - Variable selection complete!\n")
  
  # Summary
  cat("\n=== ANALYSIS SUMMARY ===\n")
  cat("   Total variables tested:", length(filter_results$W_stats), "\n")
  cat("   Variables selected:", length(filter_results$selected_vars), "\n")
  cat("   Selection proportion:", round(length(filter_results$selected_vars) / 
                                        length(filter_results$W_stats) * 100, 2), "%\n")
  cat("   Threshold used:", round(filter_results$threshold, 4), "\n")
  if (save_gds && !is.null(gds_output)) {
    cat("   Knockoff data saved to:", gds_output, "\n")
  }
  
  # Combine results
  results <- list(
    selected_vars = filter_results$selected_vars,
    W_stats = filter_results$W_stats,
    threshold = filter_results$threshold,
    gds_file = gds_output,
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

#' Fit Cox Regression Model
#'
#' Fits a Cox proportional hazards model with optional covariates.
#'
#' @param X Predictor matrix
#' @param time Survival times
#' @param status Event indicators
#' @param covariates Optional covariate matrix/data.frame
#' @return Fitted Cox model object
#' @export
fit_cox_model <- function(X, time, status, covariates = NULL) {
  
  surv_obj <- Surv(time, status)
  
  if (!is.null(covariates)) {
    data_df <- data.frame(X, covariates, surv_obj)
    var_names <- c(paste0("X", seq_len(ncol(X))), names(covariates))
    names(data_df)[seq_len(ncol(X))] <- var_names[seq_len(ncol(X))]
    formula_str <- paste("surv_obj ~", paste(var_names, collapse = " + "))
  } else {
    data_df <- data.frame(X, surv_obj)
    var_names <- paste0("X", seq_len(ncol(X)))
    names(data_df)[seq_len(ncol(X))] <- var_names
    formula_str <- paste("surv_obj ~", paste(var_names, collapse = " + "))
  }
  
  coxph(as.formula(formula_str), data = data_df)
}
