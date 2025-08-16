#' Perform Association Testing
#'
#' Performs association testing with genotype data using a fitted null model.
#' This function handles both SPACox and standard Cox regression methods
#' for genetic variant association analysis.
#'
#' @param X Genotype matrix (samples × SNPs)
#' @param null_model Fitted null Cox model from fit_null_cox_model()
#' @return List with test statistics and p-values
#' @export
#' @examples
#' \dontrun{
#' # Fit null model first
#' null_model <- fit_null_cox_model(time, status, covariates)
#' 
#' # Perform association testing
#' results <- perform_association_testing(genotype_matrix, null_model)
#' 
#' # Extract results
#' test_stats <- results$test_stats
#' p_values <- results$p_values
#' }
perform_association_testing <- function(X, null_model) {
  
  if (!is.matrix(X) && !inherits(X, "Matrix")) {
    X <- as.matrix(X)
  }
  
  n_snps <- ncol(X)
  test_stats <- numeric(n_snps)
  p_values <- numeric(n_snps)
  
  # Check model type
  model_type <- attr(null_model, "model_type")
  if (is.null(model_type)) {
    # Try to infer model type
    if (inherits(null_model, "SPACox_Null_Model")) {
      model_type <- "SPACox"
    } else {
      model_type <- "Standard Cox"
    }
  }
  
  if (model_type == "SPACox" && requireNamespace("SPACox", quietly = TRUE)) {
    # Use SPACox for testing
    tryCatch({
      
      for (j in seq_len(n_snps)) {
        if (j %% 1000 == 0) {
          cat("    Processed", j, "/", n_snps, "SNPs\n")
        }
        
        # SPACox score test
        result <- SPACox::SPACox_Score_Test(
          geno = X[, j],
          obj_nullmodel = null_model
        )
        
        test_stats[j] <- result$Score^2 / result$Var
        p_values[j] <- result$p.value
      }
      
      method_used <- "SPACox"
      
    }, error = function(e) {
      warning("SPACox testing failed: ", e$message, "\nFalling back to standard method")
      # Fallback handled below
      model_type <- "Standard Cox"
    })
  }
  
  if (model_type == "Standard Cox") {
    # Standard Cox regression testing
    for (j in seq_len(n_snps)) {
      if (j %% 1000 == 0) {
        cat("    Processed", j, "/", n_snps, "SNPs\n")
      }
      
      tryCatch({
        # Extract original data from null model
        null_data <- null_model$model
        
        # Add SNP to the data
        test_data <- null_data
        test_data$snp <- X[, j]
        
        # Get original formula and add SNP
        original_terms <- attr(terms(null_model), "term.labels")
        if (length(original_terms) > 0) {
          new_formula <- paste("Surv(time, status) ~ snp +", paste(original_terms, collapse = " + "))
        } else {
          new_formula <- "Surv(time, status) ~ snp"
        }
        
        # Fit model with SNP
        cox_fit <- survival::coxph(as.formula(new_formula), data = test_data)
        
        # Extract test statistics
        coef_summary <- summary(cox_fit)$coefficients
        
        if ("snp" %in% rownames(coef_summary)) {
          test_stats[j] <- coef_summary["snp", "z"]^2
          p_values[j] <- coef_summary["snp", "Pr(>|z|)"]
        } else {
          test_stats[j] <- 0
          p_values[j] <- 1
        }
        
      }, error = function(e) {
        test_stats[j] <- 0
        p_values[j] <- 1
      })
    }
    
    method_used <- "Standard Cox"
  }
  
  return(list(
    test_stats = test_stats,
    p_values = p_values,
    method = method_used,
    n_snps = n_snps
  ))
}

#' Calculate W Statistics
#'
#' Calculate W statistics for knockoff filter based on test statistics
#' from original and knockoff variables.
#'
#' @param Z_orig Test statistics for original variables
#' @param Z_ko Test statistics for knockoff variables  
#' @param method Method for combining statistics ("median" or "difference")
#' @return Vector of W statistics
#' @noRd
calculate_w_statistics <- function(Z_orig, Z_ko, method = "median") {
  
  if (length(Z_orig) != length(Z_ko)) {
    stop("Original and knockoff statistics must have same length")
  }
  
  if (method == "median") {
    W <- pmax(Z_orig, Z_ko) * sign(Z_orig - Z_ko)
  } else if (method == "difference") {
    W <- Z_orig - Z_ko
  } else {
    stop("Method must be 'median' or 'difference'")
  }
  
  return(W)
}

#' Apply SPA Analysis Workflow
#'
#' Complete SPA analysis workflow that performs association testing
#' on both original and knockoff variables, then calculates W statistics.
#'
#' @param X_orig Original genotype matrix
#' @param X_ko Knockoff genotype matrix
#' @param null_model Fitted null Cox model
#' @param method Method for W statistics ("median" or "difference")
#' @return List with W statistics and intermediate results
#' @noRd
spa_analysis_workflow <- function(X_orig, X_ko, null_model, method = "median") {
  
  cat("Starting SPA analysis workflow...\n")
  
  # Test original variables
  cat("Testing original variables...\n")
  orig_results <- perform_association_testing(X_orig, null_model)
  
  # Test knockoff variables  
  cat("Testing knockoff variables...\n")
  ko_results <- perform_association_testing(X_ko, null_model)
  
  # Calculate W statistics
  cat("Calculating W statistics...\n")
  W_stats <- calculate_w_statistics(
    orig_results$test_stats, 
    ko_results$test_stats, 
    method = method
  )
  
  cat("SPA analysis completed!\n")
  
  return(list(
    W_stats = W_stats,
    orig_results = orig_results,
    ko_results = ko_results,
    method = method
  ))
}
