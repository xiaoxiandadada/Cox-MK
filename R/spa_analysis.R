# Perform Association Testing
#
# Performs association testing with genotype data using a fitted null model.
# This function handles both SPACox and standard Cox regression methods
# for genetic variant association analysis.
# This is an internal function.
#
# @param X Genotype matrix (samples × SNPs)
# @param null_model Fitted null Cox model from fit_null_cox_model()
# @return List with test statistics and p-values
# @examples
# \dontrun{
# # Fit null model first
# null_model <- fit_null_cox_model(time, status, covariates)
# 
# # Perform association testing
# results <- perform_association_testing(genotype_matrix, null_model)
# 
# # Extract results
# test_stats <- results$test_stats
# p_values <- results$p_values
# }
perform_association_testing <- function(X, null_model, time = NULL, status = NULL, covariates = NULL) {
  
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
        if (!requireNamespace("SPACox", quietly = TRUE)) {
          stop("SPACox package not available for score test")
        }
        result <- get("SPACox_Score_Test", asNamespace("SPACox"))(
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
    if (!is.null(time) && !is.null(status)) {
      base_data <- data.frame(
        time = as.numeric(time),
        status = as.numeric(status),
        check.names = FALSE
      )
    } else {
      base_data <- null_model$model
      if (is.null(base_data)) {
        base_data <- model.frame(null_model)
      }
      if (is.null(base_data)) {
        stop("Null model does not contain model frame. Provide time/status inputs.")
      }
      if (!("time" %in% names(base_data)) || !("status" %in% names(base_data))) {
        if (!is.null(null_model$y)) {
          base_data$time <- null_model$y[, 1]
          base_data$status <- null_model$y[, 2]
        } else {
          stop("Unable to locate time/status columns for association testing.")
        }
      }
      surv_col <- "Surv(time, status)"
      if (surv_col %in% names(base_data)) {
        base_data[[surv_col]] <- NULL
      }
    }
    if (!is.null(covariates)) {
      covar_df <- as.data.frame(covariates)
      covar_df <- covar_df[, setdiff(names(covar_df), c("time", "status", "snp")), drop = FALSE]
      missing_cols <- setdiff(names(covar_df), names(base_data))
      if (length(missing_cols) > 0) {
        base_data <- cbind(base_data, covar_df[, missing_cols, drop = FALSE])
      }
    }

    for (j in seq_len(n_snps)) {
      if (j %% 1000 == 0) {
        cat("    Processed", j, "/", n_snps, "SNPs\n")
      }
      
      tryCatch({
        # Add SNP to the data
        test_data <- base_data
        test_data$snp <- X[, j]
        
        # Get original formula and add SNP
        original_terms <- attr(terms(null_model), "term.labels")
        rhs_terms <- c("snp", original_terms)
        rhs_str <- if (length(rhs_terms) > 0) paste(rhs_terms, collapse = " + ") else "1"
        new_formula <- as.formula(paste("Surv(time, status) ~", rhs_str))
        
        # Fit model with SNP
        cox_fit <- survival::coxph(new_formula, data = test_data)
        
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
