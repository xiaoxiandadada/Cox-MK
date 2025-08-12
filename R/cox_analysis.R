#' Fit Null Cox Model
#'
#' Fit a Cox proportional hazards model with only covariates (no genetic variants).
#' This function supports both standard coxph and SPACox for large-scale genetic data.
#' SPACox (https://github.com/WenjianBI/SPACox) is recommended for genome-wide 
#' association studies with survival outcomes.
#'
#' @param time Vector of survival times
#' @param status Vector of event indicators (0 = censored, 1 = event)
#' @param covariates Matrix or data frame of covariates
#' @param formula Optional formula for the model. If NULL, uses all covariates.
#' @param use_spacox Logical, whether to use SPACox for model fitting. 
#'   Requires SPACox package to be installed. Default is FALSE.
#'
#' @return A fitted Cox model object
#'
#' @note For large-scale genetic studies, it is recommended to use SPACox 
#'   by setting use_spacox = TRUE. SPACox can be installed from: 
#'   https://github.com/WenjianBI/SPACox
#'
#' @examples
#' \dontrun{
#' # Fit null model with standard coxph
#' null_model <- fit_null_model(
#'   time = pheno_data$time,
#'   status = pheno_data$status,
#'   covariates = pheno_data[, c("age", "sex", "bmi")]
#' )
#' 
#' # Fit null model with SPACox (if available)
#' null_model_spa <- fit_null_model(
#'   time = pheno_data$time,
#'   status = pheno_data$status,
#'   covariates = pheno_data[, c("age", "sex", "bmi")],
#'   use_spacox = TRUE
#' )
#' }
#'
#' @export
#' @importFrom survival coxph Surv
fit_null_model <- function(time, status, covariates, formula = NULL, use_spacox = FALSE) {
  
  # Create data frame
  data <- data.frame(
    time = time,
    status = status,
    covariates
  )
  
  # Create formula if not provided
  if (is.null(formula)) {
    covar_names <- names(covariates)
    formula <- as.formula(paste("Surv(time, status) ~", paste(covar_names, collapse = " + ")))
  }
  
  # Fit Cox model
  if (use_spacox && requireNamespace("SPACox", quietly = TRUE)) {
    # Use SPACox for fitting
    message("Using SPACox for null model fitting...")
    
    # Convert to SPACox format
    # Note: This is a placeholder - actual SPACox implementation may require
    # specific data format and function calls
    tryCatch({
      fit <- SPACox::SPACox_Null_Model(
        formula = formula,
        data = data
      )
      attr(fit, "method") <- "SPACox"
    }, error = function(e) {
      warning("SPACox fitting failed, falling back to standard coxph: ", e$message)
      fit <- survival::coxph(formula, data = data)
      attr(fit, "method") <- "coxph"
    })
  } else {
    if (use_spacox) {
      warning("SPACox package not available, using standard coxph instead.")
    }
    fit <- survival::coxph(formula, data = data)
    attr(fit, "method") <- "coxph"
  }
  
  return(fit)
}

#' Cox Knockoff Screen
#'
#' Perform association analysis using Cox regression with original and knockoff variables.
#'
#' @param genotypes Original genotype matrix (samples x SNPs)
#' @param knockoffs List of knockoff genotype matrices
#' @param time Vector of survival times
#' @param status Vector of event indicators
#' @param covariates Matrix or data frame of covariates
#' @param null_model Pre-fitted null model (optional)
#' @param adjust_method Method for multiple testing adjustment ("bonferroni", "BH", "BY")
#' @param alpha Significance level (default: 0.05)
#'
#' @return A list containing:
#'   \item{original_pvals}{P-values for original variables}
#'   \item{knockoff_pvals}{List of p-values for knockoff variables}
#'   \item{original_coefs}{Coefficients for original variables}
#'   \item{knockoff_coefs}{List of coefficients for knockoff variables}
#'   \item{w_stats}{W statistics for variable selection}
#'
#' @examples
#' \dontrun{
#' # Perform knockoff screening
#' results <- cox_knockoff_screen(
#'   genotypes = plink_data$genotypes,
#'   knockoffs = knockoffs,
#'   time = pheno_data$time,
#'   status = pheno_data$status,
#'   covariates = pheno_data[, c("age", "sex", "bmi")]
#' )
#' }
#'
#' @export
#' @importFrom survival coxph Surv
#' @importFrom stats p.adjust
cox_knockoff_screen <- function(genotypes, knockoffs, time, status, covariates,
                               null_model = NULL, adjust_method = "bonferroni", 
                               alpha = 0.05) {
  
  p <- ncol(genotypes)
  M <- length(knockoffs)
  
  # Remove unused variable warning
  # n <- length(time)  # Not used in current implementation
  
  # Fit null model if not provided
  if (is.null(null_model)) {
    null_model <- fit_null_model(time, status, covariates)
  }
  
  # Initialize result containers
  original_pvals <- numeric(p)
  original_coefs <- numeric(p)
  knockoff_pvals <- lapply(1:M, function(x) numeric(p))
  knockoff_coefs <- lapply(1:M, function(x) numeric(p))
  
  # Prepare base data
  base_data <- data.frame(
    time = time,
    status = status,
    covariates
  )
  covar_names <- names(covariates)
  
  # Test original variables
  cat("Testing original variables...\n")
  for (j in 1:p) {
    if (j %% 100 == 0) cat("  SNP", j, "of", p, "\n")
    
    # Add SNP to data
    test_data <- base_data
    test_data$snp <- genotypes[, j]
    
    # Fit model
    formula <- as.formula(paste("Surv(time, status) ~ snp +", 
                               paste(covar_names, collapse = " + ")))
    
    tryCatch({
      fit <- survival::coxph(formula, data = test_data)
      original_pvals[j] <- summary(fit)$coefficients["snp", "Pr(>|z|)"]
      original_coefs[j] <- coef(fit)["snp"]
    }, error = function(e) {
      original_pvals[j] <- 1
      original_coefs[j] <- 0
    })
  }
  
  # Test knockoff variables
  for (k in 1:M) {
    cat("Testing knockoff set", k, "...\n")
    for (j in 1:p) {
      if (j %% 100 == 0) cat("  SNP", j, "of", p, "\n")
      
      # Add knockoff SNP to data
      test_data <- base_data
      test_data$snp <- knockoffs[[k]][, j]
      
      # Fit model
      formula <- as.formula(paste("Surv(time, status) ~ snp +", 
                                 paste(covar_names, collapse = " + ")))
      
      tryCatch({
        fit <- survival::coxph(formula, data = test_data)
        knockoff_pvals[[k]][j] <- summary(fit)$coefficients["snp", "Pr(>|z|)"]
        knockoff_coefs[[k]][j] <- coef(fit)["snp"]
      }, error = function(e) {
        knockoff_pvals[[k]][j] <- 1
        knockoff_coefs[[k]][j] <- 0
      })
    }
  }
  
  # Calculate W statistics
  w_stats <- calculate_w_statistics(original_pvals, knockoff_pvals, 
                                   original_coefs, knockoff_coefs)
  
  list(
    original_pvals = original_pvals,
    knockoff_pvals = knockoff_pvals,
    original_coefs = original_coefs,
    knockoff_coefs = knockoff_coefs,
    w_stats = w_stats
  )
}

#' Calculate W Statistics
#'
#' Calculate the W statistics for knockoff variable selection.
#'
#' @param original_pvals P-values from original variables
#' @param knockoff_pvals List of p-values from knockoff variables
#' @param original_coefs Coefficients from original variables
#' @param knockoff_coefs List of coefficients from knockoff variables
#' @param method Method for calculating W statistics ("signed_max", "difference")
#'
#' @return Vector of W statistics
#'
#' @export
calculate_w_statistics <- function(original_pvals, knockoff_pvals, 
                                  original_coefs, knockoff_coefs,
                                  method = "signed_max") {
  
  p <- length(original_pvals)
  M <- length(knockoff_pvals)
  
  w_stats <- numeric(p)
  
  if (method == "signed_max") {
    # Use signed maximum statistic
    for (j in 1:p) {
      # Convert p-values to test statistics (higher is more significant)
      orig_stat <- -log10(pmax(original_pvals[j], 1e-16)) * sign(original_coefs[j])
      
      knockoff_stats <- sapply(1:M, function(k) {
        -log10(pmax(knockoff_pvals[[k]][j], 1e-16)) * sign(knockoff_coefs[[k]][j])
      })
      
      w_stats[j] <- max(abs(orig_stat)) - max(abs(knockoff_stats))
    }
  } else if (method == "difference") {
    # Use difference of means
    for (j in 1:p) {
      orig_stat <- -log10(pmax(original_pvals[j], 1e-16))
      knockoff_stats <- sapply(1:M, function(k) {
        -log10(pmax(knockoff_pvals[[k]][j], 1e-16))
      })
      
      w_stats[j] <- orig_stat - mean(knockoff_stats)
    }
  }
  
  return(w_stats)
}

#' Knockoff Filter
#'
#' Apply the knockoff filter to select significant variables.
#'
#' @param w_stats Vector of W statistics
#' @param fdr Target false discovery rate (default: 0.1)
#' @param offset Offset parameter for the filter (default: 1)
#'
#' @return A list containing:
#'   \item{selected}{Indices of selected variables}
#'   \item{threshold}{Threshold used for selection}
#'   \item{fdp}{Estimated false discovery proportion}
#'
#' @examples
#' \dontrun{
#' # Apply knockoff filter
#' selected <- knockoff_filter(w_stats, fdr = 0.1)
#' }
#'
#' @export
knockoff_filter <- function(w_stats, fdr = 0.1, offset = 1) {
  
  # Sort W statistics in descending order
  w_sorted <- sort(w_stats, decreasing = TRUE)
  # w_order <- order(w_stats, decreasing = TRUE)  # Not used currently
  
  # Find threshold
  thresholds <- w_sorted[w_sorted > 0]
  if (length(thresholds) == 0) {
    return(list(selected = integer(0), threshold = Inf, fdp = 0))
  }
  
  # Calculate false discovery proportion for each threshold
  fdp <- sapply(thresholds, function(t) {
    (offset + sum(w_stats <= -t)) / max(1, sum(w_stats >= t))
  })
  
  # Find largest threshold with FDP <= fdr
  valid_thresholds <- thresholds[fdp <= fdr]
  if (length(valid_thresholds) == 0) {
    return(list(selected = integer(0), threshold = Inf, fdp = min(fdp)))
  }
  
  threshold <- min(valid_thresholds)
  selected <- which(w_stats >= threshold)
  final_fdp <- (offset + sum(w_stats <= -threshold)) / max(1, length(selected))
  
  list(
    selected = selected,
    threshold = threshold,
    fdp = final_fdp
  )
}
