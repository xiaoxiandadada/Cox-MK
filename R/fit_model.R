#Fit Cox Model from Files
#'
#Fits a null Cox model by reading phenotype and covariate data from files.
#This function is designed for batch processing and large-scale analysis
#where data is stored in separate files.
#'
#@param phenotype_file Path to CSV file with columns: IID, time, status
#@param covariate_file Path to CSV file with columns: IID, covar1, covar2, ...
#@param output_file Path to RDS file to save the fitted null model
#@param use_spacox Whether to try using SPACox package (default: TRUE)
#@return Invisible path to the output file
#@keywords internal
#@examples
#\dontrun{
## Prepare example data files
#pheno_data <- data.frame(
#  IID = paste0("ID", 1:100),
#  time = rexp(100, 0.1),
#  status = rbinom(100, 1, 0.3)
#)
#covar_data <- data.frame(
#  IID = paste0("ID", 1:100),
#  age = rnorm(100, 50, 10),
#  sex = rbinom(100, 1, 0.5)
#)
#
#write.csv(pheno_data, "phenotype.csv", row.names = FALSE)
#write.csv(covar_data, "covariates.csv", row.names = FALSE)
#
## Fit model from files
#fit_cox_model_from_files(
#  phenotype_file = "phenotype.csv",
#  covariate_file = "covariates.csv", 
#  output_file = "null_model.rds"
#)
#
## Load the fitted model
#model_info <- readRDS("null_model.rds")
#}
fit_cox_model_from_files <- function(phenotype_file, covariate_file, output_file, 
                                   use_spacox = TRUE) {
  
  #
  if (!file.exists(phenotype_file)) {
    stop("Phenotype file not found: ", phenotype_file)
  }
  
  if (!file.exists(covariate_file)) {
    stop("Covariate file not found: ", covariate_file)
  }
  
  #
  cat("Loading phenotype data from:", phenotype_file, "\n")
  pheno_data <- read.csv(phenotype_file, stringsAsFactors = FALSE)
  
  # Check required columns
  required_cols <- c("IID", "time", "status")
  missing_cols <- setdiff(required_cols, colnames(pheno_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in phenotype file: ", paste(missing_cols, collapse = ", "))
  }
  
  cat("  - Loaded", nrow(pheno_data), "samples\n")
  cat("  - Events:", sum(pheno_data$status), "/", nrow(pheno_data), 
      "(", round(100 * sum(pheno_data$status) / nrow(pheno_data), 1), "%)\n")
  
  #
  cat("Loading covariate data from:", covariate_file, "\n")
  covar_data <- read.csv(covariate_file, stringsAsFactors = FALSE)
  
  # Check IID column exists
  if (!"IID" %in% colnames(covar_data)) {
    stop("IID column not found in covariate file")
  }
  
  cat("  - Loaded", nrow(covar_data), "samples\n")
  cat("  - Covariates:", paste(setdiff(colnames(covar_data), "IID"), collapse = ", "), "\n")
  
  #
  cat("Merging phenotype and covariate data...\n")
  merged_data <- merge(pheno_data, covar_data, by = "IID", all.x = TRUE)
  
  # Check for missing data after merge
  n_complete <- sum(complete.cases(merged_data))
  if (n_complete < nrow(merged_data)) {
    cat("  Warning:", nrow(merged_data) - n_complete, "samples have missing covariate data\n")
    merged_data <- merged_data[complete.cases(merged_data), ]
    cat("  Using", nrow(merged_data), "complete samples\n")
  }
  
  #
  cat("Fitting null Cox model...\n")
  
  # Fit the model
  if (use_spacox && requireNamespace("SPACox", quietly = TRUE)) {
    cat("  Using SPACox for null model fitting...\n")
    
    # Prepare formula
    covar_names <- setdiff(colnames(merged_data), c("IID", "time", "status"))
    if (length(covar_names) > 0) {
      formula_str <- paste("Surv(time, status) ~", paste(covar_names, collapse = " + "))
    } else {
      formula_str <- "Surv(time, status) ~ 1"
    }
    
    cat("  Formula:", formula_str, "\n")
    
    # Fit SPACox null model
    tryCatch({
      if (!requireNamespace("SPACox", quietly = TRUE)) {
        stop("SPACox package not available")
      }
      null_model <- SPACox::SPACox_Null_Model(
        formula = as.formula(formula_str),
        data = merged_data,
        pIDs = merged_data$IID,
        gIDs = merged_data$IID  # For unrelated samples
      )
      
      cat("  - SPACox null model fitted successfully!\n")
      model_type <- "SPACox"
      
    }, error = function(e) {
      cat("  SPACox fitting failed:", e$message, "\n")
      cat("  Falling back to standard Cox regression...\n")
      
      # Fallback to standard Cox model
      null_model <- coxph(as.formula(formula_str), data = merged_data)
      model_type <- "Standard Cox"
    })
    
  } else {
    if (use_spacox) {
      cat("  SPACox not available, using standard Cox regression...\n")
    } else {
      cat("  Using standard Cox regression...\n")
    }
    
    # Standard Cox model
    covar_names <- setdiff(colnames(merged_data), c("IID", "time", "status"))
    if (length(covar_names) > 0) {
      formula_str <- paste("Surv(time, status) ~", paste(covar_names, collapse = " + "))
    } else {
      formula_str <- "Surv(time, status) ~ 1"
    }
    
    null_model <- coxph(as.formula(formula_str), data = merged_data)
    model_type <- "Standard Cox"
  }
  
  #
  cat("Saving fitted model to:", output_file, "\n")
  
  # Add metadata to the model
  model_info <- list(
    model = null_model,
    model_type = model_type,
    formula = formula_str,
    n_samples = nrow(merged_data),
    n_events = sum(merged_data$status),
    n_covariates = length(setdiff(colnames(merged_data), c("IID", "time", "status"))),
    covariate_names = setdiff(colnames(merged_data), c("IID", "time", "status")),
    fit_time = Sys.time(),
    r_version = R.version.string
  )
  
  # Save to RDS file
  saveRDS(model_info, file = output_file)
  
  #
  cat("\n=== Model Fitting Summary ===\n")
  cat("Model type:", model_type, "\n")
  cat("Formula:", formula_str, "\n")
  cat("Samples used:", nrow(merged_data), "\n")
  cat("Events:", sum(merged_data$status), "\n")
  cat("Covariates:", length(model_info$covariate_names), "\n")
  if (length(model_info$covariate_names) > 0) {
    cat("  -", paste(model_info$covariate_names, collapse = ", "), "\n")
  }
  cat("Output file:", output_file, "\n")
  cat("Fit completed at:", format(Sys.time()), "\n")
  
  # Print model summary if available
  if (model_type == "Standard Cox" && inherits(null_model, "coxph")) {
    cat("\nModel Summary:\n")
    print(summary(null_model))
  }
  
  cat("\nModel fitting completed successfully!\n")
  
  invisible(output_file)
}

# Fit Null Cox Model
#
# Fits a null Cox model with only covariates (no genetic variants).
# This function uses SPACox when available for large-scale analysis,
# and falls back to standard Cox regression otherwise.
# This is an internal function.
#
# @param time Survival times
# @param status Event indicators
# @param covariates Covariate data frame (optional)
# @return Fitted Cox model object (SPACox or coxph)
# @examples
# \dontrun{
# # Example with covariates
# data(example_phenotype)
# data(example_covariates)
# 
# null_model <- fit_null_cox_model(
#   time = example_phenotype$time,
#   status = example_phenotype$status,
#   covariates = example_covariates
# )
# 
# # Example without covariates
# null_model <- fit_null_cox_model(
#   time = example_phenotype$time,
#   status = example_phenotype$status
# )
# }
fit_null_cox_model <- function(time, status, covariates = NULL) {
  
  # Input validation
  if (length(time) != length(status)) {
    stop("time and status must have the same length")
  }
  
  if (sum(status) < 2) {
    stop("Need at least 2 events to fit Cox model")
  }
  
  # Prepare data
  if (!is.null(covariates)) {
    if (nrow(covariates) != length(time)) {
      stop("Number of rows in covariates must match length of time/status")
    }
    
    # Create merged data
    model_data <- data.frame(
      time = time,
      status = status,
      covariates
    )
    
    # Remove incomplete cases
    model_data <- model_data[complete.cases(model_data), ]
    
    if (nrow(model_data) < 10) {
      stop("Too few complete cases for model fitting")
    }
    
    covar_names <- colnames(covariates)
    formula_str <- paste("Surv(time, status) ~", paste(covar_names, collapse = " + "))
    
  } else {
    # No covariates
    model_data <- data.frame(time = time, status = status)
    formula_str <- "Surv(time, status) ~ 1"
  }
  
  # Try SPACox first if available
  if (requireNamespace("SPACox", quietly = TRUE)) {
    tryCatch({
      
      if (!requireNamespace("SPACox", quietly = TRUE)) {
        stop("SPACox package not available")
      }
      
      # Add sample IDs for SPACox
      model_data$sample_id <- paste0("S", seq_len(nrow(model_data)))
      
      null_model <- SPACox::SPACox_Null_Model(
        formula = as.formula(formula_str),
        data = model_data,
        pIDs = model_data$sample_id,
        gIDs = model_data$sample_id  # For unrelated samples
      )
      
      # Add metadata
      attr(null_model, "model_type") <- "SPACox"
      attr(null_model, "formula") <- formula_str
      attr(null_model, "n_samples") <- nrow(model_data)
      attr(null_model, "n_events") <- sum(model_data$status)
      
      return(null_model)
      
    }, error = function(e) {
      warning("SPACox fitting failed: ", e$message, "\nFalling back to standard Cox regression")
    })
  }
  
  # Fallback to standard Cox regression
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("survival package is required")
  }
  
  null_model <- survival::coxph(as.formula(formula_str), data = model_data)
  
  # Add metadata
  attr(null_model, "model_type") <- "Standard Cox"
  attr(null_model, "formula") <- formula_str
  attr(null_model, "n_samples") <- nrow(model_data)
  attr(null_model, "n_events") <- sum(model_data$status)
  
  return(null_model)
}
