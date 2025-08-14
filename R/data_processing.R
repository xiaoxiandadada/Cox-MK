#' @title Data Processing Utilities  
#' @name data-processing
#' @description
#' Utility functions for data loading, processing, and validation specifically
#' designed for genetic association studies.
#' @importFrom utils read.table flush.console
#' @importFrom stats sd quantile residuals
#' @importFrom gdsfmt openfn.gds read.gdsn closefn.gds index.gdsn
NULL

#' Load PLINK Format Genetic Data
#'
#' Loads genetic data from PLINK binary format (.bed/.bim/.fam files) and
#' converts it to a sparse matrix format suitable for knockoff analysis.
#'
#' @param bed_path Path to PLINK .bed file (without extension)
#' @param verbose Whether to print progress messages (default: TRUE)
#' @return List containing:
#'   \item{genotypes}{Sparse genotype matrix (samples × SNPs)}
#'   \item{fam}{Sample information from .fam file}
#'   \item{bim}{SNP information from .bim file}
#'   \item{positions}{SNP positions}
#'   \item{sample_ids}{Sample IDs}
#' @export
#' @examples
#' \dontrun{
#' # Load PLINK data
#' extdata_path <- system.file('extdata', package = 'CoxMK')
#' plink_data <- load_plink_data(file.path(extdata_path, 'sample'))
#' }
load_plink_data <- function(bed_path, verbose = TRUE) {
  
  if (verbose) cat("Loading PLINK data from:", bed_path, "\n")
  
  # Check if all required files exist
  files_needed <- paste0(bed_path, c(".bed", ".bim", ".fam"))
  missing_files <- !file.exists(files_needed)
  
  if (any(missing_files)) {
    stop("Missing PLINK files: ", paste(files_needed[missing_files], collapse = ", "))
  }
  
  # Load .fam file (sample information)
  if (verbose) cat("  Reading .fam file...\n")
  fam <- read.table(paste0(bed_path, ".fam"), stringsAsFactors = FALSE,
                    col.names = c("FID", "IID", "PID", "MID", "SEX", "PHENOTYPE"))
  
  # Load .bim file (SNP information)
  if (verbose) cat("  Reading .bim file...\n")
  bim <- read.table(paste0(bed_path, ".bim"), stringsAsFactors = FALSE,
                    col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"))
  
  # Load .bed file (genotype data)
  if (verbose) cat("  Reading .bed file...\n")
  
  if (!requireNamespace("BEDMatrix", quietly = TRUE)) {
    stop("BEDMatrix package is required for reading PLINK .bed files. Please install it with: install.packages('BEDMatrix')")
  }
  
  if (verbose) cat("    Using BEDMatrix for efficient loading...\n")
  genotypes <- Matrix(as.matrix(BEDMatrix::BEDMatrix(paste0(bed_path, ".bed"))), sparse = TRUE)
  
  if (verbose) {
    cat("  Loaded", nrow(fam), "samples and", nrow(bim), "SNPs\n")
    cat("  Genotype matrix:", class(genotypes), "with dimensions", paste(dim(genotypes), collapse = " x "), "\n")
  }
  
  list(
    genotypes = genotypes,
    fam = fam,
    bim = bim,
    positions = bim$BP,
    sample_ids = fam$IID
  )
}

#' Prepare Phenotype Data for Analysis
#'
#' Loads and processes phenotype data for survival analysis, handling various
#' input formats and performing basic validation.
#'
#' @param phenotype_file Path to phenotype file
#' @param time_col Column name or index for survival times (default: "time")
#' @param status_col Column name or index for event status (default: "status")
#' @param status_file Optional separate file for event status
#' @param sep Field separator (default: tab-delimited)
#' @param header Whether file has header (default: TRUE)
#' @return Data frame with columns: time, status, and optional sample IDs
#' @export
#' @examples
#' \dontrun{
#' # Load phenotype data
#' extdata_path <- system.file('extdata', package = 'CoxMK')
#' pheno_data <- prepare_phenotype(file.path(extdata_path, 'tte_phenotype.txt'))
#' }
prepare_phenotype <- function(phenotype_file, time_col = "time", status_col = "status",
                            status_file = NULL, sep = "\t", header = TRUE) {
  
  if (!file.exists(phenotype_file)) {
    stop("Phenotype file not found: ", phenotype_file)
  }
  
  cat("Loading phenotype data from:", phenotype_file, "\n")
  
  # Load phenotype data
  pheno_data <- read.table(phenotype_file, sep = sep, header = header, 
                          stringsAsFactors = FALSE, check.names = FALSE)
  
  cat("  Loaded", nrow(pheno_data), "samples\n")
  
  # Load separate status file if provided
  if (!is.null(status_file)) {
    if (!file.exists(status_file)) {
      stop("Status file not found: ", status_file)
    }
    
    status_data <- read.table(status_file, sep = sep, header = header, stringsAsFactors = FALSE)
    
    if (nrow(pheno_data) != nrow(status_data)) {
      stop("Number of samples in phenotype and status files must match")
    }
    
    # Combine data
    pheno_data <- cbind(pheno_data, status_data)
  }
  
  # Extract time and status columns
  if (is.character(time_col)) {
    if (!time_col %in% names(pheno_data)) {
      stop("Time column '", time_col, "' not found in phenotype data")
    }
    time_values <- pheno_data[[time_col]]
  } else {
    time_values <- pheno_data[, time_col]
  }
  
  if (is.character(status_col)) {
    if (!status_col %in% names(pheno_data)) {
      stop("Status column '", status_col, "' not found in phenotype data")
    }
    status_values <- pheno_data[[status_col]]
  } else {
    status_values <- pheno_data[, status_col]
  }
  
  # Create result data frame
  result_data <- data.frame(
    time = as.numeric(time_values),
    status = as.numeric(status_values),
    stringsAsFactors = FALSE
  )
  
  # Add sample IDs if available
  if ("FID" %in% names(pheno_data) && "IID" %in% names(pheno_data)) {
    result_data$sample_id <- pheno_data$IID
  } else if ("sample_id" %in% names(pheno_data)) {
    result_data$sample_id <- pheno_data$sample_id
  }
  
  # Validation
  if (any(is.na(result_data$time)) || any(is.na(result_data$status))) {
    stop("Missing values detected in time or status columns")
  }
  
  if (any(result_data$time <= 0)) {
    stop("All survival times must be positive")
  }
  
  if (!all(result_data$status %in% c(0, 1))) {
    stop("Status values must be 0 (censored) or 1 (event)")
  }
  
  cat("  Validation passed:", nrow(result_data), "complete cases\n")
  cat("  Events:", sum(result_data$status), "(",
      round(mean(result_data$status) * 100, 1), "%)\n")
  
  return(result_data)
}

#' Load Covariate Data
#'
#' Loads and processes covariate data for inclusion in survival models.
#'
#' @param covariate_file Path to covariate file
#' @param exclude_cols Column names to exclude from covariate matrix, typically
#'   the analysis (default: c("FID", "IID")).
#' @return Data frame of covariates
#' @export
load_covariates <- function(covariate_file, exclude_cols = c("FID", "IID")) {
  
  if (!file.exists(covariate_file)) {
    stop("Covariate file not found: ", covariate_file)
  }
  
  cat("Loading covariate data from:", covariate_file, "\n")
  
  # Load covariate data
  covar_data <- read.table(covariate_file, header = TRUE, stringsAsFactors = FALSE,
                          check.names = FALSE)
  
  cat("  Loaded", nrow(covar_data), "samples with", ncol(covar_data), "variables\n")
  
  # Remove excluded columns
  if (!is.null(exclude_cols)) {
    keep_cols <- !names(covar_data) %in% exclude_cols
    if (sum(keep_cols) == 0) {
      stop("All columns excluded. Check exclude_cols parameter.")
    }
    covar_data <- covar_data[, keep_cols, drop = FALSE]
  }
  
  # Convert character variables to factors if appropriate
  for (col in names(covar_data)) {
    if (is.character(covar_data[[col]])) {
      # Check if it looks like a factor (limited unique values)
      unique_vals <- length(unique(covar_data[[col]]))
      if (unique_vals <= 10 && unique_vals < nrow(covar_data) / 10) {
        covar_data[[col]] <- as.factor(covar_data[[col]])
        cat("    Converted", col, "to factor\n")
      }
    }
  }
  
  # Check for missing values
  missing_counts <- sapply(covar_data, function(x) sum(is.na(x)))
  if (any(missing_counts > 0)) {
    cat("  Missing values detected:\n")
    for (col in names(missing_counts)[missing_counts > 0]) {
      cat("    ", col, ":", missing_counts[col], "missing\n")
    }
  }
  
  cat("  Final covariate matrix:", ncol(covar_data), "variables\n")
  
  return(covar_data)
}

#' Load Knockoff Data from GDS File
#'
#' Load original genotypes and knockoff variables from a GDS file.
#'
#' @param gds_file Path to the GDS file containing knockoff data
#' @return List containing:
#'   \item{original}{Original genotype matrix}
#'   \item{knockoffs}{List of knockoff matrices}
#'   \item{sample_ids}{Sample IDs}
#'   \item{positions}{SNP positions}
#' @export
load_knockoff_gds <- function(gds_file) {
  
  if (!file.exists(gds_file)) {
    stop("GDS file does not exist: ", gds_file)
  }
  
  if (!requireNamespace("gdsfmt", quietly = TRUE)) {
    stop("gdsfmt package is required to read GDS files")
  }
  
  cat("Loading knockoff data from GDS file:", gds_file, "\n")
  
  # Open GDS file
  g <- gdsfmt::openfn.gds(gds_file)
  
  tryCatch({
    # Read sample IDs
    sample_ids <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(g, "sample.id"))
    
    # Read SNP positions
    positions <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(g, "snp.pos"))
    
    # Read original genotype data
    original <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(g, "original"))
    
    # Read knockoff data - they are saved as separate knockoff1, knockoff2, etc.
    # First, determine how many knockoffs exist
    all_nodes <- gdsfmt::ls.gdsn(g)
    knockoff_nodes <- grep("^knockoff[0-9]+$", all_nodes, value = TRUE)
    M <- length(knockoff_nodes)
    
    # Read each knockoff matrix
    knockoffs <- vector("list", M)
    for (k in seq_len(M)) {
      node_name <- paste0("knockoff", k)
      knockoffs[[k]] <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(g, node_name))
    }
    
    cat("  Loaded data:\n")
    cat("    Samples:", length(sample_ids), "\n")
    cat("    SNPs:", length(positions), "\n") 
    cat("    Knockoff copies:", length(knockoffs), "\n")
    
    return(list(
      original = original,
      knockoffs = knockoffs,
      sample_ids = sample_ids,
      positions = positions,
      gds_file = gds_file
    ))
    
  }, finally = {
    gdsfmt::closefn.gds(g)
  })
}