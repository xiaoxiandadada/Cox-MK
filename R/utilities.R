#' Sparse Matrix Correlation Computation
#'
#' Efficiently computes correlation and covariance matrices for sparse matrices
#' using optimized sparse matrix operations. This is an internal utility function.
#'
#' @param x A sparse matrix (can be dgCMatrix or similar)
#' @return List containing:
#'   \item{cov}{The covariance matrix}
#'   \item{cor}{The correlation matrix}
#' @keywords internal
sparse_cor <- function(x) {
  # Ensure x is a matrix-like object
  if (!is.matrix(x) && !inherits(x, "Matrix")) {
    x <- as.matrix(x)
  }
  
  # Compute covariance and correlation matrices
  cov_matrix <- var(x, na.rm = TRUE)
  cor_matrix <- stats::cov2cor(cov_matrix)
  
  return(list(
    cov = cov_matrix,
    cor = cor_matrix
  ))
}

#' Sparse Cross-Covariance Computation
#'
#' Computes cross-covariance between two sparse matrices efficiently.
#' This is an internal utility function for knockoff generation.
#'
#' @param x A numeric matrix (can be sparse)
#' @param y A numeric matrix (can be sparse)
#' @return List containing:
#'   \item{cov}{The cross-covariance matrix}
#' @keywords internal
sparse_cov_cross <- function(x, y) {
  # Compute cross-covariance matrix
  n <- nrow(x)
  cMeans_x <- colMeans(x)
  cMeans_y <- colMeans(y)
  
  # Center the matrices
  x_centered <- sweep(x, 2, cMeans_x, "-")
  y_centered <- sweep(y, 2, cMeans_y, "-")
  
  # Compute cross-covariance
  cross_cov <- crossprod(x_centered, y_centered) / (n - 1)
  
  return(list(cov = cross_cov))
}

#' Find N-th Largest Value
#'
#' Efficiently finds the n-th largest value in a numeric vector.
#' This is primarily used internally for knockoff statistics computation.
#'
#' @param x A numeric vector
#' @param n The position (1 = largest, 2 = second largest, etc.)
#' @return The n-th largest value in the vector
#' @keywords internal
max_nth <- function(x, n) {
  #
  if (n > length(x)) {
    stop("n cannot be larger than length of x")
  }
  return(sort(x, partial = length(x) - (n - 1))[length(x) - (n - 1)])
}

#' Safe Matrix Conversion
#'
#' Safely converts various matrix-like objects to standard matrix format,
#' handling sparse matrices and potential memory issues. This is an internal
#' utility function.
#'
#' @param x Matrix-like object (matrix, Matrix, data.frame)
#' @param sparse Whether to maintain sparsity if possible (default: FALSE)
#' @return Matrix object
#' @keywords internal
safe_as_matrix <- function(x, sparse = FALSE) {
  
  #
  if (inherits(x, "Matrix")) {
    if (sparse) {
      return(x)
    } else {
      return(as.matrix(x))
    }
  } else if (is.matrix(x)) {
    if (sparse) {
      return(Matrix::Matrix(x, sparse = TRUE))
    } else {
      return(x)
    }
  } else if (is.data.frame(x)) {
    x_matrix <- as.matrix(x)
    if (sparse) {
      return(Matrix::Matrix(x_matrix, sparse = TRUE))
    } else {
      return(x_matrix)
    }
  } else {
    # Try to convert to matrix
    x_matrix <- as.matrix(x)
    if (sparse) {
      return(Matrix::Matrix(x_matrix, sparse = TRUE))
    } else {
      return(x_matrix)
    }
  }
}

#' Validate Matrix Dimensions
#'
#' Validates that two matrices have compatible dimensions for a given operation.
#' This is an internal utility function for dimension checking.
#'
#' @param x First matrix
#' @param y Second matrix  
#' @param operation Type of operation ("multiply", "add", "match_rows", "match_cols")
#' @return Logical indicating whether dimensions are compatible
#' @keywords internal
validate_dimensions <- function(x, y, operation = "multiply") {
  
  #
  dim_x <- dim(x)
  dim_y <- dim(y)
  
  if (is.null(dim_x) || is.null(dim_y)) {
    return(FALSE)
  }
  
  switch(operation,
    "multiply" = dim_x[2] == dim_y[1],
    "add" = all(dim_x == dim_y),
    "match_rows" = dim_x[1] == dim_y[1],
    "match_cols" = dim_x[2] == dim_y[2],
    FALSE
  )
}

#' Summary Statistics for Matrix
#'
#' Computes comprehensive summary statistics for a matrix, handling
#' missing values and providing both numerical and distributional summaries.
#' This is primarily used for debugging and internal diagnostics.
#'
#' @param x Numeric matrix
#' @param na.rm Whether to remove NA values (default: TRUE)
#' @return List containing various summary statistics
#' @keywords internal
matrix_summary <- function(x, na.rm = TRUE) {
  
  #
  dims <- dim(x)
  n_elements <- prod(dims)
  
  # Handle missing values
  if (na.rm) {
    x_clean <- x[!is.na(x)]
    n_missing <- n_elements - length(x_clean)
    prop_missing <- n_missing / n_elements
  } else {
    x_clean <- as.vector(x)
    n_missing <- sum(is.na(x_clean))
    prop_missing <- n_missing / n_elements
  }
  
  # Check if all values are missing
  if (length(x_clean) == 0) {
    return(list(
      dimensions = dims,
      n_elements = n_elements,
      n_missing = n_missing,
      prop_missing = prop_missing,
      summary = "All values are missing"
    ))
  }
  
  list(
    dimensions = dims,
    n_elements = n_elements,
    n_missing = n_missing,
    prop_missing = prop_missing,
    min = min(x_clean),
    max = max(x_clean),
    median = median(x_clean),
    mean = mean(x_clean),
    sd = sd(x_clean),
    range = diff(range(x_clean)),
    quantiles = quantile(x_clean, probs = c(0.05, 0.25, 0.75, 0.95))
  )
}

#' Progress Bar for Long Operations
#'
#' Displays a progress bar for tracking long-running operations.
#' This is primarily used internally by the package functions.
#'
#' @param current Current iteration number
#' @param total Total number of iterations
#' @param width Width of the progress bar in characters (default: 50)
#' @param prefix Text prefix for the progress bar (default: "Progress")
#' @return NULL (prints progress bar to console)
#' @keywords internal
progress_bar <- function(current, total, width = 50, prefix = "Progress") {
  
  #
  percent <- round((current / total) * 100, 1)
  filled <- round((current / total) * width)
  empty <- width - filled
  
  bar <- paste0(
    prefix, ": [",
    paste(rep("=", filled), collapse = ""),
    paste(rep("-", empty), collapse = ""),
    "] ", percent, "% (", current, "/", total, ")"
  )
  
  cat("\r", bar)
  if (current == total) cat("\n")
  utils::flush.console()
}