#' @title Utility Functions
#' @name utilities
#' @description
#' Helper functions for matrix operations, correlation calculations, and
#' other computational utilities used throughout the CoxMK package.
NULL

#' Efficient Sparse Correlation Matrix Computation
#'
#' Computes correlation and covariance matrices for sparse matrices efficiently,
#' particularly useful for large genetic datasets.
#'
#' @param x A numeric matrix (can be sparse)
#' @return List containing:
#'   \item{cov}{The covariance matrix}
#'   \item{cor}{The correlation matrix}
#' @export
#' @examples
#' \dontrun{
#' # Create example sparse matrix
#' X <- Matrix::Matrix(matrix(rnorm(1000), 100, 10), sparse = TRUE)
#' 
#' # Compute correlation
#' result <- sparse_cor(x)
#' cor_matrix <- result$cor
#' cov_matrix <- result$cov
#' }
sparse_cor <- function(x) {
  n <- nrow(x)
  cMeans <- colMeans(x)
  covmat <- (as.matrix(crossprod(x)) - n * tcrossprod(cMeans)) / (n - 1)
  sdvec <- sqrt(diag(covmat))
  cormat <- covmat / tcrossprod(sdvec)
  list(cov = covmat, cor = cormat)
}

#' Sparse Cross-Covariance Matrix Computation
#'
#' Computes cross-covariance matrix between two sparse matrices efficiently.
#'
#' @param x A numeric matrix (can be sparse)
#' @param y A numeric matrix (can be sparse)
#' @return List containing:
#'   \item{cov}{The cross-covariance matrix}
#' @export
#' @examples
#' \dontrun{
#' X <- Matrix::Matrix(matrix(rnorm(500), 50, 10), sparse = TRUE)
#' Y <- Matrix::Matrix(matrix(rnorm(500), 50, 10), sparse = TRUE)
#' 
#' cross_cov <- sparse_cov_cross(X, Y)
#' }
sparse_cov_cross <- function(x, y) {
  n <- nrow(x)
  cMeans.x <- colMeans(x)
  cMeans.y <- colMeans(y)
  covmat <- (as.matrix(crossprod(x, y)) - n * tcrossprod(cMeans.x, cMeans.y)) / (n - 1)
  list(cov = covmat)
}

#' Find N-th Largest Value in Vector
#'
#' Efficiently finds the n-th largest value in a vector using partial sorting.
#'
#' @param x A numeric vector
#' @param n The position (1 = largest, 2 = second largest, etc.)
#' @return The n-th largest value in the vector
#' @export
#' @examples
#' values <- c(1.5, 3.2, 2.8, 4.1, 2.3)
#' max_nth(values, 1)  # 4.1 (largest)
#' max_nth(values, 2)  # 3.2 (second largest)
max_nth <- function(x, n) {
  if (n > length(x)) {
    stop("n cannot be larger than length of x")
  }
  return(sort(x, partial = length(x) - (n - 1))[length(x) - (n - 1)])
}

#' Safe Matrix Conversion
#'
#' Safely converts various matrix-like objects to standard matrix format,
#' handling sparse matrices and potential memory issues.
#'
#' @param x Matrix-like object (matrix, Matrix, data.frame)
#' @param sparse Whether to maintain sparsity if possible (default: FALSE)
#' @return Matrix object
#' @export
safe_as_matrix <- function(x, sparse = FALSE) {
  
  if (inherits(x, "Matrix")) {
    if (sparse) {
      return(x)
    } else {
      return(as.matrix(x))
    }
  } else if (is.data.frame(x)) {
    # Convert data.frame to matrix
    x_mat <- as.matrix(x)
    if (sparse && requireNamespace("Matrix", quietly = TRUE)) {
      return(Matrix::Matrix(x_mat, sparse = TRUE))
    } else {
      return(x_mat)
    }
  } else if (is.matrix(x)) {
    if (sparse && requireNamespace("Matrix", quietly = TRUE)) {
      return(Matrix::Matrix(x, sparse = TRUE))
    } else {
      return(x)
    }
  } else {
    stop("x must be a matrix, Matrix, or data.frame")
  }
}

#' Validate Matrix Dimensions
#'
#' Checks if two matrices have compatible dimensions for specific operations.
#'
#' @param x First matrix
#' @param y Second matrix
#' @param operation Type of operation ("multiply", "add", "match_rows", "match_cols")
#' @return Logical indicating if dimensions are compatible
#' @export
validate_dimensions <- function(x, y, operation = "multiply") {
  
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
#'
#' @param x Numeric matrix
#' @param na.rm Whether to remove NA values (default: TRUE)
#' @return List with summary statistics
#' @export
matrix_summary <- function(x, na.rm = TRUE) {
  
  if (!is.numeric(x)) {
    stop("x must be numeric")
  }
  
  # Basic dimensions
  dims <- dim(x)
  n_elements <- prod(dims)
  
  # Missing values
  n_missing <- sum(is.na(x))
  prop_missing <- n_missing / n_elements
  
  # Numerical summaries
  if (na.rm && n_missing > 0) {
    x_clean <- x[!is.na(x)]
  } else {
    x_clean <- as.vector(x)
  }
  
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
    mean = mean(x_clean),
    median = median(x_clean),
    sd = sd(x_clean),
    range = diff(range(x_clean)),
    quantiles = quantile(x_clean, probs = c(0.05, 0.25, 0.75, 0.95))
  )
}

#' Progress Bar for Long Operations
#'
#' Simple text-based progress bar for tracking long-running operations.
#'
#' @param current Current iteration number
#' @param total Total number of iterations
#' @param width Width of progress bar in characters (default: 50)
#' @param prefix Optional prefix text
#' @return None (prints progress bar)
#' @export
progress_bar <- function(current, total, width = 50, prefix = "Progress") {
  
  percent <- current / total
  filled <- round(width * percent)
  bar <- paste0(rep("=", filled), rep("-", width - filled), collapse = "")
  
  cat(sprintf("\r%s: [%s] %d%% (%d/%d)", 
              prefix, bar, round(percent * 100), current, total))
  
  if (current == total) {
    cat("\n")
  }
  
  flush.console()
}
