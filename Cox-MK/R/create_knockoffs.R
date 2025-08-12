#' Create Model-X Knockoffs
#'
#' Generate knockoff variables for genotype data using the Model-X knockoff 
#' method with leveraging scores and clustering.
#'
#' @param X A sparse matrix (n x p) of genotype data where n is the number of 
#'   samples and p is the number of SNPs. Should be of class dgCMatrix.
#' @param pos A numeric vector of SNP positions (in base pairs).
#' @param M Number of knockoff copies to generate. Default is 5, but users can 
#'   specify any positive integer. More knockoff copies may improve power but 
#'   increase computational cost.
#' @param corr_max Maximum correlation threshold for clustering (default: 0.75).
#' @param maxN.neighbor Maximum number of neighbors to consider (default: Inf).
#' @param maxBP.neighbor Maximum base pair distance for neighbors (default: 1e5).
#' @param n.AL Number of samples for adaptive lasso (default: automatic).
#' @param thres.ultrarare Threshold for ultra-rare variants (default: 25).
#' @param R2.thres R-squared threshold (default: 1).
#' @param prob.eps Minimum probability value (default: 1e-12).
#' @param irlba.maxit Maximum iterations for irlba (default: 1500).
#' @param bigmemory Whether to use big.matrix for storage (default: FALSE).
#'
#' @return A list of M matrices, each containing knockoff variables.
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(example_genotypes)
#' data(example_positions)
#' 
#' # Create knockoffs with different numbers of copies
#' knockoffs_3 <- create_knockoffs(example_genotypes, example_positions, M = 3)
#' knockoffs_5 <- create_knockoffs(example_genotypes, example_positions, M = 5)
#' knockoffs_10 <- create_knockoffs(example_genotypes, example_positions, M = 10)
#' }
#'
#' @export
#' @importFrom Matrix Matrix
#' @importFrom irlba irlba
#' @importFrom stats cutree hclust as.dist coef lm.fit
create_knockoffs <- function(
    X, pos, M = 5,
    corr_max = 0.75,
    maxN.neighbor = Inf,
    maxBP.neighbor = 1e5,
    n.AL = floor(10 * nrow(X)^(1/3) * log(nrow(X))),
    thres.ultrarare = 25,
    R2.thres = 1,
    prob.eps = 1e-12,
    irlba.maxit = 1500,
    bigmemory = FALSE) {
  
  # Input validation
  if (!inherits(X, "dgCMatrix")) X <- Matrix::Matrix(X, sparse = TRUE)
  if (length(pos) != ncol(X)) {
    stop("Length of pos must equal number of columns in X")
  }
  if (ncol(X) == 0 || nrow(X) == 0) {
    stop("X must have positive dimensions")
  }
  
  n <- nrow(X)
  
  ## ---- leverage probability calculation ----
  safe_irlba <- function(A, nv) {
    nv <- max(1, min(nv, min(dim(A)) - 1L))
    # For small matrices, use regular SVD instead of irlba
    if (min(dim(A)) <= 10 || nv >= min(dim(A)) * 0.5) {
      tryCatch({
        svd_result <- svd(as.matrix(A), nu = nv, nv = nv)
        # Convert to irlba format
        list(u = svd_result$u, d = svd_result$d[1:nv], v = svd_result$v)
      }, error = function(e) NULL)
    } else {
      tryCatch(irlba::irlba(A, nv = nv, maxit = irlba.maxit),
               error = function(e) NULL,
               warning = function(w) {
                 if (grepl("did not converge", w$message, ignore.case = TRUE))
                   invokeRestart("muffleWarning")
                 else warning(w)
                 NULL
               })
    }
  }
  
  fit <- safe_irlba(X, floor(sqrt(ncol(X) * log(ncol(X)))))
  if (is.null(fit)) {
    prob <- rep(1 / n, n)
  } else {
    h <- rowSums(fit$u^2)
    h[!is.finite(h)] <- 0
    prob <- 0.5 * h / sum(h) + 0.5 / n
  }
  prob <- pmax(prob, prob.eps)
  prob <- prob / sum(prob)
  
  ## ---- subsample rows ----
  n.AL <- min(n.AL, n)
  idx.AL <- sample.int(n, n.AL, FALSE, prob)
  w <- 1 / sqrt(n.AL * prob[idx.AL])
  X.AL <- w * X[idx.AL, , drop = FALSE]
  
  # Ensure X.AL is still a matrix-like object
  if (!is.matrix(X.AL) && !inherits(X.AL, "Matrix")) {
    X.AL <- Matrix::Matrix(X.AL, sparse = TRUE)
  }
  
  ## ---- correlation & clusters ----
  cor.X <- sparse_cor(X)$cor
  clusters <- if (ncol(X) > 1)
    cutree(hclust(as.dist(1 - abs(cor.X)), "single"), h = 1 - corr_max)
  else 1L
  
  # Ensure X.AL has correct dimensions for colSums
  if (ncol(X.AL) == 0) {
    skip <- logical(0)
  } else {
    # Ensure X.AL is 2D for colSums
    if (is.vector(X.AL) || length(dim(X.AL)) != 2) {
      X.AL <- as.matrix(X.AL)
    }
    skip <- colSums(as.matrix(X.AL) != 0) <= thres.ultrarare
  }
  
  ## ---- container ----
  X_k <- lapply(1:M, function(i) {
    if (bigmemory && requireNamespace("bigmemory", quietly = TRUE)) {
      bigmemory::big.matrix(n, ncol(X), init = 0)
    } else {
      matrix(0, n, ncol(X))
    }
  })
  
  ## ---- loop over clusters ----
  for (cl in unique(clusters)) {
    ids <- which(clusters == cl)
    fitv <- res <- matrix(0, n, length(ids))
    
    for (jj in seq_along(ids)) {
      j <- ids[jj]
      y <- X[, j]
      
      idx.pos <- which(pos >= pos[j] - maxBP.neighbor &
                       pos <= pos[j] + maxBP.neighbor)
      tmp <- abs(cor.X[j, ])
      tmp[clusters == cl] <- 0
      tmp[-idx.pos] <- 0
      idx <- order(tmp, decreasing = TRUE)[1:100]
      idx <- setdiff(idx, j)
      idx <- na.omit(idx)
      
      if (length(idx) && !skip[j]) {
        x <- as.matrix(X[, idx, drop = FALSE])
        beta <- coef(lm.fit(cbind(1, x), y))
        fitv[, jj] <- cbind(1, x) %*% beta
      }
      res[, jj] <- y - fitv[, jj]
    }
    
    samp <- replicate(M, sample.int(n))
    for (k in seq_len(M))
      X_k[[k]][, ids] <- round(fitv + res[samp[, k], ], 1)
  }
  
  X_k
}

#' Helper function for sparse correlation calculation
#'
#' @param x A sparse matrix
#' @return A list containing covariance and correlation matrices
#' @keywords internal
sparse_cor <- function(x) {
  # Ensure x is a matrix-like object
  if (is.vector(x)) {
    x <- matrix(x, ncol = 1)
  } else if (!is.matrix(x) && !inherits(x, "Matrix")) {
    # Try to convert to matrix
    x <- as.matrix(x)
  }
  
  # Check dimensions after conversion
  if (is.null(dim(x)) || length(dim(x)) != 2) {
    stop("Input must be convertible to a 2D matrix")
  }
  
  # Handle edge case of single column
  if (ncol(x) == 1) {
    var_val <- if (nrow(x) > 1) var(x[,1]) else 0
    return(list(cov = matrix(var_val, 1, 1), 
                cor = matrix(1, 1, 1)))
  }
  
  n <- nrow(x)
  
  # Use as.matrix for colMeans to ensure compatibility
  if (inherits(x, "Matrix")) {
    cm <- colMeans(as.matrix(x))
  } else {
    cm <- colMeans(x)
  }
  
  cov <- (as.matrix(crossprod(x)) - n * tcrossprod(cm)) / (n - 1)
  sdv <- sqrt(diag(cov))
  
  # Handle zero variance
  sdv[sdv == 0] <- 1
  
  list(cov = cov, cor = cov / tcrossprod(sdv))
}

#' Helper function for sparse cross-covariance calculation
#'
#' @param x A sparse matrix
#' @param y A sparse matrix
#' @return Cross-covariance matrix
#' @keywords internal
sparse_cov_cross <- function(x, y) {
  n <- nrow(x)
  (as.matrix(crossprod(x, y)) -
      n * tcrossprod(colMeans(x), colMeans(y))) / (n - 1)
}
