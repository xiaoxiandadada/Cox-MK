#' Model-X Knockoff Statistics and Fil#' @return The n-th largest value in the vector
#' @export
max.nth <- function(x, n) {
  return(sort(x, partial = length(x) - (n - 1))[length(x) - (n - 1)])
}

#' Sparse Correlation Matrix Computation
#'
#' Computes correlation and covariance matrices for sparse matrices efficiently.
#'
#' @param x A numeric matrix (can be sparse)
#' @return A list containing:
#'   \item{cov}{The covariance matrix}
#'   \item{cor}{The correlation matrix}
#' @export
sparse.cor <- function(x) {
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
#' @return A list containing:
#'   \item{cov}{The cross-covariance matrix}
#' @export
sparse.cov.cross <- function(x, y) {
  n <- nrow(x)
  cMeans.x <- colMeans(x)
  cMeans.y <- colMeans(y)
  covmat <- (as.matrix(crossprod(x, y)) - n * tcrossprod(cMeans.x, cMeans.y)) / (n - 1)
  list(cov = covmat)
}

#' @title Knockoff Statistics Functions
#' @name knockoff-statistics
#' @description Core statistical functions for Model-X knockoffs methodology
#' @importFrom stats median pchisq var
NULL

#' Calculate Model-X Knockoff Statistics
#'
#' @param T_0 Test statistics for original variables
#' @param T_k Test statistics for knockoff variables (matrix)
#' @param method Method for calculating statistics ("median" or "max")
#' @return Matrix with kappa and tau statistics
#' @keywords internal
MK.statistic <- function(T_0, T_k, method = 'median') {
  T_0 <- as.matrix(T_0)
  T_k <- as.matrix(T_k)
  T.temp <- cbind(T_0, T_k)
  T.temp[is.na(T.temp)] <- 0

  which.max.alt <- function(x) {
    temp.index <- which(x == max(x))
    if (length(temp.index) != 1) {
      return(temp.index[2])
    } else {
      return(temp.index[1])
    }
  }
  kappa <- apply(T.temp, 1, which.max.alt) - 1

  if (method == 'max') {
    tau <- apply(T.temp, 1, max) - apply(T.temp, 1, max_nth, n = 2)
  }
  if (method == 'median') {
    tau <- apply(T.temp, 1, max) - apply(T.temp, 1, median)
  }
  return(cbind(kappa, tau))
}

#' Find N-th Largest Value
#'
#' Returns the n-th largest value in a vector
#'
#' @param x A numeric vector
#' @param n The position (1 = largest, 2 = second largest, etc.)
#' @return The n-th largest value in the vector
#' @export
max_nth <- function(x, n) {
  return(sort(x, partial = length(x) - (n - 1))[length(x) - (n - 1)])
}

#' Calculate Model-X Knockoff Threshold by Statistics
#'
#' @param kappa Vector of kappa statistics
#' @param tau Vector of tau statistics  
#' @param M Number of knockoff copies
#' @param fdr Target false discovery rate
#' @param Rej.Bound Maximum number of rejections
#' @return Threshold value
#' @keywords internal
MK.threshold.byStat <- function(kappa, tau, M, fdr = 0.1, Rej.Bound = 10000) {
  b <- order(tau, decreasing = T)
  c_0 <- kappa[b] == 0
  ratio <- c()
  temp_0 <- 0
  for (i in seq_along(b)) {
    if (c_0[i] == T) {
      temp_0 <- temp_0 + 1
    }
    temp_1 <- i - temp_0
    if (temp_1 == 0) {
      ratio <- c(ratio, 1)
    } else {
      ratio <- c(ratio, (temp_0 + 1) / temp_1)
    }
    if (i >= Rej.Bound) break
  }
  ok <- which(ratio <= fdr)
  if (length(ok) > 0) {
    return(tau[b][ok[length(ok)]])
  } else {
    return(Inf)
  }
}

#' Calculate Model-X Knockoff Threshold
#'
#' @param T_0 Test statistics for original variables
#' @param T_k Test statistics for knockoff variables
#' @param fdr Target false discovery rate
#' @param method Method for calculating statistics
#' @param Rej.Bound Maximum number of rejections
#' @return Threshold value
#' @export
MK.threshold <- function(T_0, T_k, fdr = 0.1, method = 'median', Rej.Bound = 10000) {
  stat <- MK.statistic(T_0, T_k, method = method)
  kappa <- stat[, 1]
  tau <- stat[, 2]
  t <- MK.threshold.byStat(kappa, tau, M = ncol(T_k), fdr = fdr, Rej.Bound = Rej.Bound)
  return(t)
}

#' Calculate Model-X Knockoff Q-values by Statistics
#'
#' @param kappa Vector of kappa statistics
#' @param tau Vector of tau statistics
#' @param M Number of knockoff copies
#' @param Rej.Bound Maximum number of rejections
#' @return Vector of q-values
#' @export
MK.q.byStat <- function(kappa, tau, M, Rej.Bound = 10000) {
  b <- order(tau, decreasing = T)
  c_0 <- kappa[b] == 0
  ratio <- c()
  temp_0 <- 0
  
  for (i in seq_along(b)) {
    if (c_0[i] == T) {
      temp_0 <- temp_0 + 1
    }
    temp_1 <- i - temp_0
    if (temp_1 == 0) {
      ratio <- c(ratio, 1)
    } else {
      ratio <- c(ratio, (temp_0 + 1) / temp_1)
    }
    if (i >= Rej.Bound) break
  }
  
  # Calculate cumulative minimum from right to left
  if (length(ratio) > 0) {
    for (i in (length(ratio) - 1):1) {
      ratio[i] <- min(ratio[i], ratio[i + 1])
    }
  }
  
  q <- rep(1, length(tau))
  if (length(ratio) > 0) {
    q[b[seq_along(ratio)]] <- ratio
  }
  
  return(q)
}

#' Get P-values using Score Test
#'
#' @param X Genotype matrix
#' @param result.prelim Preliminary results from KS.prelim
#' @return Vector of p-values
#' @keywords internal
Get.p <- function(X, result.prelim) {
  mu <- result.prelim$nullglm$fitted.values
  Y.res <- result.prelim$Y - mu
  outcome <- result.prelim$out_type
  
  if (outcome == 'D') {
    # For binary outcomes, use SPA test if available
    if (requireNamespace("SPAtest", quietly = TRUE)) {
      p <- SPAtest::ScoreTest_SPA(t(X), result.prelim$Y, result.prelim$X, 
                                  method = c("fastSPA"), minmac = -Inf)$p.value
    } else {
      # Fallback to score test
      v <- mu * (1 - mu)
      A <- (t(X) %*% Y.res)^2
      B <- colSums(v * X^2)
      C <- t(X) %*% (v * result.prelim$X0) %*% result.prelim$inv.X0
      D <- t(t(result.prelim$X0) %*% as.matrix(v * X))
      p <- pchisq(as.numeric(A / (B - rowSums(C * D))), df = 1, lower.tail = F)
    }
  } else {
    # For continuous outcomes
    v <- rep(as.numeric(var(Y.res)), nrow(X))
    A <- (t(X) %*% Y.res)^2
    B <- colSums(v * X^2)
    C <- t(X) %*% (v * result.prelim$X0) %*% result.prelim$inv.X0
    D <- t(t(result.prelim$X0) %*% as.matrix(v * X))
    p <- pchisq(as.numeric(A / (B - rowSums(C * D))), df = 1, lower.tail = F)
  }
  
  return(as.matrix(p))
}

#' Cauchy Combination Function
#'
#' @param p Vector of p-values
#' @return Combined p-value using Cauchy combination
#' @keywords internal
Get.cauchy <- function(p) {
  # Remove NA values
  p <- p[!is.na(p)]
  if (length(p) == 0) return(NA)
  
  # Cauchy combination
  C <- mean(tan((0.5 - p) * pi))
  p.cauchy <- 0.5 - atan(C) / pi
  
  return(p.cauchy)
}

#' Cauchy Combination for Sliding Windows
#'
#' @param p Vector of p-values
#' @param window.matrix Window indicator matrix
#' @return Vector of combined p-values
#' @keywords internal
Get.cauchy.scan <- function(p, window.matrix) {
  if (is.vector(window.matrix)) {
    window.matrix <- matrix(window.matrix, ncol = 1)
  }
  
  result <- apply(window.matrix, 2, function(w) {
    if (sum(w) == 0) return(NA)
    selected.p <- p[w == 1]
    Get.cauchy(selected.p)
  })
  
  return(result)
}

#' Imputation Function
#'
#' @param Z Matrix with missing values
#' @param impute.method Imputation method ("fixed" or "mean")
#' @return Imputed matrix
#' @keywords internal
Impute <- function(Z, impute.method) {
  if (impute.method == 'fixed') {
    # Fixed value imputation (use mean)
    for (j in seq_len(ncol(Z))) {
      missing.index <- is.na(Z[, j])
      if (sum(missing.index) > 0) {
        Z[missing.index, j] <- mean(Z[!missing.index, j], na.rm = TRUE)
      }
    }
  }
  return(Z)
}
