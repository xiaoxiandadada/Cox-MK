#' @title Knockoff Filter and Variable Selection
#' @name knockoff-filter
#' @description
#' Functions for computing W statistics and applying knockoff filters
#' for variable selection with false discovery rate control.
NULL

#' Calculate W Statistics for Knockoff Analysis
#'
#' Computes W statistics by comparing test statistics from original variables
#' with those from their knockoff counterparts. These statistics are used
#' for variable selection with FDR control.
#'
#' @param t_orig Vector of test statistics for original variables
#' @param t_knock Vector or list of test statistics for knockoff variables.
#'   If a list, should contain M vectors of the same length as t_orig.
#' @param method Method for computing W statistics:
#'   \itemize{
#'     \item "difference": W_j = T_j - max(T_{j,k}) (default)
#'     \item "mk_median": Uses Model-X knockoff median-based statistics
#'     \item "ratio": W_j = T_j / max(T_{j,k})
#'   }
#' @return Vector of W statistics for variable selection
#' @export
#' @examples
#' \dontrun{
#' # Example with difference method
#' t_orig <- c(5.2, 3.1, 8.7, 2.4, 6.9)
#' t_knock <- list(
#'   c(2.1, 4.2, 3.3, 1.8, 2.9),
#'   c(1.9, 3.8, 4.1, 2.2, 3.1)
#' )
#' 
#' w_diff <- calculate_w_statistics(t_orig, t_knock, method = "difference")
#' w_mk <- calculate_w_statistics(t_orig, t_knock, method = "mk_median")
#' }
calculate_w_statistics <- function(t_orig, t_knock, method = "difference") {
  
  # Input validation
  if (!is.numeric(t_orig)) {
    stop("t_orig must be numeric")
  }
  
  if (is.list(t_knock)) {
    # Multiple knockoff copies
    t_knock_matrix <- do.call(cbind, t_knock)
  } else if (is.vector(t_knock)) {
    # Single knockoff copy
    t_knock_matrix <- matrix(t_knock, ncol = 1)
  } else {
    t_knock_matrix <- as.matrix(t_knock)
  }
  
  if (length(t_orig) != nrow(t_knock_matrix)) {
    stop("Length of t_orig must match number of rows in t_knock")
  }
  
  # Calculate W statistics based on method
  if (method == "difference") {
    # Standard knockoff: W_j = T_j - max(T_j^knockoff)
    t_knock_max <- apply(t_knock_matrix, 1, max, na.rm = TRUE)
    W <- t_orig - t_knock_max
    
  } else if (method == "mk_median") {
    # Model-X knockoff with median method
    mk_result <- MK.statistic(t_orig, t_knock_matrix, method = "median")
    W <- mk_result[, "tau"]
    
  } else if (method == "ratio") {
    # Ratio-based statistics
    t_knock_max <- apply(t_knock_matrix, 1, max, na.rm = TRUE)
    W <- ifelse(t_knock_max > 0, t_orig / t_knock_max, t_orig)
    
  } else {
    stop("Method must be one of: 'difference', 'mk_median', 'ratio'")
  }
  
  return(W)
}

#' Apply Knockoff Filter for Variable Selection
#'
#' Applies the knockoff filter to select variables while controlling the
#' false discovery rate (FDR) at a specified level.
#'
#' @param W Vector of W statistics from \code{\link{calculate_w_statistics}}
#' @param fdr Target false discovery rate (default: 0.1)
#' @param offset Offset parameter for knockoff filter (default: 1)
#' @return Vector of indices of selected variables
#' @export
#' @examples
#' \dontrun{
#' # Generate some example W statistics
#' W <- c(2.1, -0.5, 3.8, -1.2, 4.5, 0.3, -2.1, 1.9)
#' 
#' # Apply knockoff filter
#' selected <- knockoff_filter(W, fdr = 0.1)
#' print(selected)  # Indices of selected variables
#' }
knockoff_filter <- function(W, fdr = 0.1, offset = 1) {
  
  if (!is.numeric(W)) {
    stop("W must be numeric")
  }
  
  if (fdr <= 0 || fdr >= 1) {
    stop("fdr must be between 0 and 1")
  }
  
  # Sort W statistics in decreasing order
  W_sorted_idx <- order(W, decreasing = TRUE)
  W_sorted <- W[W_sorted_idx]
  
  # Calculate knockoff threshold
  n_positive <- cumsum(W_sorted > 0)
  n_negative <- cumsum(W_sorted <= 0)
  
  # FDR control: find largest k such that (offset + n_negative[k]) / max(1, n_positive[k]) <= fdr
  ratio <- (offset + n_negative) / pmax(1, n_positive)
  k_max <- max(which(ratio <= fdr), 0)
  
  if (k_max > 0) {
    threshold <- W_sorted[k_max]
    selected <- which(W >= threshold & W > 0)
  } else {
    threshold <- Inf
    selected <- integer(0)
  }
  
  # Add threshold as attribute for reference
  attr(selected, "threshold") <- threshold
  attr(selected, "fdr") <- fdr
  
  return(selected)
}

#' Model-X Knockoff Statistics Computation
#'
#' Computes kappa and tau statistics for Model-X knockoff methodology.
#' This is an internal function called by calculate_w_statistics when
#' method = "mk_median".
#'
#' @param T_0 Test statistics for original variables
#' @param T_k Matrix of test statistics for knockoff variables
#' @param method Method for tau calculation ("median" or "max")
#' @return Matrix with kappa and tau columns
#' @keywords internal
MK.statistic <- function(T_0, T_k, method = "median") {
  T_0 <- as.matrix(T_0)
  T_k <- as.matrix(T_k)
  T.temp <- cbind(T_0, T_k)
  T.temp[is.na(T.temp)] <- 0

  # Find which knockoff has maximum value (kappa statistic)
  which.max.alt <- function(x) {
    temp.index <- which(x == max(x))
    if (length(temp.index) != 1) {
      return(temp.index[2])  # If tie, use second occurrence
    } else {
      return(temp.index[1])
    }
  }
  kappa <- apply(T.temp, 1, which.max.alt) - 1

  # Calculate tau statistic
  if (method == "max") {
    tau <- apply(T.temp, 1, max) - apply(T.temp, 1, max_nth, n = 2)
  } else if (method == "median") {
    tau <- apply(T.temp, 1, max) - apply(T.temp, 1, median)
  }
  
  return(cbind(kappa = kappa, tau = tau))
}

#' Calculate Model-X Knockoff Threshold
#'
#' Computes the threshold for variable selection using Model-X knockoff methodology.
#'
#' @param T_0 Test statistics for original variables
#' @param T_k Test statistics for knockoff variables
#' @param fdr Target false discovery rate
#' @param method Method for calculating statistics
#' @param Rej.Bound Maximum number of rejections
#' @return Threshold value
#' @export
MK.threshold <- function(T_0, T_k, fdr = 0.1, method = "median", Rej.Bound = 10000) {
  stat <- MK.statistic(T_0, T_k, method = method)
  kappa <- stat[, "kappa"]
  tau <- stat[, "tau"]
  t <- MK.threshold.byStat(kappa, tau, M = ncol(as.matrix(T_k)), fdr = fdr, Rej.Bound = Rej.Bound)
  return(t)
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
  b <- order(tau, decreasing = TRUE)
  c_0 <- kappa[b] == 0
  ratio <- c()
  temp_0 <- 0
  for (i in seq_along(b)) {
    if (c_0[i] == TRUE) {
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

#' Calculate Model-X Knockoff Q-values by Statistics
#'
#' @param kappa Vector of kappa statistics
#' @param tau Vector of tau statistics
#' @param M Number of knockoff copies
#' @param Rej.Bound Maximum number of rejections
#' @return Vector of q-values
#' @export
MK.q.byStat <- function(kappa, tau, M, Rej.Bound = 10000) {
  b <- order(tau, decreasing = TRUE)
  c_0 <- kappa[b] == 0
  ratio <- c()
  temp_0 <- 0
  
  for (i in seq_along(b)) {
    if (c_0[i] == TRUE) {
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
