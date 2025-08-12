###############################################################################
#  KnockoffScreen.r
###############################################################################
suppressPackageStartupMessages({
  library(Matrix)
  library(bigmemory)
  library(irlba)
})

## ---------- helpers ---------------------------------------------------------
sparse.cor <- function(x) {
  n <- nrow(x)
  cm <- colMeans(x)
  cov <- (as.matrix(crossprod(x)) - n * tcrossprod(cm)) / (n - 1)
  sdv <- sqrt(diag(cov))
  list(cov = cov, cor = cov / tcrossprod(sdv))
}

sparse.cov.cross <- function(x, y) {
  n <- nrow(x)
  (as.matrix(crossprod(x, y)) -
      n * tcrossprod(colMeans(x), colMeans(y))) / (n - 1)
}

## ---------- shrinkage create.MK ---------------------------------------
create.MK <- function(
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

  if (!inherits(X, "dgCMatrix")) X <- Matrix(X, sparse = TRUE)
  n <- nrow(X)

  ## ---- leverage prob -------------------------------------------------------
  safe_irlba <- function(A, nv) {
    nv <- max(1, min(nv, min(dim(A)) - 1L))
    tryCatch(irlba(A, nv = nv, maxit = irlba.maxit),
             error = function(e) NULL,
             warning = function(w) {
               if (grepl("did not converge", w$message, ignore.case = TRUE))
                 invokeRestart("muffleWarning")
               else warning(w)
               NULL
             })
  }
  fit <- safe_irlba(X, floor(sqrt(ncol(X) * log(ncol(X)))))
  if (is.null(fit)) {
    prob <- rep(1 / n, n)
  } else {
    h <- rowSums(fit$u^2); h[!is.finite(h)] <- 0
    prob <- 0.5 * h / sum(h) + 0.5 / n
  }
  prob <- pmax(prob, prob.eps); prob <- prob / sum(prob)

  ## ---- subsample rows ------------------------------------------------------
  n.AL <- min(n.AL, n)
  idx.AL <- sample.int(n, n.AL, FALSE, prob)
  w <- 1 / sqrt(n.AL * prob[idx.AL])
  X.AL <- w * X[idx.AL, , drop = FALSE]

  ## ---- correlation & clusters ---------------------------------------------
  cor.X <- sparse.cor(X)$cor
  clusters <- if (ncol(X) > 1)
      cutree(hclust(as.dist(1 - abs(cor.X)), "single"), h = 1 - corr_max)
    else 1L
  skip <- colSums(X.AL != 0) <= thres.ultrarare

  ## ---- container -----------------------------------------------------------
  X_k <- lapply(1:M, function(i)
    if (bigmemory) big.matrix(n, ncol(X), init = 0)
    else           matrix(0, n, ncol(X)))

  ## ---- loop over clusters --------------------------------------------------
  for (cl in unique(clusters)) {
    ids  <- which(clusters == cl)
    fitv <- res  <- matrix(0, n, length(ids))

    for (jj in seq_along(ids)) {
      j <- ids[jj]; y <- X[, j]

      idx.pos <- which(pos >= pos[j] - maxBP.neighbor &
                       pos <= pos[j] + maxBP.neighbor)
      tmp <- abs(cor.X[j, ]); tmp[clusters == cl] <- 0; tmp[-idx.pos] <- 0
      idx <- order(tmp, decreasing = TRUE)[1:100]
      idx <- setdiff(idx, j); idx <- na.omit(idx)

      if (length(idx) && !skip[j]) {
        x <- as.matrix(X[, idx, drop = FALSE])
        beta <- coef(lm.fit(cbind(1, x), y))
        fitv[, jj] <- cbind(1, x) %*% beta
      }
      res[, jj]  <- y - fitv[, jj]
    }

    samp <- replicate(M, sample.int(n))
    for (k in seq_len(M))
      X_k[[k]][, ids] <- round(fitv + res[samp[, k], ], 1)
  }

  X_k
}