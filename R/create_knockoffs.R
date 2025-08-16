#' @title Knockoff Variable Generation
#' @name knockoff-generation
#' @description
#' Functions for generating Model-X knockoff variables using leveraging scores,
#' clustering, and adaptive lasso methods specifically designed for genetic data.
#' @importFrom Matrix Matrix crossprod
#' @importFrom irlba irlba
#' @importFrom stats cutree hclust as.dist coef lm.fit cov2cor cor
#' @importFrom gdsfmt createfn.gds add.gdsn closefn.gds
NULL

#' Create Model-X Knockoffs for Genetic Data
#'
#' Generate knockoff variables for genotype data using the Model-X knockoff 
#' method with leveraging scores and clustering specifically optimized for
#' genetic variant data.
#'
#' @param X A sparse matrix (n x p) of genotype data where n is the number of 
#'   samples and p is the number of SNPs. Typically coded as 0, 1, 2 for
#'   genotype dosages.
#' @param pos A numeric vector of SNP positions (in base pairs) for linkage
#'   disequilibrium-aware knockoff generation.
#' @param chr_info Optional data frame with chromosome information from BIM file.
#'   If provided, will extract chromosome number automatically.
#' @param sample_ids A character vector of sample IDs (default: NULL, will generate)
#' @param M Number of knockoff copies to generate (default: 5). More copies
#'   can improve statistical power but increase computational cost.
#' @param save_gds Whether to save knockoffs to GDS format (default: TRUE)
#' @param output_dir Directory to save GDS files (default: extdata folder)
#' @param start Start position for file naming (default: min(pos))
#' @param end End position for file naming (default: max(pos))
#' @param corr_max Maximum correlation threshold for clustering variants
#'   (default: 0.75). Higher values create fewer, larger clusters.
#' @param maxN.neighbor Maximum number of neighboring variants to consider 
#'   for each variant (default: Inf).
#' @param maxBP.neighbor Maximum base pair distance to consider variants as
#'   neighbors (default: 100,000 bp).
#' @param n.AL Number of samples to use for adaptive lasso fitting 
#'   (default: automatically determined based on sample size).
#' @param thres.ultrarare Minimum minor allele count threshold for variant
#'   inclusion (default: 25).
#' @param R2.thres R-squared threshold for model fitting (default: 1).
#' @param prob.eps Minimum probability value to prevent numerical issues
#'   (default: 1e-12).
#' @param irlba.maxit Maximum iterations for truncated SVD (default: 1500).
#'
#' @return If save_gds is TRUE, returns the path to the saved GDS file.
#'   Otherwise, returns a list of M matrices, each of the same dimensions as X, 
#'   containing knockoff variables.
#'
#' @export
create_knockoffs <- function(
    X, pos, chr_info = NULL, sample_ids = NULL, M = 5,
    save_gds = TRUE,
    output_dir = NULL,
    start = NULL,
    end = NULL,
    corr_max = 0.75,
    maxN.neighbor = Inf,
    maxBP.neighbor = 1e5,
    n.AL = floor(10 * nrow(X)^(1/3) * log(nrow(X))),
    thres.ultrarare = 25,
    R2.thres = 1,
    prob.eps = 1e-12,
    irlba.maxit = 1500) {

  ## ---------- helpers ---------------------------------------------------------
  sparse.cor <- function(x) {
    if (is.null(dim(x)) || length(dim(x)) < 2) {
      stop("Input must be a matrix with at least 2 dimensions")
    }
    n <- nrow(x)
    if (n < 2) {
      stop("Matrix must have at least 2 rows")
    }
    # Convert to regular matrix for colMeans
    if (inherits(x, "Matrix")) {
      cm <- Matrix::colMeans(x)
    } else {
      cm <- colMeans(x)
    }
    cov <- (as.matrix(crossprod(x)) - n * tcrossprod(cm)) / (n - 1)
    sdv <- sqrt(diag(cov))
    list(cov = cov, cor = cov / tcrossprod(sdv))
  }

  ## ---------- Input validation ------------------------------------------------
  if (!inherits(X, "dgCMatrix")) X <- Matrix(X, sparse = TRUE)
  n <- nrow(X)
  p <- ncol(X)

  if (length(pos) != p) {
    stop("Length of pos must equal number of columns in X")
  }

  # Extract chromosome number from chr_info if provided
  chr <- 1  # default
  if (!is.null(chr_info) && "CHR" %in% names(chr_info)) {
    unique_chrs <- unique(chr_info$CHR)
    if (length(unique_chrs) == 1) {
      chr <- unique_chrs[1]
    } else {
      chr <- unique_chrs[1]  # Use first chromosome if multiple
      warning("Multiple chromosomes found, using chr ", chr)
    }
  }

  # Set default values
  if (is.null(start)) start <- min(pos)
  if (is.null(end)) end <- max(pos)
  if (is.null(sample_ids)) sample_ids <- paste0("SAMPLE_", seq_len(n))
  if (is.null(output_dir)) {
    # Always save to inst/extdata in the package source directory
    output_dir <- file.path(getwd(), "inst", "extdata")
    if (!dir.exists(output_dir)) {
      # If we're in installed package, use the installed extdata
      installed_dir <- system.file("extdata", package = "CoxMK")
      if (installed_dir != "" && dir.exists(dirname(installed_dir))) {
        output_dir <- installed_dir
      }
    }
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

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
  # Remove X.AL creation from here since it's moved below

  ## ---- correlation & clusters ---------------------------------------------
  if (ncol(X) > 1 && nrow(X) > 1) {
    cor.X <- sparse.cor(X)$cor
    clusters <- cutree(hclust(as.dist(1 - abs(cor.X)), "single"), h = 1 - corr_max)
  } else {
    cor.X <- matrix(1, ncol(X), ncol(X))
    clusters <- rep(1L, ncol(X))
  }
  
  # Create X.AL for ultra-rare threshold calculation
  if (n.AL < n) {
    X.AL <- X[idx.AL, , drop = FALSE]
  } else {
    X.AL <- X
  }
  
  # Calculate skip using Matrix-aware functions
  if (inherits(X.AL, "Matrix")) {
    skip <- Matrix::colSums(X.AL != 0) <= thres.ultrarare
  } else {
    skip <- colSums(X.AL != 0) <= thres.ultrarare
  }

  ## ---- container -----------------------------------------------------------
  X_k <- lapply(1:M, function(i) matrix(0, n, ncol(X)))

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

  cat("Knockoff generation complete!\n")

  # Helper function to convert to dense integer matrix
  to_dense <- function(x) { 
    m <- as.matrix(x)
    mode(m) <- "integer"
    m 
  }

  # Save to GDS format if requested
  if (save_gds) {
    cat("Saving knockoffs to GDS format...\n")
    
    # Create GDS filename
    gds_filename <- sprintf("chr%d_%d_%d_knockoff.gds", chr, start, end)
    gds_path <- file.path(output_dir, gds_filename)
    
    # Convert matrices to dense integer format
    geno_orig <- to_dense(X)
    geno_k <- lapply(X_k, to_dense)
    
    # Create GDS file
    if (requireNamespace("gdsfmt", quietly = TRUE)) {
      g <- gdsfmt::createfn.gds(gds_path)
      
      # Add sample IDs
      gdsfmt::add.gdsn(g, "sample.id", sample_ids, "string")
      
      # Add SNP positions
      gdsfmt::add.gdsn(g, "snp.pos", pos, "int32")
      
      # Add original genotypes
      gdsfmt::add.gdsn(g, "original", geno_orig,
                      storage = "bit2", valdim = dim(geno_orig), 
                      compress = "LZMA_RA")
      
      # Add knockoff genotypes
      for (i in seq_len(M)) {
        gdsfmt::add.gdsn(g, paste0("knockoff", i), geno_k[[i]],
                        storage = "bit2", valdim = dim(geno_k[[i]]), 
                        compress = "LZMA_RA")
      }
      
      gdsfmt::closefn.gds(g)
      cat("GDS written:", gds_path, "\n")
      
      return(list(knockoffs = X_k, gds_file = gds_path))
    } else {
      warning("gdsfmt package not available, returning knockoff matrices instead")
      return(list(knockoffs = X_k, gds_file = NULL))
    }
  }
  
  return(list(knockoffs = X_k, gds_file = NULL))
}
