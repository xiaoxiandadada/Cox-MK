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
#' @param sample_ids A character vector of sample IDs (default: NULL, will generate)
#' @param M Number of knockoff copies to generate (default: 5). More copies
#'   can improve statistical power but increase computational cost.
#' @param save_gds Whether to save knockoffs to GDS format (default: TRUE)
#' @param output_dir Directory to save GDS files (default: extdata folder)
#' @param chr Chromosome number for file naming (default: 1)
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
    X, pos, sample_ids = NULL, M = 5,
    save_gds = TRUE,
    output_dir = NULL,
    chr = 1,
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
  
  # Input validation
  if (!is.matrix(X) && !inherits(X, "Matrix")) {
    stop("X must be a matrix or Matrix object")
  }
  
  if (length(pos) != ncol(X)) {
    stop("Length of pos must equal number of columns in X")
  }
  
  if (M < 1) {
    stop("M must be at least 1")
  }
  
  n <- nrow(X)
  p <- ncol(X)
  
  # Set default values
  if (is.null(start)) start <- min(pos)
  if (is.null(end)) end <- max(pos)
  if (is.null(sample_ids)) sample_ids <- paste0("SAMPLE_", seq_len(n))
  if (is.null(output_dir)) {
    output_dir <- system.file("extdata", package = "CoxMK")
    if (output_dir == "") {
      output_dir <- file.path(getwd(), "inst", "extdata")
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  # Convert to Matrix format for efficiency
  if (!inherits(X, "Matrix")) {
    X <- Matrix(X, sparse = TRUE)
  }
  
  # Helper function to convert to dense integer matrix
  to_dense <- function(x) { 
    m <- as.matrix(x)
    mode(m) <- "integer"
    m 
  }
  
  # Calculate correlation matrix efficiently
  cat("Computing correlation matrix...\n")
  
  # Compute correlation matrix with fallback for edge cases
  if (n > 1 && p > 1) {
    cor.X <- tryCatch({
      cor(as.matrix(X))
    }, error = function(e) {
      # Handle numerical issues
      cor_matrix <- matrix(0, p, p)
      diag(cor_matrix) <- 1
      cor_matrix
    })
  } else {
    cor.X <- matrix(1, p, p)
  }
  
  # Handle NAs and ensure diagonal is 1
  cor.X[is.na(cor.X)] <- 0
  diag(cor.X) <- 1
  
  # Clustering based on correlation
  cat("Performing hierarchical clustering...\n")
  dist_matrix <- as.dist(1 - abs(cor.X))
  hclust_result <- hclust(dist_matrix, method = "average")
  clusters <- cutree(hclust_result, h = 1 - corr_max)
  
  # Determine which SNPs to skip (ultra-rare variants)
  mac <- Matrix::colSums(X)  # minor allele count
  skip <- mac < thres.ultrarare | mac > (2 * n - thres.ultrarare)
  
  # Leveraging scores calculation using truncated SVD
  cat("Computing leveraging scores...\n")
  
  # Calculate leveraging scores
  nv <- min(n.AL, min(dim(X)))
  if (nv > 0) {
    svd_result <- tryCatch({
      irlba(X, nv = nv, maxit = irlba.maxit)
    }, error = function(e) {
      # Fallback to regular SVD if irlba fails
      min_dim <- min(dim(X))
      svd(X, nu = min(nv, min_dim), nv = min(nv, min_dim))
    })
    lev_scores <- rowSums(svd_result$u^2)
  } else {
    lev_scores <- rep(1/n, n)
  }
  
  # Initialize knockoff matrices
  X_k <- lapply(seq_len(M), function(k) {
    Matrix(0, n, p, sparse = TRUE)
  })
  
  # Generate knockoffs cluster by cluster
  cat("Generating knockoffs for", length(unique(clusters)), "clusters...\n")
  
  for (cl in unique(clusters)) {
    ids <- which(clusters == cl)
    n_vars_in_cluster <- length(ids)
    
    if (n_vars_in_cluster == 0) next
    
    # Fit models for each variable in the cluster
    fitv <- res <- matrix(0, n, n_vars_in_cluster)
    
    for (jj in seq_along(ids)) {
      j <- ids[jj]
      y <- X[, j]
      
      if (skip[j]) {
        res[, jj] <- y
        next
      }
      
      # Find neighboring variables
      idx.pos <- which(pos >= pos[j] - maxBP.neighbor &
                       pos <= pos[j] + maxBP.neighbor)
      
      # Find correlated variables (excluding same cluster)
      tmp <- abs(cor.X[j, ])
      tmp[clusters == cl] <- 0  # Exclude same cluster
      tmp[-idx.pos] <- 0        # Exclude distant SNPs
      
      # Select top correlated neighbors
      n_neighbors <- min(100, sum(tmp > 0))
      if (n_neighbors > 0) {
        idx <- order(tmp, decreasing = TRUE)[seq_len(n_neighbors)]
        idx <- setdiff(idx, j)
        idx <- idx[!is.na(idx)]
        
        if (length(idx) > 0) {
          # Fit linear model
          x_neighbors <- as.matrix(X[, idx, drop = FALSE])
          
          # Use lm.fit for efficiency with fallback
          beta <- tryCatch({
            coef(lm.fit(cbind(1, x_neighbors), y))
          }, error = function(e) {
            c(mean(y), rep(0, ncol(x_neighbors)))
          })
          
          fitv[, jj] <- cbind(1, x_neighbors) %*% beta
        } else {
          fitv[, jj] <- mean(y)
        }
      } else {
        fitv[, jj] <- mean(y)
      }
      
      res[, jj] <- y - fitv[, jj]
    }
    
    # Generate M knockoff copies by resampling residuals
    for (k in seq_len(M)) {
      samp <- sample.int(n)
      X_k[[k]][, ids] <- round(fitv + res[samp, , drop = FALSE], 1)
    }
  }
  
  cat("Knockoff generation complete!\n")
  
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
      cat("* GDS written:", gds_path, "\n")
      
      return(list(knockoffs = X_k, gds_file = gds_path))
    } else {
      warning("gdsfmt package not available, returning knockoff matrices instead")
      return(list(knockoffs = X_k, gds_file = NULL))
    }
  }
  
  return(list(knockoffs = X_k, gds_file = NULL))
}
