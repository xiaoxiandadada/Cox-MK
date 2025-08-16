library(testthat)
library(CoxMK)

test_that("Package loads correctly", {
  expect_true("CoxMK" %in% loadedNamespaces())
})

test_that("Example data loads correctly", {
  extdata_path <- system.file("extdata", package = "CoxMK")
  expect_true(dir.exists(extdata_path))
  
  # Check PLINK files exist
  expect_true(file.exists(file.path(extdata_path, "sample.bed")))
  expect_true(file.exists(file.path(extdata_path, "sample.bim")))  
  expect_true(file.exists(file.path(extdata_path, "sample.fam")))
  
  # Check phenotype and covariate files exist
  expect_true(file.exists(file.path(extdata_path, "tte_phenotype.txt")))
  expect_true(file.exists(file.path(extdata_path, "covariates.txt")))
})

test_that("Core functions work", {
  # Load example data from extdata
  extdata_path <- system.file("extdata", package = "CoxMK")
  
  # Load PLINK data
  plink_data <- load_plink_data(file.path(extdata_path, "sample"))
  
  # Test knockoff creation with small subset
  n_samples <- min(10, nrow(plink_data$genotypes))
  n_snps <- min(5, ncol(plink_data$genotypes))
  
  subset_X <- plink_data$genotypes[seq_len(n_samples), seq_len(n_snps)]
  subset_pos <- plink_data$positions[seq_len(n_snps)]
  
  knockoffs <- create_knockoffs(
    X = subset_X,
    pos = subset_pos,
    M = 2,
    save_gds = FALSE  # Don't save GDS for testing
  )
  
  expect_true("knockoffs" %in% names(knockoffs))
  expect_equal(length(knockoffs$knockoffs), 2)
  expect_equal(dim(knockoffs$knockoffs[[1]]), dim(subset_X))
})
