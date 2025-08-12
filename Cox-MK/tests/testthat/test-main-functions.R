library(testthat)
library(CoxMK)

test_that("create_knockoffs works with example data and different M values", {
  data(example_genotypes)
  data(example_positions)
  
  # Test with M = 2
  knockoffs_2 <- create_knockoffs(
    X = example_genotypes,
    pos = example_positions,
    M = 2,
    corr_max = 0.75
  )
  
  expect_equal(length(knockoffs_2), 2)
  expect_equal(dim(knockoffs_2[[1]]), dim(example_genotypes))
  expect_equal(dim(knockoffs_2[[2]]), dim(example_genotypes))
  
  # Test with default M = 5
  knockoffs_5 <- create_knockoffs(
    X = example_genotypes,
    pos = example_positions,
    M = 5,
    corr_max = 0.75
  )
  
  expect_equal(length(knockoffs_5), 5)
  
  # Test with M = 10
  knockoffs_10 <- create_knockoffs(
    X = example_genotypes,
    pos = example_positions,
    M = 10,
    corr_max = 0.75
  )
  
  expect_equal(length(knockoffs_10), 10)
})

test_that("fit_null_model works with SPACox option", {
  data(example_phenotype)
  data(example_covariates)
  
  pheno_data <- merge(example_phenotype, example_covariates, by = c("FID", "IID"))
  covariates <- pheno_data[, c("age", "sex", "bmi")]
  
  # Test standard coxph
  model <- fit_null_model(
    time = pheno_data$time,
    status = pheno_data$status,
    covariates = covariates,
    use_spacox = FALSE
  )
  
  expect_s3_class(model, "coxph")
  
  # Test SPACox option (should fall back to coxph if SPACox not available)
  model_spa <- fit_null_model(
    time = pheno_data$time,
    status = pheno_data$status,
    covariates = covariates,
    use_spacox = TRUE
  )
  
  expect_s3_class(model_spa, "coxph")
})

test_that("calculate_w_statistics works", {
  p <- 10
  original_pvals <- runif(p)
  original_coefs <- rnorm(p)
  knockoff_pvals <- list(runif(p), runif(p))
  knockoff_coefs <- list(rnorm(p), rnorm(p))
  
  w_stats <- calculate_w_statistics(
    original_pvals = original_pvals,
    knockoff_pvals = knockoff_pvals,
    original_coefs = original_coefs,
    knockoff_coefs = knockoff_coefs
  )
  
  expect_equal(length(w_stats), p)
  expect_true(all(is.finite(w_stats)))
})

test_that("knockoff_filter works", {
  w_stats <- c(3, 1, -2, 0.5, -1, 2)
  
  result <- knockoff_filter(w_stats, fdr = 0.2)
  
  expect_is(result, "list")
  expect_true("selected" %in% names(result))
  expect_true("threshold" %in% names(result))
  expect_true("fdp" %in% names(result))
  expect_true(all(result$selected %in% seq_along(w_stats)))
})
