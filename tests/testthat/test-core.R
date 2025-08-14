library(testthat)
library(CoxMK)

test_that("Package loads correctly", {
  expect_true("CoxMK" %in% loadedNamespaces())
})

test_that("Example data loads correctly", {
  data(example_genotypes)
  data(example_positions)
  data(example_phenotype)
  data(example_covariates)
  
  expect_true(is.matrix(example_genotypes))
  expect_true(is.numeric(example_positions))
  expect_true(is.data.frame(example_phenotype))
  expect_true(is.data.frame(example_covariates))
})

test_that("Core functions work", {
  data(example_genotypes)
  data(example_positions)
  
  # Test knockoff creation (simple test)
  knockoffs <- create_knockoffs(
    X = example_genotypes[1:10, 1:5],  # Small subset for testing
    pos = example_positions[1:5],
    M = 2
  )
  
  expect_equal(length(knockoffs), 2)
  expect_equal(dim(knockoffs[[1]]), c(10, 5))
})
