# Create example data for the package
library(Matrix)

# Load the sample data
extdata_path <- system.file("extdata", package = "CoxKnockoff")

# Create small example datasets for documentation
set.seed(123)

# Example genotype matrix (20 samples, 50 SNPs)
example_genotypes <- Matrix(
  matrix(sample(0:2, 20 * 50, replace = TRUE, prob = c(0.25, 0.5, 0.25)),
         nrow = 20, ncol = 50),
  sparse = TRUE
)

# Example SNP positions
example_positions <- sort(sample(1:1000000, 50))

# Example phenotype data
example_phenotype <- data.frame(
  FID = paste0("FAM", 1:20),
  IID = paste0("IND", 1:20),
  time = round(rexp(20, rate = 1/500) + 30),
  status = rbinom(20, 1, 0.3)
)

# Example covariate data
example_covariates <- data.frame(
  FID = paste0("FAM", 1:20),
  IID = paste0("IND", 1:20),
  age = round(rnorm(20, 50, 10)),
  sex = sample(0:1, 20, replace = TRUE),
  bmi = round(rnorm(20, 25, 3), 1)
)

# Save data objects
usethis::use_data(example_genotypes, overwrite = TRUE)
usethis::use_data(example_positions, overwrite = TRUE)
usethis::use_data(example_phenotype, overwrite = TRUE)
usethis::use_data(example_covariates, overwrite = TRUE)
