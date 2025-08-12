#' Example Genotype Data
#'
#' A sparse matrix containing example genotype data for 20 samples and 50 SNPs.
#' Values are coded as 0, 1, 2 representing the number of minor alleles.
#'
#' @format A 20 x 50 sparse matrix (dgCMatrix)
#' @source Simulated data for package examples
"example_genotypes"

#' Example SNP Positions
#'
#' A vector of genomic positions (in base pairs) for the example SNPs.
#'
#' @format A numeric vector of length 50
#' @source Simulated data for package examples
"example_positions"

#' Example Phenotype Data
#'
#' Time-to-event phenotype data for the example samples.
#'
#' @format A data frame with 20 rows and 4 variables:
#' \describe{
#'   \item{FID}{Family ID}
#'   \item{IID}{Individual ID}
#'   \item{time}{Survival time in days}
#'   \item{status}{Event status (0 = censored, 1 = event)}
#' }
#' @source Simulated data for package examples
"example_phenotype"

#' Example Covariate Data
#'
#' Covariate data for the example samples.
#'
#' @format A data frame with 20 rows and 5 variables:
#' \describe{
#'   \item{FID}{Family ID}
#'   \item{IID}{Individual ID}
#'   \item{age}{Age in years}
#'   \item{sex}{Sex (0 = male, 1 = female)}
#'   \item{bmi}{Body mass index}
#' }
#' @source Simulated data for package examples
"example_covariates"
