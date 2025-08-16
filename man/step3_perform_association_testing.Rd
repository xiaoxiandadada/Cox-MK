\name{step3_perform_association_testing}
\alias{step3_perform_association_testing}
\alias{perform_association_testing}
\title{Step 3: Perform Association Testing}
\description{
Perform association testing with genotype data using a fitted null model.
Handles both SPACox and standard Cox regression methods.
}
\usage{
perform_association_testing(X, null_model)
}
\arguments{
\item{X}{Genotype matrix (samples × SNPs)}
\item{null_model}{Fitted null Cox model from \code{\link{step2_fit_null_cox_model}}}
}
\value{
List containing:
\item{test_stats}{Test statistics for each SNP}
\item{p_values}{P-values for each SNP}
\item{method}{Method used for testing}
\item{n_snps}{Number of SNPs tested}
}
\details{
This is Step 3 of the 4-step CoxMK workflow. It performs association testing
for both original and knockoff variables using the null model from Step 2.
}
\examples{
\dontrun{
# Perform association testing
test_results <- perform_association_testing(
  X = genotype_matrix,
  null_model = fitted_null_model
)
}
}
\seealso{
\code{\link{step1_create_knockoffs}}, \code{\link{step2_fit_null_cox_model}}, \code{\link{step4_knockoff_filter}}
}
\keyword{association}
\keyword{step3}
