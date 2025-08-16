\name{step1_create_knockoffs}
\alias{step1_create_knockoffs}
\alias{create_knockoffs}
\title{Step 1: Generate Knockoff Variables}
\description{
Generate knockoff variables using the Model-X approach for variable selection
in Cox regression analysis.
}
\usage{
create_knockoffs(X, pos = NULL, M = 5, method = "sdp")
}
\arguments{
\item{X}{Genotype matrix (samples × SNPs)}
\item{pos}{Genomic positions (optional)}
\item{M}{Number of knockoff copies (default: 5)}
\item{method}{Method for knockoff generation ("sdp", "equi", "asdp")}
}
\value{
List containing:
\item{knockoffs}{List of M knockoff matrices}
\item{method}{Method used}
\item{Q}{Correlation matrix (if applicable)}
}
\details{
This is Step 1 of the 4-step CoxMK workflow. It generates knockoff variables
that preserve the correlation structure of the original variables while being
conditionally independent given the original variables.
}
\examples{
\dontrun{
# Generate knockoffs for genotype data
knockoff_result <- create_knockoffs(X = genotype_matrix, M = 3)
}
}
\seealso{
\code{\link{step2_fit_null_cox_model}}, \code{\link{step3_perform_association_testing}}, \code{\link{step4_knockoff_filter}}
}
\keyword{knockoffs}
\keyword{step1}
