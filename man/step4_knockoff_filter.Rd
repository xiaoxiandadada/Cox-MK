\name{step4_knockoff_filter}
\alias{step4_knockoff_filter}
\alias{knockoff_filter}
\title{Step 4: Apply Knockoff Filter}
\description{
Apply the knockoff filter to select variables while controlling the
false discovery rate (FDR) at a specified level.
}
\usage{
knockoff_filter(W, fdr = 0.1, offset = 1)
}
\arguments{
\item{W}{Vector of W statistics computed from original and knockoff test statistics}
\item{fdr}{Target false discovery rate (default: 0.1)}
\item{offset}{Offset parameter for knockoff filter (default: 1)}
}
\value{
Vector of indices of selected variables. The \code{threshold} attribute
contains the threshold value used for selection.
}
\details{
This is Step 4 of the 4-step CoxMK workflow. It applies the knockoff filter
to control FDR while selecting truly associated variables.

The W statistics should be computed by comparing test statistics from Step 3
for original variables with those from their knockoff counterparts.
}
\examples{
\dontrun{
# Apply knockoff filter
selected_vars <- knockoff_filter(W_statistics, fdr = 0.1)

# Check threshold used
threshold <- attr(selected_vars, "threshold")
}
}
\seealso{
\code{\link{step1_create_knockoffs}}, \code{\link{step2_fit_null_cox_model}}, \code{\link{step3_perform_association_testing}}
}
\keyword{selection}
\keyword{step4}
