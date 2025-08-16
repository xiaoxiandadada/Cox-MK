\name{CoxMK-package}
\alias{CoxMK-package}
\docType{package}
\title{CoxMK: Cox Regression with Model-X Knockoffs}
\description{
This package provides a complete workflow for Cox regression analysis with Model-X knockoffs
for variable selection in survival analysis.

The workflow follows four main steps:
\enumerate{
  \item Generate knockoff variables using \code{\link{create_knockoffs}}
  \item Fit null Cox model using \code{\link{fit_null_cox_model}}
  \item Perform association testing using \code{\link{perform_association_testing}}
  \item Apply knockoff filter using \code{\link{knockoff_filter}}
}
}
\section{Main Functions}{
\itemize{
  \item \code{\link{create_knockoffs}} - Step 1: Generate knockoff variables
  \item \code{\link{fit_null_cox_model}} - Step 2: Fit null Cox model
  \item \code{\link{perform_association_testing}} - Step 3: Association testing
  \item \code{\link{knockoff_filter}} - Step 4: Variable selection
}
}
\author{CoxMK Development Team}
\keyword{package}
\keyword{survival}
\keyword{knockoffs}
