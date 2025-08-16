\name{CoxMK-package}
\alias{CoxMK-package}
\docType{package}
\title{CoxMK: Cox Regression with Multiple Knockoffs}
\description{
This package provides a complete workflow for Cox regression analysis with Multiple knockoffs
for variable selection in survival analysis.

The workflow follows four main steps:
\enumerate{
  \item Generate multiple knockoff variables
  \item Fit null Cox model  
  \item Perform association testing
  \item Apply knockoff filter for variable selection
}
}
\section{Main Functions}{
\itemize{
  \item \code{\link{cox_knockoff_analysis}} - Complete analysis workflow
  \item \code{\link{create_knockoffs}} - Generate multiple knockoff variables  
  \item \code{\link{knockoff_filter}} - Apply knockoff filter for variable selection
  \item \code{\link{calculate_w_statistics}} - Calculate W statistics
}
}
\author{CoxMK Development Team}
\keyword{package}
\keyword{survival}
\keyword{knockoffs}
