\name{step2_fit_null_cox_model}
\alias{step2_fit_null_cox_model}
\alias{fit_null_cox_model}
\title{Step 2: Fit Null Cox Model}
\description{
Fit a null Cox model using only covariates (no genetic variants).
This function uses SPACox when available for large-scale analysis.
}
\usage{
fit_null_cox_model(time, status, covariates = NULL)
}
\arguments{
\item{time}{Survival times}
\item{status}{Event indicators (1=event, 0=censored)}
\item{covariates}{Covariate data frame (optional)}
}
\value{
Fitted Cox model object (SPACox or coxph) with additional metadata:
\item{model_type}{Type of model fitted ("SPACox" or "Standard Cox")}
\item{formula}{Formula used for fitting}
\item{n_samples}{Number of samples}
\item{n_events}{Number of events}
}
\details{
This is Step 2 of the 4-step CoxMK workflow. The fitted null model is used
for efficient association testing in Step 3.
}
\examples{
\dontrun{
# Fit null model with covariates
null_model <- fit_null_cox_model(
  time = survival_time,
  status = event_status,
  covariates = covariate_data
)
}
}
\seealso{
\code{\link{step1_create_knockoffs}}, \code{\link{step3_perform_association_testing}}, \code{\link{step4_knockoff_filter}}
}
\keyword{cox}
\keyword{step2}
