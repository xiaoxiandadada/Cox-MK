\name{CoxMK-workflow}
\alias{CoxMK-workflow}
\title{CoxMK: Four-Step Cox Regression with Model-X Knockoffs Workflow}
\description{
Complete workflow for Cox regression analysis with Model-X knockoffs for survival analysis.
This package implements a four-step methodology for variable selection with false discovery rate control.
}

\section{Overview}{
The CoxMK workflow consists of four main steps:

\subsection{Step 1: Generate Knockoff Variables}{
\strong{Main Function:} \code{\link{create_knockoffs}}

\strong{Supporting Functions:}
\itemize{
  \item \code{\link{load_plink_data}} - Load genotype data from PLINK format
}

\strong{Purpose:} Generate knockoff variables that preserve the correlation structure 
of the original variables while being conditionally independent.

\strong{Input:} Genotype matrix or PLINK files
\strong{Output:} Knockoff matrices stored in GDS format
}

\subsection{Step 2: Fit Null Cox Model}{
\strong{Main Function:} \code{\link{fit_null_cox_model}}

\strong{Supporting Functions:}
\itemize{
  \item \code{\link{prepare_phenotype}} - Prepare survival time and status data
  \item \code{\link{load_covariates}} - Load covariate data
}

\strong{Purpose:} Fit null Cox proportional hazards model using SPACox for efficient 
large-scale analysis.

\strong{Input:} Phenotype data (time, status) and covariates
\strong{Output:} Fitted null model object
}

\subsection{Step 3: Perform Association Testing}{
\strong{Main Function:} \code{\link{perform_association_testing}}

\strong{Purpose:} Conduct association testing using both original and knockoff variables
with the fitted null model.

\strong{Input:} Original and knockoff genotype data, fitted null model
\strong{Output:} Test statistics for original and knockoff variables
}

\subsection{Step 4: Apply Knockoff Filter}{
\strong{Main Function:} \code{\link{knockoff_filter}}

\strong{Supporting Functions:}
\itemize{
  \item \code{\link{calculate_w_statistics}} - Calculate W statistics
}

\strong{Purpose:} Apply the knockoff filter to select significant variables 
with false discovery rate control.

\strong{Input:} Test statistics from original and knockoff variables
\strong{Output:} Selected variables with FDR control
}
}

\section{Complete Workflow Function}{
\strong{One-Step Analysis:} \code{\link{cox_knockoff_analysis}}

For users who want to run the complete analysis in one step, the 
\code{cox_knockoff_analysis} function provides a wrapper that executes 
all four steps automatically.
}

\section{Example Workflow}{
\preformatted{
# Step-by-step analysis
library(CoxMK)

# Step 1: Generate knockoffs
knockoffs <- create_knockoffs(X, M = 5)

# Step 2: Fit null model  
null_model <- fit_null_cox_model(time, status, covariates)

# Step 3: Perform association testing
results <- perform_association_testing(X_orig, X_knockoffs, null_model)

# Step 4: Apply knockoff filter
selected <- knockoff_filter(results$W_stats, fdr = 0.05)

# Alternative: Complete analysis in one step
results <- cox_knockoff_analysis(
  plink_prefix = "data/genotypes",
  time = survival_time,
  status = event_status,
  covariates = covar_data,
  fdr = 0.05
)
}
}

\seealso{
\strong{Step 1:} \code{\link{create_knockoffs}}

\strong{Step 2:} \code{\link{fit_null_cox_model}}

\strong{Step 3:} \code{\link{perform_association_testing}}

\strong{Step 4:} \code{\link{knockoff_filter}}

\strong{Complete Workflow:} \code{\link{cox_knockoff_analysis}}
}

\keyword{survival}
\keyword{knockoffs}
\keyword{variable selection}
