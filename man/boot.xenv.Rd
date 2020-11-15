\name{boot.xenv}
\alias{boot.xenv}
\title{Bootstrap for xenv}
\description{
 Compute bootstrap standard error for the predictor envelope estimator.
}
\usage{
boot.xenv(X, Y, u, B)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors and n is number of observations.  The predictors must be continuous variables.}
  \item{Y}{Responses. An n by r matrix, r is the number of responses. The response can be univariate or multivariate and must be continuous variable.}
  \item{u}{Dimension of the envelope. An integer between 0 and p.}
  \item{B}{Number of bootstrap samples. A positive integer.}
}
\details{
 This function computes the bootstrap standard errors for the regression coefficients in the envelope model in predictor space by bootstrapping the residuals.
}
\value{The output is a p by r matrix.
\item{bootse}{The standard error for elements in beta computed by bootstrap.} 
}

\examples{
data(wheatprotein)
X <- wheatprotein[, 1:6]
Y <- wheatprotein[, 7]

\dontrun{B <- 100}
\dontrun{bootse <- boot.xenv(X, Y, 2, B)}
\dontrun{bootse}
}

