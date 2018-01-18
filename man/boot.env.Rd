\name{boot.env}
\alias{boot.env}
\title{Bootstrap for env}
\description{
 Compute bootstrap standard error for the envelope estimator.
}
\usage{
boot.env(X, Y, u, B)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
  \item{u}{Dimension of the envelope. An integer between 0 and r.}
  \item{B}{Number of bootstrap samples. A positive integer.}
}
\details{
 This function computes the bootstrap standard errors for the regression coefficients in the envelope model by bootstrapping the residuals.
}
\value{The output is an r by p matrix.
\item{bootse}{The standard error for elements in beta computed by bootstrap.} 
}

\examples{
data(wheatprotein)
X <- wheatprotein[, 8]
Y <- wheatprotein[, 1:6]

u <- u.env(X, Y)
u

B <- 100
bootse <- boot.env(X, Y, 1, B)
bootse
}

