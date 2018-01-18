\name{boot.stenv}
\alias{boot.stenv}
\title{Bootstrap for stenv}
\description{
 Compute bootstrap standard error for the simultaneous envelope estimator.
}
\usage{
boot.stenv(X, Y, q, u, B)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
    \item{q}{Dimension of the X-envelope. An integer between 0 and p.}
  \item{u}{Dimension of the Y-envelope. An integer between 0 and r.}
  \item{B}{Number of bootstrap samples. A positive integer.}
}
\details{
 This function computes the bootstrap standard errors for the regression coefficients in the envelope model by bootstrapping the residuals.
}
\value{The output is an p by r matrix.
\item{bootse}{The standard error for elements in beta computed by bootstrap.} 
}

\examples{
data(fiberpaper)
X <- fiberpaper[, 5:7]
Y <- fiberpaper[, 1:4]
u <- u.stenv(X, Y)
u

B <- 100
bootse <- boot.stenv(X, Y, 2, 3, B)
bootse
}
