\name{boot.penv}
\alias{boot.penv}
\title{Bootstrap for penv}
\description{
 Compute bootstrap standard error for the partial envelope estimator.
}
\usage{
boot.penv(X1, X2, Y, u, B)
}
\arguments{
  \item{X1}{Predictors of main interest. An n by p1 matrix, n is the number of observations, and p1 is the number of main predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{X2}{Covariates, or predictors not of main interest.  An n by p2 matrix, p2 is the number of covariates.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
  \item{u}{Dimension of the partial envelope. An integer between 0 and r.}
  \item{B}{Number of bootstrap samples. A positive integer.}
}
\details{
 This function computes the bootstrap standard errors for the regression coefficients beta1 in the partial envelope model by bootstrapping the residuals.
}
\value{The output is an r by p1 matrix.
\item{bootse}{The standard error for elements in beta1 computed by bootstrap.} 
}

\examples{
data(fiberpaper)
X1 <- fiberpaper[, 7]
X2 <- fiberpaper[, 5:6]
Y <- fiberpaper[, 1:4]

B <- 100
bootse <- boot.penv(X1, X2, Y, 1, B)
bootse
}

