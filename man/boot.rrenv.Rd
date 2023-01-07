\name{boot.rrenv}
\alias{boot.rrenv}
\title{Bootstrap for rrenv}
\description{
 Compute bootstrap standard error for the reduced rank envelope estimator.
}
\usage{
boot.rrenv(X, Y, u, d, B)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
  \item{u}{Dimension of the envelope. An integer between 0 and r.}
  \item{d}{The rank of the coefficient matrix.  An integer between 0 and u.}
  \item{B}{Number of bootstrap samples. A positive integer.}
}
\details{
 This function computes the bootstrap standard errors for the regression coefficients in the reduced rank envelope model.
}
\value{The output is an r by p matrix.
\item{bootse}{The standard error for elements in beta computed by bootstrap.} 
}

\examples{
data(vehicles)
X <- vehicles[, 1:11]  
Y <- vehicles[, 12:15]
X <- scale(X)
Y <- scale(Y)  # The scales of Y are vastly different, so scaling is reasonable here

B <- 100
\dontrun{bootse <- boot.rrenv(X, Y, 4, 2, B)}
\dontrun{bootse}
}

