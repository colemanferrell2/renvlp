\name{boot.senv}
\alias{boot.senv}
\title{Bootstrap for senv}
\description{
 Compute bootstrap standard error for the scaled envelope estimator.
}
\usage{
boot.senv(X, Y, u, B)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
  \item{u}{Dimension of the scaled envelope. An integer between 0 and r.}
  \item{B}{Number of bootstrap samples. A positive integer.}
}
\details{
 This function computes the bootstrap standard errors for the regression coefficients in the scaled envelope model by bootstrapping the residuals.
}
\value{The output is an r by p matrix.
\item{bootse}{The standard error for elements in beta computed by bootstrap.} 
}

\examples{
data(sales)
X <- sales[, 1:3]
Y <- sales[, 4:7]

u <- u.senv(X, Y)
u

\dontrun{B <- 100}
\dontrun{bootse <- boot.senv(X, Y, 2, B)}
\dontrun{bootse}
}

