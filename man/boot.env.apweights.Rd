\name{boot.env.apweights}
\alias{boot.env.apweights}
\title{Bootstrap for env.apweights}
\description{
 Compute bootstrap standard error for the response envelope estimator that accommodates nonconstant variance.
}
\usage{
boot.env.apweights(X, Y, u, B)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
  \item{u}{Dimension of the envelope. An integer between 0 and r.}
  \item{B}{Number of bootstrap samples. A positive integer.}
}
\details{
 This function computes the bootstrap standard errors for the regression coefficients in the envelope model with nonconstant variance by bootstrapping the residuals.
}
\value{The output is an r by p matrix.
\item{bootse}{The standard error for elements in beta computed by bootstrap.} 
}

\examples{
data(concrete)
X <- concrete[, 1:7]
Y <- concrete[, 8:10]

\dontrun{u <- u.env.apweights(X, Y)}
\dontrun{u}

B <- 100
\dontrun{bootse <- boot.env.apweights(X, Y, 1, B)}
\dontrun{bootse}
}

