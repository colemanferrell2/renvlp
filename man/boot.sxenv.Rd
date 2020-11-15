\name{boot.sxenv}
\alias{boot.sxenv}
\title{Bootstrap for sxenv}
\description{
 Compute bootstrap standard error for the scaled predictor envelope estimator.
}
\usage{
boot.sxenv(X, Y, u, R, B)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
  \item{u}{Dimension of the scaled envelope in the predictor space. An integer between 0 and p.}
  \item{R}{The number of replications of the scales. A vector, the sum of all elements of R must be p.}
  \item{B}{Number of bootstrap samples. A positive integer.}
}
\details{
 This function computes the bootstrap standard errors for the regression coefficients in the scaled envelope model in the predictor space by bootstrapping the residuals.
}
\value{The output is an p by r matrix.
\item{bootse}{The standard error for elements in beta computed by bootstrap.} 
}

\examples{
data(sales)
Y <- sales[, 1:3]
X <- sales[, 4:7]
R <- rep(1, 4)

u <- u.sxenv(X, Y, R)
u

B <- 100
\dontrun{bootse <- boot.sxenv(X, Y, 2, R, B)}
\dontrun{bootse}
}
