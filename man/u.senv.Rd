\name{u.senv}
\alias{u.senv}
\title{Select the dimension of senv}
\description{
 This function outputs dimensions selected by Akaike information criterion (AIC), Bayesian information criterion (BIC) for the scaled response envelope model.
}
\usage{
u.senv(X, Y)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
}
\value{
\item{u.aic}{Dimension of the scaled envelope subspace selected by AIC.}
\item{u.bic}{Dimension of the scaled envelope subspace selected by BIC.}
\item{aic.seq}{AIC value for dimension from 0 to r.}
\item{bic.seq}{BIC value for dimension from 0 to r.}
}

\examples{
data(sales)
X <- sales[, 1:3]
Y <- sales[, 4:7]

u <- u.senv(X, Y)
u
}

