\name{u.sxenv}
\alias{u.sxenv}
\title{Select the dimension of sxenv}
\description{
 This function outputs dimensions selected by Akaike information criterion (AIC), Bayesian information criterion (BIC) for the scaled envelope model in the predictor space. 
}
\usage{
u.sxenv(X, Y, R)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
  \item{R}{The number of replications of the scales. A vector, the sum of all elements of R must be p.}
}
\value{
\item{u.aic}{Dimension of the scaled envelope subspace in the predictor space selected by AIC.}
\item{u.bic}{Dimension of the scaled envelope subspace in the predictor space selected by BIC.}
\item{aic.seq}{AIC value for dimension from 0 to p.}
\item{bic.seq}{BIC value for dimension from 0 to p.}
}

\examples{
data(sales)
Y <- sales[, 1:3]
X <- sales[, 4:7]
R <- rep(1, 4)

u <- u.sxenv(X, Y, R)
u
}

