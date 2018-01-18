\name{u.env}
\alias{u.env}
\title{Select the dimension of env}
\description{
 This function outputs dimensions selected by Akaike information criterion (AIC), Bayesian information criterion (BIC) and likelihood ratio testing with specified significance level for the envelope model.
}
\usage{
u.env(X, Y, alpha = 0.01)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
  \item{alpha}{Significance level for testing.  The default is 0.01.}
}
\value{
\item{u.aic}{Dimension of the envelope subspace selected by AIC.}
\item{u.bic}{Dimension of the envelope subspace selected by BIC.}
\item{u.lrt}{Dimension of the envelope subspace selected by the likelihood ratio testing procedure.}
\item{loglik.seq}{Log likelihood for dimension from 0 to r.}
\item{aic.seq}{AIC value for dimension from 0 to r.}
\item{bic.seq}{BIC value for dimension from 0 to r.}
}

\examples{
data(wheatprotein)
X <- wheatprotein[, 8]
Y <- wheatprotein[, 1:6]
u <- u.env(X, Y)
u
}

