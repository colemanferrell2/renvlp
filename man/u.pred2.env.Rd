\name{u.pred2.env}
\alias{u.pred2.env}
\title{Select the dimension of the constructed partial envelope for prediction based on envelope model}
\description{
 This function outputs dimensions selected by Akaike information criterion (AIC), Bayesian information criterion (BIC) and likelihood ratio testing with specified significance level for the constructed partial envelope model.
}
\usage{
u.pred2.env(X, Y, Xnew, alpha = 0.01)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
  \item{Xnew}{The value of X with which to estimate or predict Y.  A p dimensional vector.}
  \item{alpha}{Significance level for testing.  The default is 0.01.}
}

\value{
\item{u.aic}{Dimension of the constructed partial envelope subspace selected by AIC.}
\item{u.bic}{Dimension of the constructed partial envelope subspace selected by BIC.}
\item{u.lrt}{Dimension of the constructed partial envelope subspace selected by the likelihood ratio testing procedure.}
\item{loglik.seq}{Log likelihood for dimension from 0 to r.}
\item{aic.seq}{AIC value for dimension from 0 to r.}
\item{bic.seq}{BIC value for dimension from 0 to r.}
}

\examples{
data(fiberpaper)
X <- fiberpaper[, 5:7]
Y <- fiberpaper[, 1:4]

u <- u.pred2.env(X, Y, X[10, ])
u
}

