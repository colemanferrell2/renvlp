\name{u.stenv}
\alias{u.stenv}
\title{Select the dimension of stenv}
\description{
 This function outputs dimensions selected by Akaike information criterion (AIC), Bayesian information criterion (BIC) and likelihood ratio testing with specified significance level for the simultaneous envelope model. 
}
\usage{
u.stenv(X, Y, alpha = 0.01)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
  \item{alpha}{Significance level for testing. The default is 0.01.}
}
\value{
\item{d}{Rank of beta selected by the Bura-Cook estimator.}	
\item{u.aic}{Dimension of the simultaneous envelope subspace selected by AIC.}
\item{u.bic}{Dimension of the simultaneous envelope subspace selected by BIC.}
\item{u.lrt}{Dimension of the simultaneous envelope subspace selected by the likelihood ratio testing procedure.}
\item{loglik.mat}{Log likelihood for dimension from (1, 1) to (r, p).}
\item{aic.mat}{AIC value for dimension from (1, 1) to (r, p).}
\item{bic.mat}{BIC value for dimension from (1, 1) to (r, p).}
}

\examples{
data(fiberpaper)
X <- fiberpaper[, 5:7]
Y <- fiberpaper[, 1:4]

u <- u.stenv(X, Y)
u
}

