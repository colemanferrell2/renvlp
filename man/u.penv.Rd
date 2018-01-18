\name{u.penv}
\alias{u.penv}
\title{Select the dimension of penv}
\description{
 This function outputs dimensions selected by Akaike information criterion (AIC), Bayesian information criterion (BIC) and likelihood ratio testing with specified significance level for the partial envelope model.
}
\usage{
u.penv(X1, X2, Y, alpha = 0.01)
}
\arguments{
  \item{X1}{Predictors of main interest. An n by p1 matrix, n is the number of observations, and p1 is the number of main predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{X2}{Covariates, or predictors not of main interest.  An n by p2 matrix, p2 is the number of covariates.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
  \item{alpha}{Significance level for testing.  The default is 0.01.}
}
\value{
\item{u.aic}{Dimension of the partial envelope subspace selected by AIC.}
\item{u.bic}{Dimension of the partial envelope subspace selected by BIC.}
\item{u.lrt}{Dimension of the partial envelope subspace selected by the likelihood ratio testing procedure.}
\item{loglik.seq}{Log likelihood for dimension from 0 to r.}
\item{aic.seq}{AIC value for dimension from 0 to r.}
\item{bic.seq}{BIC value for dimension from 0 to r.}
}

\examples{
data(fiberpaper)
X1 <- fiberpaper[, 7]
X2 <- fiberpaper[, 5:6]
Y <- fiberpaper[, 1:4]

u <- u.penv(X1, X2, Y)
u
}

