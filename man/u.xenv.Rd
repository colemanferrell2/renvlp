\name{u.xenv}
\alias{u.xenv}
\title{Select the dimension of xenv}
\description{
This function outputs dimensions selected by Akaike information criterion (AIC), Bayesian information criterion (BIC) and likelihood ratio testing with specified significance level for the envelope model.
}
\usage{
u.xenv(X, Y, alpha = 0.01)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors and n is number of observations.  The predictors must be continuous variables.}
  \item{Y}{Responses. An n by r matrix, r is the number of responses. The response can be univariate or multivariate and must be continuous variable.}
  \item{alpha}{Significance level for testing.  The default is 0.01.}
}
\value{
\item{u.aic}{Dimension of the envelope subspace selected by AIC.}
\item{u.bic}{Dimension of the envelope subspace selected by BIC.}
\item{u.lrt}{Dimension of the envelope subspace selected by the likelihood ratio testing procedure.}
\item{loglik.seq}{Log likelihood for dimension from 0 to p.}
\item{aic.seq}{AIC value for dimension from 0 to p.}
\item{bic.seq}{BIC value for dimension from 0 to p.}
}

\examples{
data(wheatprotein)
X <- wheatprotein[, 1:6]
Y <- wheatprotein[, 7]

u <- u.xenv(X, Y)
u
}

