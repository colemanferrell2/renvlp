\name{u.rrenv}
\alias{u.rrenv}
\title{Select the dimension of rrenv}
\description{
 This function outputs the envelope dimension selected by Akaike information criterion (AIC), Bayesian information criterion (BIC) and likelihood ratio testing with specified significance level for the reduced rank envelope model.
}
\usage{
u.rrenv(X, Y, d, alpha = 0.01)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
  \item{d}{The rank of the coefficient matrix.  An integer between 0 and r.}
  \item{alpha}{Significance level for testing.  The default is 0.01.}
}
\value{
\item{u.aic}{Dimension of the envelope subspace selected by AIC.}
\item{u.bic}{Dimension of the envelope subspace selected by BIC.}
\item{u.lrt}{Dimension of the envelope subspace selected by the likelihood ratio testing procedure.}
\item{loglik.seq}{Log likelihood for dimension from d to r.}
\item{aic.seq}{AIC value for dimension from d to r.}
\item{bic.seq}{BIC value for dimension from d to r.}
}

\examples{
data(vehicles)
X <- vehicles[, 1:11]  
Y <- vehicles[, 12:15]
X <- scale(X)
Y <- scale(Y)  # The scales of Y are vastly different, so scaling is reasonable here
d <- d.select(X, Y, 0.01)
d

\dontrun{u <- u.rrenv(X, Y, 2)}
\dontrun{u}

}

