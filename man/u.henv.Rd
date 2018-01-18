\name{u.henv}
\alias{u.henv}
\title{Select the dimension of henv}
\description{
 This function outputs dimensions selected by Akaike information criterion (AIC), Bayesian information criterion (BIC) and likelihood ratio testing with specified significance level for the heteroscedastic envelope model.
}

\usage{
u.henv(X, Y, alpha = 0.01)
}

\arguments{
  \item{X}{
  A group indicator vector of length \eqn{n}, where \eqn{n} 
  denotes the number of observations.
  }
  \item{Y}{
  Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.
  }
  \item{alpha}{
  Significance level for testing.  The default is 0.01.
  }
}
\value{
\item{u.aic}{Dimension of the heteroscedastic envelope subspace selected by AIC.}
\item{u.bic}{Dimension of the heteroscedastic envelope subspace selected by BIC.}
\item{u.lrt}{Dimension of the heteroscedastic envelope subspace selected by the likelihood ratio testing procedure.}
\item{loglik.seq}{Log likelihood for dimension from 0 to r.}
\item{aic.seq}{AIC value for dimension from 0 to r.}
\item{bic.seq}{BIC value for dimension from 0 to r.}
}

\examples{
data(waterstrider)
X <- waterstrider[ , 1]
Y <- waterstrider[ , 2:5]

\dontrun{u <- u.henv(X, Y)}
\dontrun{u}
}

