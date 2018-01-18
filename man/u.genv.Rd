\name{u.genv}
\alias{u.genv}
\title{Select the dimension of genv}
\description{
 This function outputs dimensions selected by Akaike information criterion (AIC), Bayesian information criterion (BIC) and likelihood ratio testing with specified significance level for the groupwise envelope model.
}

\usage{
u.genv(X, Y, Z, alpha = 0.01)
}

\arguments{
  \item{X}{
  Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.
  }
  \item{Y}{
  Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.
  }
  \item{Z}{
   A group indicator vector of length \eqn{n}, where \eqn{n} 
  denotes the number of observations.
  }
  \item{alpha}{
  Significance level for testing.  The default is 0.01.
  }
}
\value{
\item{u.aic}{Dimension of the groupwise envelope subspace selected by AIC.}
\item{u.bic}{Dimension of the groupwise envelope subspace selected by BIC.}
\item{u.lrt}{Dimension of the groupwise envelope subspace selected by the likelihood ratio testing procedure.}
\item{loglik.seq}{Log likelihood for dimension from 0 to r.}
\item{aic.seq}{AIC value for dimension from 0 to r.}
\item{bic.seq}{BIC value for dimension from 0 to r.}
}

\examples{
data(fiberpaper)
X <- fiberpaper[ , c(5, 7)]
Y <- fiberpaper[ , 1:3]
Z <- as.numeric(fiberpaper[ , 6] > mean(fiberpaper[ , 6]))

u <- u.genv(X, Y, Z)
u
}

