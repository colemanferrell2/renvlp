\name{logit.envMU}
\alias{logit.envMU}
\title{Estimate the envelope subspace in logistic regression}
\description{
 Estimate the envelope subspace with specified dimension in logistic regression. 
}
\usage{
logit.envMU(X, Y, u, initial = NULL)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{Y}{Response. An n by 1 matrix. The univariate response must be binary.}
  \item{u}{Dimension of the envelope. An integer between 0 and p.}
  \item{initial}{The user-specified value of Gamma for the envelope subspace.}
}
\details{
 This function estimate the envelope subspace in logistic regression using an non-Grassmann optimization algorithm.  The starting value and optimization algorithm is described in Cook et al. (2016).
}
\value{
\item{Gammahat}{The orthonormal basis of the envelope subspace.}
\item{Gamma0hat}{The orthonormal basis of the complement of the envelope subspace.}
\item{muhat}{The estimated intercept of the canonical parameter.}
\item{etahat}{The estimated beta of the canonical parameter with respect to Gamma.}
\item{weighthat}{The estimated weight defined as C"(theta) / E(C"(theta)) where C(theta) is the conditional log likelihood.}
\item{Vhat}{The estimated V defined as V = theta + (Y - mu (theta) / W).}
\item{avar}{The asympotic covariance of vec(beta).}
\item{objfun}{The minimized objective function.}
}
\references{
Cook, R. D., Forzani, L. and Su, Z. (2016) A Note on Fast Envelope Estimation. Journal of Multivariate Analysis. 150, 42-54.
}
