\name{stenvMU}
\alias{stenvMU}
\title{Estimate the simultaneous envelope subspace}
\description{
 Estimate the simulatneous envelope subspace with specified dimension. 
}
\usage{
stenvMU(X, Y, q, u)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
  \item{q}{Dimension of the X-envelope. An integer between 0 and p.}
  \item{u}{Dimension of the Y-envelope. An integer between 0 and r.}
}
\details{
 This function estimate the simultaneous envelope subspace using an non-Grassmann optimization algorithm. 
}
\value{
\item{Gammahat}{The orthonormal basis of the Y-envelope subspace.}
\item{Gamma0hat}{The orthonormal basis of the complement of the Y-envelope subspace.}
\item{Phihat}{The orthonormal basis of the X-envelope subspace.}
\item{Phi0hat}{The orthonormal basis of the complement of the X-envelope subspace.}
\item{objfun}{The minimized objective function.}
}
\references{
Cook, R. D., Forzani, L. and Su, Z. (2016) A Note on Fast Envelope Estimation. Journal of Multivariate Analysis. 150, 42-54.
}
