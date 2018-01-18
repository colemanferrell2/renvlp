\name{senvMU}
\alias{senvMU}
\title{Estimate the scaled envelope subspace}
\description{
 Estimate the scaled envelope subspace with specified dimension. 
}
\usage{
senvMU(X, Y, u)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
  \item{u}{Dimension of the envelope. An integer between 0 and r.}
}
\details{
 This function estimate the scaled envelope subspace using an non-Grassmann optimization algorithm and nonlinear optimization using augmented Lagrange method. 
}
\value{
\item{Gammahat}{The orthonormal basis of the scaled envelope subspace.}
\item{Gamma0hat}{The orthonormal basis of the complement of the scaled envelope subspace.}
\item{Lambdahat}{The matrix of estimated scales.}
\item{objfun}{The minimized objective function.}
}
\references{
Cook, R. D., Forzani, L. and Su, Z. (2016) A Note on Fast Envelope Estimation. Journal of Multivariate Analysis. 150, 42-54.

Ye, Y., Interior algorithms for linear, quadratic, and linearly constrained non linear programming, PhD Thesis, Departments of EES stanford University, Stanford CA. 
}

