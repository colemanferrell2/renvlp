\name{envMU}
\alias{envMU}
\title{Estimate the envelope subspace}
\description{
 Estimate the envelope subspace with specified dimension. 
}
\usage{
envMU(M, U, u)
}
\arguments{
  \item{M}{M matrix in the envelope objective function. An r by r semi-positive definite matrix.}
  \item{U}{U matrix in the envelope objective function. An r by r semi-positive definite matrix.}
  \item{u}{Dimension of the envelope. An integer between 0 and r.}
}
\details{
 This function estimate the envelope subspace using an non-Grassmann optimization algorithm.  The starting value and optimization algorithm is described in Cook et al. (2016).
}
\value{
\item{Gammahat}{The orthonormal basis of the envelope subspace.}
\item{Gamma0hat}{The orthonormal basis of the complement of the envelope subspace.}
\item{objfun}{The minimized objective function.}
}
\references{
Cook, R. D., Forzani, L. and Su, Z. (2016) A Note on Fast Envelope Estimation. Journal of Multivariate Analysis. 150, 42-54.
}
