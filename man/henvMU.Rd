\name{henvMU}
\alias{henvMU}
\title{Estimate the heteroscedastic envelope subspace}
\description{
 Estimate the heteroscedastic envelope subspace with specified dimension. 
}
\usage{
henvMU(M, U, MU, u, n, ng, L)
}

\arguments{
  \item{M}{
  A matrix M for the non-Grassmann manifold optimization problem in 
  Cook et al. (2016)
  }
  \item{U}{
  A matrix U for the non-Grassmann manifold optimization problem in 
  Cook et al. (2016)
  }
  \item{MU}{
  Sum of matrix M and U.
  }
  \item{u}{
  A given dimension of the heteroscedastic envelope space.  It should 
  be an interger between \eqn{0} and \eqn{r}.
  }
  \item{n}{
  The number of observations.
  }
  \item{ng}{
  A \eqn{L} by \eqn{1} vector of the number of observations in 
  each group.
  }
  \item{L}{
  The number of groups.  
  }
}

\value{
  \item{Gamma}{
  An \eqn{r} by \eqn{u} matrix for the orthonormnal basis of the  
  heteroscedastic envelope subspace.
  }
  \item{Gamma0}{
  An \eqn{r} by \eqn{(r-u)} matrix for the orthonornal basis of 
  the complement of the heteroscedastic envelope subspace.
  }
}

\references{
  Cook, R. D., Forzani, L. and Su, Z. (2016) A Note on Fast Envelope Estimation. Journal of Multivariate Analysis. 150, 42-54.
}


