\name{boot.henv}
\alias{boot.henv}
\title{Bootstrap for henv}
\description{
 Compute bootstrap standard error for the heteroscedastic envelope.
}

\usage{
boot.henv(X, Y, u, B)
}
\arguments{
  \item{X}{
  A group indicator vector of length \eqn{n}, where \eqn{n} denotes the number of observations.
  }
  \item{Y}{
  Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.
  }
  \item{u}{
  Dimension of the heteroscedastic envelope. An integer between 0 and 
  r.
  }
  \item{B}{
  Number of bootstrap samples. A positive integer.
  }
}
\details{
 This function computes the bootstrap standard errors for the regression coefficients in the heteroscedastic envelope model by bootstrapping the residuals.
}
\value{The output is an r by p matrix.
\item{bootse}{The standard error for elements in beta computed by bootstrap.} 
}
\examples{
data(waterstrider)
X <- waterstrider[ , 1]
Y <- waterstrider[ , 2:5]

B <- 100
\dontrun{res <- boot.henv(X, Y, 2, B)}
\dontrun{res}

}
