\name{boot.genv}
\alias{boot.genv}
\title{Bootstrap for genv}
\description{
 Compute bootstrap standard error for the groupwise envelope.
}

\usage{
boot.genv(X, Y, Z, u, B)
}
\arguments{
  \item{X}{
  Predictors. An n by p matrix, p is the number of predictors. The 
  predictors can be univariate or multivariate, discrete or 
  continuous.  
  }
  \item{Y}{
  Multivariate responses. An n by r matrix, r is the number of 
  responses and n is number of observations. The responses must be 
  continuous variables.
  }
  \item{Z}{
   A group indicator vector of length \eqn{n}, where \eqn{n} 
  denotes the number of observations.	
  }
  \item{u}{
  Dimension of the groupwise envelope. An integer between 0 and r.
  }
  \item{B}{
  Number of bootstrap samples. A positive integer.
  }
}
\details{
 This function computes the bootstrap standard errors for the regression coefficients in the groupwise envelope model by bootstrapping the residuals.
}
\value{The output is an p by r matrix.
\item{bootse}{The standard error for elements in beta computed by bootstrap.} 
}
\examples{
data(fiberpaper)
X <- fiberpaper[ , c(5, 7)]
Y <- fiberpaper[ , 1:3]
Z <- as.numeric(fiberpaper[ , 6] > mean(fiberpaper[ , 6]))

B <- 100
\dontrun{res <- boot.genv(X, Y, Z, 2, B)}
\dontrun{res$bootse[[1]]}
\dontrun{res$bootse[[2]]}

}
