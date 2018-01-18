\name{cv.henv}
\alias{cv.henv}
\title{Cross validation for henv}
\description{
 Compute the prediction error for the heteroscedastic envelope estimator using m-fold cross validation.
}

\usage{
cv.henv(X, Y, u, m, nperm)
}

\arguments{
  \item{X}{
  A group indicator vector of length \eqn{n}, where \eqn{n}  
  denotes the number of observations.
  }
  \item{Y}{
  Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.
  }
  \item{u}{
  Dimension of the heteroscedastic envelope. An integer between 0 and 
  r.
  }
  \item{m}{
  A positive integer that is used to indicate \code{m}-fold cross 
  validation.
  }
  \item{nperm}{
  A positive integer indicating number of permutations of the observations, m-fold cross validation is run on each permutation.
  }
}
\details{
This function computes prediction errors using m-fold cross validation.  For a fixed dimension u, the data is randomly partitioned into m parts, each part is in turn used for testing for the prediction performance while the rest m-1 parts are used for training.  This process is repeated for \code{nperm} times, and average prediction error is reported. As Y is multivariate, the identity inner product is used for computing the prediction errors.
}
\value{The output is a real nonnegative number.
\item{cvPE}{The prediction error estimated by m-fold cross validation.} 
}

\examples{
data(waterstrider)
X <- waterstrider[ , 1]
Y <- waterstrider[ , 2:5]

m <- 5
nperm <- 50

\dontrun{cvPE <- cv.henv(X, Y, 2, m, nperm)}
\dontrun{cvPE}
}
