\name{cv.rrenv.apweights}
\alias{cv.rrenv.apweights}
\title{Cross validation for rrenv.apweights}
\description{
 Compute the prediction error using m-fold cross validation for the reduced rank envelope estimator that accommodates nonconstant error covariance.
}
\usage{
cv.rrenv.apweights(X, Y, u, d, m, nperm)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
  \item{u}{Dimension of the envelope. An integer between 0 and r.}
  \item{d}{The rank of the coefficient matrix.  An integer between 0 and u.}
  \item{m}{A positive integer that is used to indicate m-fold cross validation.}
  \item{nperm}{A positive integer indicating number of permutations of the observations, m-fold cross validation is run on each permutation.}
}
\details{
This function computes prediction errors using m-fold cross validation.  For a fixed dimension u, the data is randomly partitioned into m parts, each part is in turn used for testing for the prediction performance while the rest m-1 parts are used for training.  This process is repeated for \code{nperm} times, and average prediction error is reported. As Y is multivariate, the identity inner product is used for computing the prediction errors.
}
\value{The output is a real nonnegative number.
\item{cvPE}{The prediction error estimated by m-fold cross validation.} 
}

\examples{
data(vehicles)
X <- vehicles[, 1:11]  
Y <- vehicles[, 12:15]
X <- scale(X)
Y <- scale(Y)  # The scales of Y are vastly different, so scaling is reasonable here

m <- 5
nperm <- 50
\dontrun{cvPE <- cv.rrenv.apweights(X, Y, 3, 2, m, nperm)}
\dontrun{cvPE}
}

