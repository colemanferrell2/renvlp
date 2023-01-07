\name{cv.env.tcond}
\alias{cv.env.tcond}
\title{Cross validation for env.tcond}
\description{
 Compute the prediction error using m-fold cross validation for the response envelope estimator where the errors follow a t-distribution.
}
\usage{
cv.env.tcond(X, Y, u, df, m, nperm)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
  \item{u}{Dimension of the envelope. An integer between 0 and r.}
  \item{df}{Degrees of freedom of the t-distribution. A positive number that is greater than 2.}
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
data(concrete)
X <- concrete[1:78, 1:7]  # The first 78 observations are training data
Y <- concrete[1:78, 8:10]
\dontrun{u <- u.env.tcond(X, Y, 6)}
\dontrun{u}

m <- 5
nperm <- 50
\dontrun{cvPE <- cv.env.tcond(X, Y, 2, 6, m, nperm)}
\dontrun{cvPE}
}

