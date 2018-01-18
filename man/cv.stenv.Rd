\name{cv.stenv}
\alias{cv.stenv}
\title{Cross validation for stenv}
\description{
 Compute the prediction error for the simultaneous envelope estimator using m-fold cross validation.
}
\usage{
cv.stenv(X, Y, q, u, m, nperm)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors and n is number of observations.  The predictors must be continuous variables.}
  \item{Y}{Responses. An n by r matrix, r is the number of responses. The response can be univariate or multivariate and must be continuous variable.}
  \item{q}{Dimension of the X-envelope. An integer between 0 and p.}
  \item{u}{Dimension of the Y-envelope. An integer between 0 and r.}
  \item{m}{A positive integer that is used to indicate m-fold cross validation.}
  \item{nperm}{A positive integer indicating number of permutations of the observations, m-fold cross validation is run on each permutation.}
}
\details{
This function computes prediction errors using m-fold cross validation.  For a fixed dimension (q, u), the data is randomly partitioned into m parts, each part is in turn used for testing for the prediction performance while the rest m-1 parts are used for training.  This process is repeated for \code{nperm} times, and average prediction error is reported. If Y is multivariate, the identity inner product is used for computing the prediction errors.
}
\value{The output is a real nonnegative number.
\item{cvPE}{The prediction error estimated by m-fold cross validation.} 
}

\examples{
data(fiberpaper)
X <- fiberpaper[, 5:7]
Y <- fiberpaper[, 1:4]

m <- 5
nperm <- 50
\dontrun{cvPE <- cv.stenv(X, Y, 2, 3, m, nperm)}
\dontrun{cvPE}
}

