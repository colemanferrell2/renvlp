\name{cv.xenv}
\alias{cv.xenv}
\title{Cross validation for xenv}
\description{
 Compute the prediction error for the predictor envelope estimator using m-fold cross validation.
}
\usage{
cv.xenv(X, Y, u, m, nperm)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors and n is number of observations.  The predictors must be continuous variables.}
  \item{Y}{Responses. An n by r matrix, r is the number of responses. The response can be univariate or multivariate and must be continuous variable.}
  \item{u}{Dimension of the envelope. An integer between 0 and p.}
  \item{m}{A positive integer that is used to indicate m-fold cross validation.}
  \item{nperm}{A positive integer indicating number of permutations of the observations, m-fold cross validation is run on each permutation.}
}
\details{
This function computes prediction errors using m-fold cross validation.  For a fixed dimension u, the data is randomly partitioned into m parts, each part is in turn used for testing for the prediction performance while the rest m-1 parts are used for training.  This process is repeated for \code{nperm} times, and average prediction error is reported. If Y is multivariate, the identity inner product is used for computing the prediction errors.
}
\value{The output is a real nonnegative number.
\item{cvPE}{The prediction error estimated by m-fold cross validation.} 
}

\examples{
data(wheatprotein)
X <- wheatprotein[, 1:6]
Y <- wheatprotein[, 7]

m <- 5
nperm <- 50
\dontrun{cvPE <- cv.xenv(X, Y, 2, m, nperm)}
\dontrun{cvPE}
}

