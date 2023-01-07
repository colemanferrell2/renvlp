\name{cv.eppls}
\alias{cv.eppls}
\title{Cross validation for peplos}
\description{
 Compute the prediction error for the Envelope-based Partial Partial Least Squares estimator using m-fold cross validation.
}
\usage{
cv.eppls(X1, X2, Y, u, m, nperm)
}
\arguments{
  \item{X1}{Predictors of main interest. An n by p1 matrix, n is the number of observations, and p1 is the number of main predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{X2}{Covariates, or predictors not of main interest.  An n by p2 matrix, p2 is the number of covariates.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
  \item{u}{Dimension of the envelope. An integer between 0 and r.}
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
data(amitriptyline)
  
Y <- amitriptyline[ , 1:2]
X1 <- amitriptyline[ , 4:7]
X2 <- amitriptyline[ , 3]

m <- 5
nperm <- 50
\dontrun{cvPE <- cv.eppls(X1, X2, Y, 2, m, nperm)}
\dontrun{cvPE}
}

