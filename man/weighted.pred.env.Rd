\name{weighted.pred.env}
\alias{weighted.pred.env}
\title{Estimation or prediction using weighted partial envelope}
\description{
 Perform estimation or prediction through weighted partial envelope model.
}
\usage{
weighted.pred.env(X, Y, Xnew)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
  \item{Xnew}{The value of X with which to estimate or predict Y.  A p dimensional vector.}
}
\details{
This function evaluates the envelope model at new value Xnew.  It can perform estimation: find the fitted value when X = Xnew, or prediction: predict Y when X = Xnew.  But it does not provide the estimation or prediction error.  This function performs prediction using the same procedure as in pred2.env, except that the partial envelope estimator with dimension u is replaced by a weighted partial envelope estimator.  The weights are decided based on BIC values.  
}
\value{
\item{value}{The fitted value or the predicted value evaluated at Xnew.} 
}

\examples{
data(fiberpaper)
X <- fiberpaper[, 5:7]
Y <- fiberpaper[, 1:4]

\dontrun{pred.res <- weighted.pred.env(X, Y, X[10, ])}

}

