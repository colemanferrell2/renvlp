\name{pred.sxenv}
\alias{pred.sxenv}
\title{Estimation or prediction for sxenv}
\description{
 Perform estimation or prediction under the scaled envelope model in the predictor space.
}
\usage{
pred.sxenv(m, Xnew)
}
\arguments{
  \item{m}{A list containing estimators and other statistics inherited from stenv.}
  \item{Xnew}{The value of X with which to estimate or predict Y.  A p dimensional vector.}
}
\details{
This function evaluates the scaled envelope model in the predictor space at new value Xnew.  It can perform estimation: find the fitted value when X = Xnew, or prediction: predict Y when X = Xnew.  The covariance matrix and the standard errors are also provided.
}
\value{The output is a list that contains following components.
\item{value}{The fitted value or the predicted value evaluated at Xnew.} 
\item{covMatrix.estm}{The covariance matrix of the fitted value at Xnew.} 
\item{SE.estm}{The standard error of the fitted value at Xnew.} 
\item{covMatrix.pred}{The covariance matrix of the predicted value at Xnew.} 
\item{SE.pred}{The standard error of the predicted value at Xnew.} 
}

\examples{
data(sales)
Y <- sales[, 1:3]
X <- sales[, 4:7]
R <- rep(1, 4)

m <- sxenv(X, Y, 2, R)

pred.res <- pred.sxenv(m, X[1, ])
pred.res

}

