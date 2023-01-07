\name{pred.rrenv}
\alias{pred.rrenv}
\title{Estimation or prediction for rrenv}
\description{
 Perform estimation or prediction under the reduced rank envelope model.
}
\usage{
pred.rrenv(m, Xnew)
}
\arguments{
  \item{m}{A list containing estimators and other statistics inherited from rrenv.}
  \item{Xnew}{The value of X with which to estimate or predict Y.  A p dimensional vector.}
}
\details{
This function evaluates the reduced rank envelope model at new value Xnew.  It can perform estimation: find the fitted value when X = Xnew, or prediction: predict Y when X = Xnew.  The covariance matrix and the standard errors are also provided.
}
\value{The output is a list that contains following components.
\item{value}{The fitted value or the predicted value evaluated at Xnew.} 
\item{covMatrix.estm}{The covariance matrix of the fitted value at Xnew.} 
\item{SE.estm}{The standard error of the fitted value at Xnew.} 
\item{covMatrix.pred}{The covariance matrix of the predicted value at Xnew.} 
\item{SE.pred}{The standard error of the predicted value at Xnew.}
}

\examples{
data(vehicles)
X <- vehicles[, 1:11]  
Y <- vehicles[, 12:15]
X <- scale(X)
Y <- scale(Y)  # The scales of Y are vastly different, so scaling is reasonable here

m <- rrenv(X, Y, 4, 2)
m

X <- as.matrix(X)
pred.res <- pred.rrenv(m, X[2, ])
pred.res

}

