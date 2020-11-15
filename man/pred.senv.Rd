\name{pred.senv}
\alias{pred.senv}
\title{Estimation or prediction for senv}
\description{
 Perform estimation or prediction under the scaled response envelope model.
}
\usage{
pred.senv(m, Xnew)
}
\arguments{
  \item{m}{A list containing estimators and other statistics inherited from scale.env.}
  \item{Xnew}{The value of X with which to estimate or predict Y.  A p dimensional vector.}
}
\details{
This function evaluates the scaled envelope model at new value Xnew.  It can perform estimation: find the fitted value when X = Xnew, or prediction: predict Y when X = Xnew.  The covariance matrix and the standard errors are also provided.
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
X <- sales[, 1:3]
Y <- sales[, 4:7]

m <- senv(X, Y, 2)

pred.res <- pred.senv(m, X[2, ])
pred.res

}

