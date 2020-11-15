\name{pred.env}
\alias{pred.env}
\title{Estimation or prediction for env}
\description{
 Perform estimation or prediction under the response envelope model.
}
\usage{
pred.env(m, Xnew)
}
\arguments{
  \item{m}{A list containing estimators and other statistics inherited from env.}
  \item{Xnew}{The value of X with which to estimate or predict Y.  A p dimensional vector.}
}
\details{
This function evaluates the envelope model at new value Xnew.  It can perform estimation: find the fitted value when X = Xnew, or prediction: predict Y when X = Xnew.  The covariance matrix and the standard errors are also provided.
}
\value{The output is a list that contains following components.
\item{value}{The fitted value or the predicted value evaluated at Xnew.} 
\item{covMatrix.estm}{The covariance matrix of the fitted value at Xnew.} 
\item{SE.estm}{The standard error of the fitted value at Xnew.} 
\item{covMatrix.pred}{The covariance matrix of the predicted value at Xnew.} 
\item{SE.pred}{The standard error of the predicted value at Xnew.}
}

\examples{
data(wheatprotein)
X <- wheatprotein[, 8]
Y <- wheatprotein[, 1:6]
u <- u.env(X, Y)
u

m <- env(X, Y, 1)
m

X <- as.matrix(X)
pred.res <- pred.env(m, X[2, ])
pred.res

}

