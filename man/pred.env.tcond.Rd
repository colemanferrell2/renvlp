\name{pred.env.tcond}
\alias{pred.env.tcond}
\title{Estimation or prediction for env.tcond}
\description{
 Perform estimation or prediction under the response envelope model, where the errors follow a multivariate t-distribution.
}
\usage{
pred.env.tcond(m, Xnew)
}
\arguments{
  \item{m}{A list containing estimators and other statistics inherited from env.tcond.}
  \item{Xnew}{The value of X with which to estimate or predict Y.  A p dimensional vector.}
}
\details{
This function evaluates the envelope model with t distributed errors at new value Xnew.  It can perform estimation: find the fitted value when X = Xnew, or prediction: predict Y when X = Xnew.  The covariance matrix and the standard errors are also provided.
}
\value{The output is a list that contains following components.
\item{value}{The fitted value or the predicted value evaluated at Xnew.} 
\item{covMatrix.estm}{The covariance matrix of the fitted value at Xnew.} 
\item{SE.estm}{The standard error of the fitted value at Xnew.} 
\item{covMatrix.pred}{The covariance matrix of the predicted value at Xnew.} 
\item{SE.pred}{The standard error of the predicted value at Xnew.}
}

\examples{
data(concrete)
X <- concrete[1:78, 1:7]  # The first 78 observations are training data
Y <- concrete[1:78, 8:10]
\dontrun{u <- u.env.tcond(X, Y, 6)}
\dontrun{u}

m <- env.tcond(X, Y, 2, 6)
m

X <- as.matrix(X)
pred.res <- pred.env.tcond(m, X[2, ])
pred.res

}

