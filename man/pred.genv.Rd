\name{pred.genv}
\alias{pred.genv}
\title{Estimation or prediction for genv}
\description{
 Perform estimation or prediction under the groupwise envelope model.
}
\usage{
pred.genv(m, Xnew, Znew)
}
\arguments{
  \item{m}{A list containing estimators and other statistics inherited from env.}
  \item{Xnew}{The value of X with which to estimate or predict Y.  A p dimensional vector.}
  \item{Znew}{A group indicator of X.}
}
\details{
This function evaluates the envelope model at new value Xnew.  It can perform estimation: find the fitted value when X = Xnew with a group indicator Z = Znew, or prediction: predict Y when X = Xnew with a group indicator Z = Znew.  The covariance matrix and the standard errors are also provided.
}
\value{The output is a list that contains following components.
\item{value}{The fitted value or the predicted value evaluated at Znew.} 
\item{covMatrix.estm}{The covariance matrix of the fitted value at Znew.} 
\item{SE.estm}{The standard error of the fitted value at Znew.} 
\item{covMatrix.pred}{The covariance matrix of the predicted value at Znew.} 
\item{SE.pred}{The standard error of the predicted value at Znew.}
}

\examples{
data(fiberpaper)
X <- fiberpaper[ , c(5, 7)]
Y <- fiberpaper[ , 1:3]
Z <- as.numeric(fiberpaper[ , 6] > mean(fiberpaper[ , 6]))
u <- u.genv(X, Y, Z)
u

m <- genv(X, Y, Z, 2)
m

X <- as.matrix(X)
pred.res <- pred.genv(m, X[2, ], Z[2])
pred.res

}