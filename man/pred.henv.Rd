\name{pred.henv}
\alias{pred.henv}
\title{Estimation or prediction for henv}
\description{
 Perform estimation or prediction under the heteroscedastic envelope model.
}
\usage{
pred.henv(m, Xnew)
}
\arguments{
  \item{m}{A list containing estimators and other statistics inherited from henv.}
  \item{Xnew}{The value of X with which to estimate or predict Y. An r by 1 vector.}
}
\details{
This function evaluates the heteroscedastic envelope model at new value Xnew.  It can perform estimation: find the fitted value when X = Xnew, or prediction: predict Y when X = Xnew.  The covariance matrix and the standard errors are also provided.
}
\value{The output is a list that contains following components.
\item{value}{The fitted value or the predicted value evaluated at Xnew.} 
\item{covMatrix.estm}{The covariance matrix of the fitted value at Xnew.} 
\item{SE.estm}{The standard error of the fitted value at Xnew.} 
\item{covMatrix.pred}{The covariance matrix of the predicted value at Xnew.} 
\item{SE.pred}{The standard error of the predicted value at Xnew.} 
}

\examples{
data(waterstrider)
X <- waterstrider[ , 1]
Y <- waterstrider[ , 2:5]


\dontrun{m <- henv(X, Y, 2)}

\dontrun{pred.res <- pred.henv(m, X[2])}

}