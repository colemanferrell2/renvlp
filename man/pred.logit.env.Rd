\name{pred.logit.env}
\alias{pred.logit.env}
\title{Estimation or prediction for logit.env}
\description{
 Perform estimation or prediction under the envelope model in logistic regression.
}
\usage{
pred.logit.env(m, Xnew)
}
\arguments{
  \item{m}{A list containing estimators and other statistics inherited from xenv.}
  \item{Xnew}{The value of X with which to estimate or predict Y.  A p dimensional vector.}
}
\details{
This function evaluates the envelope model at new value Xnew.  It can perform estimation: find the fitted value when X = Xnew.  The covariance matrix of estimation and the standard errors of estimation are also provided.
}
\value{The output is a list that contains following components.
\item{value}{The fitted value or predicted value evaluated at Xnew.} 
\item{covMatrix.estm}{The covariance matrix of the fitted value at Xnew.} 
\item{SE.estm}{The standard error of the fitted value at Xnew.}  
}

\examples{
data(horseshoecrab)	
X1 <- as.numeric(horseshoecrab[ , 1] == 2)
X2 <- as.numeric(horseshoecrab[ , 1] == 3)
X3 <- as.numeric(horseshoecrab[ , 1] == 4)
X4 <- as.numeric(horseshoecrab[ , 2] == 2)
X5 <- as.numeric(horseshoecrab[ , 2] == 3)
X6 <- horseshoecrab[ , 3]
X7 <- horseshoecrab[ , 5]
X <- cbind(X1, X2, X3, X4, X5, X6, X7)
Y <- as.numeric(ifelse(horseshoecrab[ , 4] > 0, 1, 0))

m <- logit.env(X, Y, 1)

pred.res <- pred.logit.env(m, X[1, ])
pred.res

}

