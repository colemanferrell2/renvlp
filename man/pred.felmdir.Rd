\name{pred.felmdir}
\alias{pred.felmdir}
\title{Estimation or prediction for felmdir}
\description{
 Perform estimation or prediction under the functional envelope linear model, using the direct estimation.
}
\usage{
pred.felmdir(X, Y, ux, uy, t1, t2, Xnew, knots = c(0, 0.25, 0.5, 0.75, 1))
}
\arguments{
  \item{X}{Predictor function. An n by T1 matrix, T1 is number of observed time points, which is the length of t1.  Here we assume that each function is observed at the same time points.}
  \item{Y}{Response function. An n by T2 matrix, T2 is number of observed time points, which is the length of t2.  Here we assume that each function is observed at the same time points.}
  \item{ux}{Dimension of the predictor envelope. An integer between 0 and number of knots +2.}
  \item{uy}{Dimension of the response envelope. An integer between 0 and number of knots +2.}
  \item{t1}{The observed time points for the predictor functions.}
  \item{t2}{The observed time points for the response functions.}
  \item{Xnew}{The value of X with which to estimate or predict Y.  A T1 dimensional vector. The observed time points should be the same as those of X.}
  \item{knots}{The location of knots of the cubic splines used for estimation.  Locations should be positive. The default location of the knots are 0, 0.25, 0.5, 0.75, 1.}
  }
\details{
This function evaluates the functional envelope linear model at new value Xnew.  It can perform estimation: find the fitted value when X = Xnew, or prediction: predict Y when X = Xnew.  The covariance matrix and the standard errors are also provided.  The estimation method uses the direct estimation in Su et al. (2022) with cubic splines.
}
\value{The output is a list that contains following components.
\item{value}{The fitted value or the predicted value evaluated at Xnew.  The fitted or predicted values are at the same observation points as Y.} 
\item{covMatrix.estm}{The covariance matrix of the fitted value at Xnew.} 
\item{SE.estm}{The standard error of the fitted value at Xnew.} 
\item{covMatrix.pred}{The covariance matrix of the predicted value at Xnew.} 
\item{SE.pred}{The standard error of the predicted value at Xnew.}
}

\references{
Su, Z., Li, B. and Cook, R. D. (2022+) Envelope model for function-on-function linear regression.
}

\examples{
data(NJdata)
dataX <- matrix(NJdata[,6], nrow = 21)
X <- as.matrix(dataX[, 32:61])
dataY <- matrix(NJdata[,3], nrow = 21)
Y <- as.matrix(dataY[, 32:61])
t1 <- 0:29
t2 <- t1

m <- pred.felmdir(X, Y, 3, 1, t1, t2, X[1,])
m$value
m$SE.estm
m$SE.pred
}

