\name{pred2.env}
\alias{pred2.env}
\title{Estimation or prediction for env}
\description{
 Perform estimation or prediction under the envelope model through partial envelope model.
}
\usage{
pred2.env(X, Y, u, Xnew)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
  \item{u}{The dimension of the constructed partial envelope model.}
  \item{Xnew}{The value of X with which to estimate or predict Y.  A p dimensional vector.}
}
\details{
This function evaluates the envelope model at new value Xnew.  It can perform estimation: find the fitted value when X = Xnew, or prediction: predict Y when X = Xnew.  The covariance matrix and the standard errors are also provided.  Compared to predict.env, this function performs prediction through partial envelope model, which can be more accurate if the partial envelope is of smaller dimension and contains less variant material information.  The constructed partial envelope model is obtained by the following: Let A0 by a p by p-1 matrix, such that A = (Xnew, A0) has full rank.  Let phi1 = beta * Xnew, phi2 = beta * A0, phi = (phi1, phi2) and X = inverse of A * X = (Z1, Z2')'.  Then the model Y = alpha + beta * X + epsilon can be reparameterized as Y = alpha + phi1 * Z1 + phi2 * Z2 + epsilon.  We then fit a partial envelope model with Z1 as the predictor of interest, and predict at (Z1, Z2')' =  inverse of A * Xnew.
}
\value{The output is a list that contains following components.
\item{value}{The fitted value or the predicted value evaluated at Xnew.} 
\item{covMatrix.estm}{The covariance matrix of the fitted value at Xnew.} 
\item{SE.estm}{The standard error of the fitted value at Xnew.} 
\item{covMatrix.pred}{The covariance matrix of the predicted value at Xnew.} 
\item{SE.pred}{The standard error of the predicted value at Xnew.}  
}

\examples{
data(fiberpaper)
X <- fiberpaper[, 5:7]
Y <- fiberpaper[, 1:4]

u <- u.pred2.env(X, Y, X[10, ])
pred.res <- pred2.env(X, Y, u$u.bic, X[10, ])
pred.res$SE.estm
pred.res$SE.pred
}

