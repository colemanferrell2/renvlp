\name{pred.eppls}
\alias{pred.eppls}
\title{Estimation or prediction for eppls}
\description{
 Perform estimation or prediction under the Envelope-based Partial Partial Least Squares.
}
\usage{
pred.eppls(m, X1new, X2new)
}
\arguments{
  \item{m}{A list containing estimators and other statistics inherited from eppls.}
  \item{X1new}{The value of X1 with which to estimate or predict Y.  A p1 dimensional vector.}
  \item{X2new}{The value of X2 with which to estimate or predict Y.  A p2 dimensional vector.}
}
\details{
This function evaluates the partial envelope model at new value Xnew.  It can perform estimation: find the fitted value when X = Xnew, or prediction: predict Y when X = Xnew.  The covariance matrix and the standard errors are also provided.
}
\value{The output is a list that contains following components.

  \item{value}{
  The fitted value or the predicted value evaluated at X1new and X2new.
  } 
  \item{covMatrix.estm}{
  The covariance matrix of the fitted value at X1new and X2new.
  } 
  \item{SE.estm}{
  The standard error of the fitted value at X1new and X2new.
  } 
  \item{covMatrix.pred}{
  The covariance matrix of the predicted value at X1new and X2new.
  } 
  \item{SE.pred}{
  The standard error of the predicted value at X1new and X2new.
  }  
}
\references{
Park, Y., Su, Z. and Chung, D. (2022+)  Envelope-based Partial Partial Least Squares with Application to Cytokine-based Biomarker Analysis for COVID-19.
}

\examples{
  data(amitriptyline)
  
  Y <- amitriptyline[ , 1:2]
  X1 <- amitriptyline[ , 4:7]
  X2 <- amitriptyline[ , 3]
  u <- u.eppls(X1, X2, Y)
  
  u
  
  m <- eppls(X1, X2, Y, 2)
  
  pred.res <- pred.eppls(m, X1[1, ], X2[1])
  pred.res
}
