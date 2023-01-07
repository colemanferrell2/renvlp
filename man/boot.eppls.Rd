\name{boot.eppls}
\alias{boot.eppls}
\title{Bootstrap for eppls}
\description{
 Compute bootstrap standard error for the Envelope-based Partial Partial Least Squares estimator.
}
\usage{
boot.eppls(X1, X2, Y, u, B)
}
\arguments{
   \item{X1}{
An \eqn{n} by \eqn{p1} matrix of continuous predictors, where \eqn{p1} is the  number of continuous predictors with \eqn{p1 < n}.
 }
  \item{X2}{
An \eqn{n} by \eqn{p2} matrix of categorical predictors, where \eqn{p2} is the  number of categorical predictors with \eqn{p2 < n}.
 }
  \item{Y}{
An \eqn{n} by \eqn{r} matrix of multivariate responses, where \eqn{r} is the number of responses.
 }
  \item{u}{
A given dimension of the Envelope-based Partial Partial Least Squares.  It should be an interger between \eqn{0} and \eqn{p1}.
 }

  \item{B}{
  Number of bootstrap samples. A positive integer.
  }
}
\details{
 This function computes the bootstrap standard errors for the regression coefficients beta1 and beta2 in the Envelope-based Partial Partial Least Squares by bootstrapping the residuals.
}
\value{The output is a list that contains the following components:
\item{bootse1}{
The standard error for elements in beta1 computed by bootstrap. The output is an p1 by r matrix.
} 
\item{bootse1}{
The standard error for elements in beta2 computed by bootstrap. The output is an p2 by r matrix.
} 
}

\examples{
data(amitriptyline)
  
Y <- amitriptyline[ , 1:2]
X1 <- amitriptyline[ , 4:7]
X2 <- amitriptyline[ , 3]

B <- 100
\dontrun{bootse <- boot.eppls(X1, X2, Y, 2, B)}
\dontrun{bootse}
}

