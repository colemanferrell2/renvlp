\name{testcoef.pois.env}
\alias{testcoef.pois.env}
\title{Hypothesis test of the coefficients of the envelope model}
\description{
This function tests the null hypothesis L * beta = A versus the alternative hypothesis L * beta ~= A, where beta is estimated under the envelope model in poisson regression.
}
\usage{
testcoef.pois.env(m, L, A)
}
\arguments{
  \item{m}{A list containing estimators and other statistics inherited from pois.env.}
    \item{L}{The matrix multiplied to beta on the left.  It is a d1 by p matrix, while d1 is less than or equal to p.}
  \item{A}{The matrix on the right hand side of the equation.  It is a d1 by 1 matrix.}

Note that inputs \code{L} and \code{A} must be matrices, if not, use \code{as.matrix} to convert them.   
  }
\details{
This function tests for hypothesis H0: L beta = A, versus Ha: L beta != A.  The beta is estimated by the envelope model in predictor space.  If L = Ip and A = 0, then the test is equivalent to the standard F test on if beta = 0.  The test statistic used is vec(L beta - A) hat{Sigma}^{-1} vec(L beta - A)^{T}, where beta is the envelope estimator and hat{Sigma} is the estimated asymptotic covariance of vec(L beta - A). The reference distribution is chi-squared distribution with degrees of freedom d1. 
}
\value{The output is a list that contains following components.
\item{chisqStatistic}{The test statistic.} 
\item{dof}{The degrees of freedom of the reference chi-squared distribution.} 
\item{pValue}{p-value of the test.} 
\item{covMatrix}{The covariance matrix of vec(L beta).} 
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
Y <- horseshoecrab[ , 4]

m <- pois.env(X, Y, 1)

L <- diag(7)
A <- matrix(0, 7, 1)

test.res <- testcoef.pois.env(m, L, A)
test.res
}
