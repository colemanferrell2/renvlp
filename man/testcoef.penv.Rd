\name{testcoef.penv}
\alias{testcoef.penv}
\title{Hypothesis test of the coefficients of the partial envelope model}
\description{
This function tests the null hypothesis L * beta1 * R = A versus the alternative hypothesis L * beta1 * R ~= A, where beta is estimated under the partial envelope model.
}
\usage{
testcoef.penv(m, L, R, A)
}
\arguments{
  \item{m}{A list containing estimators and other statistics inherited from penv.}
  \item{L}{The matrix multiplied to beta on the left.  It is a d1 by r matrix, while d1 is less than or equal to r.}
  \item{R}{The matrix multiplied to beta on the right.  It is a p1 by d2 matrix, while d2 is less than or equal to p1.}
  \item{A}{The matrix on the right hand side of the equation.  It is a d1 by d2 matrix.}
  
Note that inputs \code{L}, \code{R} and \code{A} must be matrices, if not, use \code{as.matrix} to convert them.  
  }
\details{
This function tests for hypothesis H0: L beta1 R = A, versus Ha: L beta1 R != A.  The beta is estimated by the partial envelope model.  If L = Ir, R = Ip1 and A = 0, then the test is equivalent to the standard F test on if beta1 = 0.  The test statistics used is vec(L beta1 R - A) hat{Sigma}^{-1} vec(L beta1 R - A)^{T}, where beta is the envelope estimator and hat{Sigma} is the estimated asymptotic covariance of vec(L beta1 R - A). The reference distribution is chi-squared distribution with degrees of freedom d1 * d2. 
}
\value{The output is a list that contains following components.
\item{chisqStatistic}{The test statistic.} 
\item{dof}{The degrees of freedom of the reference chi-squared distribution.} 
\item{pValue}{p-value of the test.} 
\item{covMatrix}{The covariance matrix of vec(L beta1 R).} 
}

\examples{
data(fiberpaper)
X1 <- fiberpaper[, 7]
X2 <- fiberpaper[, 5:6]
Y <- fiberpaper[, 1:4]
m <- penv(X1, X2, Y, 1)
m

L <- diag(4)
R <- as.matrix(1)
A <- matrix(0, 4, 1)

test.res <- testcoef.penv(m, L, R, A)
test.res
}

