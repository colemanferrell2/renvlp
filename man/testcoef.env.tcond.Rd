\name{testcoef.env.tcond}
\alias{testcoef.env.tcond}
\title{Hypothesis test of the coefficients of the response envelope model with t-distributed errors}
\description{
This function tests the null hypothesis L * beta * R = A versus the alternative hypothesis L * beta * R ~= A, where beta is estimated under the response envelope model with t-distributed errors.
}
\usage{
testcoef.env.tcond(m, L, R, A)
}
\arguments{
  \item{m}{A list containing estimators and other statistics inherited from env.tcond.}
  \item{L}{The matrix multiplied to beta on the left.  It is a d1 by r matrix, while d1 is less than or equal to r.}
  \item{R}{The matrix multiplied to beta on the right.  It is a p by d2 matrix, while d2 is less than or equal to p.}
  \item{A}{The matrix on the right hand side of the equation.  It is a d1 by d2 matrix.}

Note that inputs \code{L}, \code{R} and \code{A} must be matrices, if not, use \code{as.matrix} to convert them.  
  }
\details{
This function tests for hypothesis H0: L beta R = A, versus Ha: L beta R != A.  The beta is estimated by the envelope model with t-distributed errors.  If L = Ir, R = Ip and A = 0, then the test is equivalent to the standard F test on if beta = 0.  The test statistic used is vec(L beta R - A) hat{Sigma}^{-1} vec(L beta R - A)^{T}, where beta is the envelope estimator and hat{Sigma} is the estimated asymptotic covariance of vec(L beta R - A). The reference distribution is chi-squared distribution with degrees of freedom d1 * d2. 
}
\value{The output is a list that contains following components.
\item{chisqStatistic}{The test statistic.} 
\item{dof}{The degrees of freedom of the reference chi-squared distribution.} 
\item{pValue}{p-value of the test.} 
\item{covMatrix}{The covariance matrix of vec(L beta R).} 
}

\examples{
data(concrete)
X <- concrete[1:78, 1:7]  # The first 78 observations are training data
Y <- concrete[1:78, 8:10]
m <- env.tcond(X, Y, 2, 6)
m

L <- diag(3)
R <- matrix(1, 7, 1)
A <- matrix(0, 3, 1)

test.res <- testcoef.env.tcond(m, L, R, A)
test.res
}

