\name{env.tcond}
\alias{env.tcond}
\title{Fit the envelope model with t-distributed errors}
\description{
 For fixed envelope dimension u, fit the envelope model in multivariate linear regression, where the errors follow a multivariate t-distribution. 
}
\usage{
env.tcond(X, Y, u, df, asy = TRUE)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
  \item{u}{Dimension of the envelope. An integer between 0 and r.}
  \item{df}{Degrees of freedom of the t-distribution. A positive number that is greater than 2.}
  \item{asy}{Flag for computing the asymptotic variance of the envelope estimator.  The default is \code{TRUE}. When p and r are large, computing the asymptotic variance can take much time and memory.  If only the envelope estimators are needed, the flag can be set to \code{asy = FALSE}.}
}
\details{
 This function fits the envelope model to the responses and predictors, \deqn{
 Y_{i} = \alpha + \Gamma\eta X_{i}+\varepsilon_{i}, \Sigma=\Gamma\Omega\Gamma'+\Gamma_{0}\Omega_{0}\Gamma'_{0}, i=1, ..., n,
 }
 using the maximum likelihood estimation.  The errors \eqn{\varepsilon_{i}} follow a multivariate t-distribution with scale matrix \eqn{\Sigma} and degrees of freedom \code{df}.   When the dimension of the envelope is between 1 and r-1, the alternating algorithm in Forzani and Su (2021) is implemented.  When the dimension is r, then the envelope model degenerates to the standard multivariate linear regression.  When the dimension is 0, it means that X and Y are uncorrelated, and the fitting is different.
}
\value{The output is a list that contains the following components:
\item{beta}{The envelope estimator of the regression coefficients.} 
\item{Sigma}{The envelope estimator of the error covariance matrix.}
\item{Gamma}{An orthogonal basis of the envelope subspace.}
\item{Gamma0}{An orthogonal basis of the complement of the envelope subspace.}
\item{eta}{The coordinates of beta with respect to Gamma.}
\item{Omega}{The coordinates of Sigma with respect to Gamma.}
\item{Omega0}{The coordinates of Sigma with respect to Gamma0.}
\item{mu}{The estimated intercept.}
\item{loglik}{The maximized log likelihood function.}
\item{covMatrix}{The asymptotic covariance of vec(beta). The covariance matrix returned are asymptotic.  For the actual standard errors, multiply by 1 / n.}
\item{asySE}{The asymptotic standard error for elements in beta under the envelope model.  The standard errors returned are asymptotic, for actual standard errors, multiply by 1 / sqrt(n).}
\item{ratio}{The asymptotic standard error ratio of the standard multivariate linear regression estimator (with consideration of nonconstant variance) over the envelope estimator, for each element in beta.}
\item{n}{The number of observations in the data.}
}
\references{
Forzani, L. and Su, Z. (2021). Envelopes for elliptical multivariate linear regression. Statist. Sinica 31, 301-332.
}

\examples{
data(concrete)
X <- concrete[1:78, 1:7]  # The first 78 observations are training data
Y <- concrete[1:78, 8:10]
\dontrun{u <- u.env.tcond(X, Y, 6)}
\dontrun{u}

m <- env.tcond(X, Y, 2, 6)
m
m$beta
}

