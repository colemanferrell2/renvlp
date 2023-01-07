\name{rrenv.apweights}
\alias{rrenv.apweights}
\title{Fit the reduced-rank envelope model with nonconstant variance}
\description{
 For rank d and fixed envelope dimension u, fit the reduced-rank envelope model with nonconstant error variance.  
}
\usage{
rrenv.apweights(X, Y, u, d, asy = TRUE)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
  \item{u}{Dimension of the envelope. An integer between 0 and r.}
  \item{d}{The rank of the coefficient matrix.  An integer between 0 and u.}
  \item{asy}{Flag for computing the asymptotic variance of the reduced rank envelope estimator.  The default is \code{TRUE}. When p and r are large, computing the asymptotic variance can take much time and memory.  If only the envelope estimators are needed, the flag can be set to \code{asy = FALSE}.}
}
\details{
 This function fits the reduced rank envelope model to the responses and predictors, \deqn{
 Y_{i} = \alpha + \Gamma\eta BX_{i}+\varepsilon_{i}, \Sigma=c_{i}(\Gamma\Omega\Gamma'+\Gamma_{0}\Omega_{0}\Gamma'_{0}), i=1, ..., n,
 }
 using the maximum likelihood estimation.  The errors \eqn{\varepsilon_{i}} follow a normal distribution.  It allows that the error covariance matrix to be nonconstant.  When 0 < d < u < r, the estimation procedure in Cook et al. (2015) is implemented.  When d < u = r, then the model is equivalent to a reduced rank regression model.  When d = u, or d = p < r, then B can be taken as the identity matrix and the model reduces to a response envelope model.  When the dimension is d = u = r, then the envelope model degenerates to the standard multivariate linear regression.  When the u = 0, it means that X and Y are uncorrelated, and the fitting is different.
}
\value{The output is a list that contains the following components:
\item{Gamma}{An orthogonal basis of the envelope subspace.}
\item{Gamma0}{An orthogonal basis of the complement of the envelope subspace.}
\item{mu}{The estimated intercept.}
\item{beta}{The envelope estimator of the regression coefficients.} 
\item{Sigma}{The envelope estimator of the error covariance matrix.}
\item{eta}{The eta matrix in the coefficient matrix.}
\item{B}{The B matrix in the coefficient matrix.}
\item{Omega}{The coordinates of Sigma with respect to Gamma.}
\item{Omega0}{The coordinates of Sigma with respect to Gamma0.}
\item{loglik}{The maximized log likelihood function.}
\item{covMatrix}{The asymptotic covariance of vec(beta). The covariance matrix returned are asymptotic.  For the actual standard errors, multiply by 1 / n.}
\item{asySE}{The asymptotic standard error for elements in beta under the reduced rank envelope model.  The standard errors returned are asymptotic, for actual standard errors, multiply by 1 / sqrt(n).}
\item{ratio}{The asymptotic standard error ratio of the standard multivariate linear regression estimator (with consideration of nonconstant variance) over the envelope estimator, for each element in beta.}
\item{n}{The number of observations in the data.}
\item{C1}{The estimated weights \eqn{c_i}.}
}
\references{
Cook, R. D., Forzani, L. and Zhang, X. (2015). Envelopes and reduced-rank regression.  Biometrika 102, 439-456.

Forzani, L. and Su, Z. (2021). Envelopes for elliptical multivariate linear regression. Statist. Sinica 31, 301-332.
}

\examples{
data(vehicles)
X <- vehicles[, 1:11]  
Y <- vehicles[, 12:15]
X <- scale(X)
Y <- scale(Y)  # The scales of Y are vastly different, so scaling is reasonable here
d <- d.select(X, Y, 0.01)
d

\dontrun{u <- u.rrenv.apweights(X, Y, 2)}
\dontrun{u}

\dontrun{m <- rrenv.apweights(X, Y, 3, 2)}
\dontrun{m}
\dontrun{m$beta}
}

