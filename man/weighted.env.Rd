\name{weighted.env}
\alias{weighted.env}
\title{Weighted envelope estimator}
\description{
 Compute the weighted envelope estimator with weights computed from BIC. 
}
\usage{
weighted.env(X, Y, bstrpNum = 0, min.u = 1, 
max.u = ncol(as.matrix(Y)), boot.resi = "full")
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
  \item{bstrpNum}{Number of bootstrap samples. A positive integer.}
  \item{min.u}{Lower bound of the range of u to compute bootstrap error. A postive integer between 1 and p. This argument is relevant only when \code{bstrpNum}>0.}
  \item{max.u}{Upper bound of the range of u to compute bootstrap error. A postive integer between 1 and p. This argument is relevant only when \code{bstrpNum}>0.}
  \item{boot.resi}{A string that can be "full" or "weighted" indicating the model from which the residuals are calculated.  If the input is "full", then the residuals are obtained using the standard estimators; and if the input is "weighted", then the residuals are obtained using the weighted envelope estimators.  This argument is for computing residuals in residual bootstrap, and it is relevant only when \code{bstrpNum}>0.}
}
\details{
 This function computes the weighted envelope estimator in a standard multivariate linear regression.  And the weighted envelope estimator takes the form \deqn{
 \hat{\beta}_{w}=\sum_{j=1}^{r}w_{j}\hat{\beta}_{j},
 }
 where \eqn{\hat{\beta}_{j}} is the envelope estimator of \eqn{\beta} with \eqn{u=j} and \eqn{w_{j}}'s are the weights computed from BIC values
\deqn{
w_{j}=\frac{\exp(-b_{j})}{\sum_{k=1}^{r}\exp(-b_{k})},
}
where \eqn{b_{j}} is the BIC criterion evaluated at the envelope estimator \eqn{\hat{\beta}_{j}}.  For details, see Eck and Cook (2017).

The variation of the weighted envelope estimator is estimated by residual bootstrap.  The user can specify the range for bootstrap \eqn{u=}(\code{min.u}, \code{max.u}), if the weights outside of the range are small.
}
\value{The output is a list that contains the following components:
\item{beta}{The weighted envelope estimator of the regression coefficients.} 
\item{mu}{The weighted estimated intercept.}
\item{Sigma}{The weighted envelope estimator of the error covariance matrix.}
\item{w}{Weights computed based on BIC.}
\item{loglik}{The log likelihood function computed with weighted envelope estimator.}
\item{n}{The number of observations in the data.}
\item{bootse}{The standard error for elements in beta computed by residual bootstrap.  This output is available only when \code{bstrpNum}>0.}
\item{ratios}{The boostrap standard error ratio of the standard multivariate linear regression estimator over the weighted envelope estimator for each element in beta.  This output is available only when \code{bstrpNum}>0.}
\item{bic_select}{A table that lists how many times BIC selected each candidate dimension.  If BIC never selects a dimension, this dimension does not appear on the table.  This output is available only when \code{bstrpNum}>0.}
}
\references{
Eck, D. J. and Cook, R. D. (2017). Weighted Envelope Estimation to Handle Variability in Model
Selection. Biometrika. To appear.
}

\examples{
data(wheatprotein)
X <- wheatprotein[, 8]
Y <- wheatprotein[, 1:6]
m <- weighted.env(X, Y)
m$w
m$beta

\dontrun{m2 <- weighted.env(X, Y, bstrpNum = 100, min.u = 1, max.u = 6, boot.resi = "full")}
\dontrun{m2$bic_select}
\dontrun{m2$bootse}
}

