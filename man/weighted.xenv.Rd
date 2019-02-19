\name{weighted.xenv}
\alias{weighted.xenv}
\title{Weighted predictor envelope estimator}
\description{
 Compute the weighted predictor envelope estimator with weights computed from BIC. 
}
\usage{
weighted.xenv(X, Y, bstrpNum = 0, min.u = 1, 
max.u = ncol(as.matrix(X)), boot.resi = "full")
}
\arguments{
\item{X}{Predictors. An n by p matrix, p is the number of predictors and n is number of observations.  The predictors must be continuous variables.}
  \item{Y}{Responses. An n by r matrix, r is the number of responses. The response can be univariate or multivariate and must be continuous variable.}
  \item{bstrpNum}{Number of bootstrap samples. A positive integer.}
  \item{min.u}{Lower bound of the range of u to compute bootstrap error. A postive integer between 1 and p. This argument is relevant only when \code{bstrpNum}>0.}
  \item{max.u}{Upper bound of the range of u to compute bootstrap error. A postive integer between 1 and p. This argument is relevant only when \code{bstrpNum}>0.}
  \item{boot.resi}{A string that can be "full" or "weighted" indicating the model from which the residuals are calculated.  If the input is "full", then the residuals are obtained using the standard estimators; and if the input is "weighted", then the residuals are obtained using the weighted predictor envelope estimators.  This argument is for computing residuals in residual bootstrap, and it is relevant only when \code{bstrpNum}>0.}
}
\details{
 This function computes the weighted predictor envelope estimator in a standard multivariate linear regression.  And the weighted predictor envelope estimator takes the form \deqn{
 \hat{\beta}_{w}=\sum_{j=1}^{p}w_{j}\hat{\beta}_{j},
 }
 where \eqn{\hat{\beta}_{j}} is the predictor envelope estimator of \eqn{\beta} with \eqn{u=j} and \eqn{w_{j}}'s are the weights computed from BIC values
\deqn{
w_{j}=\frac{\exp(-b_{j})}{\sum_{k=1}^{p}\exp(-b_{k})},
}
where \eqn{b_{j}} is the BIC criterion evaluated at the predictor envelope estimator \eqn{\hat{\beta}_{j}}.  For details, see Eck and Cook (2017).

The variation of the weighted predictor envelope estimator is estimated by residual bootstrap.  The user can specify the range for bootstrap \eqn{u=}(\code{min.u}, \code{max.u}), if the weights outside of the range are small.
}
\value{The output is a list that contains the following components:
\item{beta}{The weighted predictor envelope estimator of the regression coefficients.} 
\item{mu}{The weighted estimated intercept.}
\item{SigmaX}{The weighted predictor envelope estimator of the covariance matrix of X.}
\item{SigmaYcX}{The weighted predictor envelope estimator of the error covariance matrix.}
\item{w}{Weights computed based on BIC.}
\item{loglik}{The log likelihood function computed with weighted predictor envelope estimator.}
\item{n}{The number of observations in the data.}
\item{bootse}{The standard error for elements in beta computed by residual bootstrap.  This output is available only when \code{bstrpNum}>0.}
\item{ratios}{The boostrap standard error ratio of the standard multivariate linear regression estimator over the weighted predictor envelope estimator for each element in beta.  This output is available only when \code{bstrpNum}>0.}
\item{bic_select}{A table that lists how many times BIC selected each candidate dimension.  If BIC never selects a dimension, this dimension does not appear on the table.  This output is available only when \code{bstrpNum}>0.}
}
\references{
Eck, D. J. and Cook, R. D. (2017). Weighted Envelope Estimation to Handle Variability in Model
Selection. Biometrika. To appear.
}

\examples{
data(wheatprotein)
X <- wheatprotein[, 1:6]
Y <- wheatprotein[, 7]
m <- weighted.xenv(X, Y)
m$w
m$beta

\dontrun{m2 <- weighted.xenv(X, Y, bstrpNum = 100, min.u = 2, max.u = 4, boot.resi = "full")}
\dontrun{m2$w}
\dontrun{m2$bootse}
}

