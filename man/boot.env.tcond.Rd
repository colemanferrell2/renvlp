\name{boot.env.tcond}
\alias{boot.env.tcond}
\title{Bootstrap for env.tcond}
\description{
 Compute bootstrap standard error for the response envelope estimator with t-distributed errors.
}
\usage{
boot.env.tcond(X, Y, u, df, B)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
  \item{u}{Dimension of the envelope. An integer between 0 and r.}
  \item{df}{Degrees of freedom of the t-distribution. A positive number that is greater than 2.}
  \item{B}{Number of bootstrap samples. A positive integer.}
}
\details{
 This function computes the bootstrap standard errors for the regression coefficients in the envelope model with t-distributed errors by bootstrapping the residuals.
}
\value{The output is an r by p matrix.
\item{bootse}{The standard error for elements in beta computed by bootstrap.} 
}

\examples{
data(concrete)
X <- concrete[1:78, 1:7]  # The first 78 observations are training data
Y <- concrete[1:78, 8:10]

\dontrun{u <- u.env.tcond(X, Y, 6)}
\dontrun{u}

B <- 100
\dontrun{bootse <- boot.env.tcond(X, Y, 2, 6, B)}
\dontrun{bootse}
}

