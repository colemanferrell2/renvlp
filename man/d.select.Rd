\name{d.select}
\alias{d.select}
\title{Select the rank of beta}
\description{
 This function outputs the rank selected by a chi-squared test developed by Bura and Cook (2003) with specified significance level for the beta. 
}
\usage{
d.select(X, Y, alpha = 0.01)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
  \item{alpha}{Significance level for testing. The default is 0.01.}
}
\details{
 This function estimate the rank of beta using a chi-squared test.  The test statistic and degrees of freedom are described in Bura and Cook (2003).
}
\value{
\item{rank.beta}{Rank of beta guided by the Bura-Cook estimator.}
}
\references{
Bura, E. and Cook, R. D. (2003). Rank estimation in reduced-rank regression. Journal of Multivariate Analysis, 87, 159 - 176.
}


