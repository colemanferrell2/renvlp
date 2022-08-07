\name{u.felmdir}
\alias{u.felmdir}
\title{Find the envelope dimensions in the functional envelope linear model}
\description{
  Fit the dimensions of the response and predictor envelopes in function-on-function linear regression, under direct estimation.
}
\usage{
  u.felmdir(X, Y, t1, t2, knots = c(0, 0.25, 0.5, 0.75, 1))
}
\arguments{
  \item{X}{Predictor function. An n by T1 matrix, T1 is number of observed time points, which is the length of t1.  Here we assume that each function is observed at the same time points.}
  \item{Y}{Response function. An n by T2 matrix, T2 is number of observed time points, which is the length of t2.  Here we assume that each function is observed at the same time points.}
  \item{t1}{The observed time points for the predictor functions.}
  \item{t2}{The observed time points for the response functions.}
  \item{knots}{The location of knots of the cubic splines used for estimation.  Locations should be positive. The default location of the knots are 0, 0.25, 0.5, 0.75, 1.}
}
\details{
  This function finds the dimension of the predictor and response envelope model by Bayesian information criterion (BIC) performed on the direct estimation.  To be more specific, consider the envelope model to the function-on-function linear regression, \deqn{
 Y = \alpha + B X + \epsilon, 
 }
where X and Y are random functions in Hilbert spaces \eqn{H_X} and \eqn{H_Y}, \eqn{\alpha} is a fixed member in \eqn{H_Y}, \eqn{\epsilon} is a random member of \eqn{H_Y}, and B: \eqn{H_X -> H_Y} is a linear operator.  We use cubic splines as the basis for both \eqn{H_X} and \eqn{H_Y}. The coefficients \eqn{[X]} and \eqn{[Y]} with respect to the basis are computed.  The predictor and response envelope model is fitted on the linear regression model of \eqn{[Y]} on \eqn{[X]}, and the dimensions of the predictor and response envelopes are calculated using BIC. The details are included in Section 7 of Su et al. (2022).  
}
\value{The output is a list that contains the following components:
  \item{ux}{The estimated dimension of the predictor envelope.}
  \item{uy}{The estimated dimension of the response envelope.}
  \item{beta}{The envelope estimator of the regression coefficients in the regression of \eqn{[Y]} on \eqn{[X]}, when the dimensions of envelopes are taken at their estimated values.} 
  \item{betafull}{The standard estimator, i.e., the OLS estimator of the regression coefficients in the regression of \eqn{[Y]} on \eqn{[X]}.}
  \item{alpha}{The envelope estimator of the intercept in the regression of \eqn{[Y]} on \eqn{[X]}, when the dimensions of envelopes are taken at their estimated values.}
  \item{alphafull}{The standard estimator of the intercept in the regression of \eqn{[Y]} on \eqn{[X]}.}
}
\references{
  Su, Z., Li, B. and Cook, R. D. (2022+) Envelope model for function-on-function linear regression.
}

\examples{
  data(NJdata)
  dataX <- matrix(NJdata[,6], nrow = 21)
  X <- as.matrix(dataX[, 32:61])
  dataY <- matrix(NJdata[,3], nrow = 21)
  Y <- as.matrix(dataY[, 32:61])
  t1 <- 0:29
  t2 <- t1
  
  \dontrun{m <- u.felmdir(X, Y, t1, t2)}
  \dontrun{m$ux}
  \dontrun{m$uy}
}

