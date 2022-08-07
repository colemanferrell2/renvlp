\name{felmKL}
\alias{felmKL}
\title{Fit the functional envelope linear model}
\description{
 Fit the response and predictor envelope model in function-on-function linear regression with dimensions ux and uy, using Karhunen-Loeve expansion based estimation. 
}
\usage{
felmKL(X, Y, ux, uy, t1, t2, knots = c(0, 0.25, 0.5, 0.75, 1))
}
\arguments{
  \item{X}{Predictor function. An n by T1 matrix, T1 is number of observed time points, which is the length of t1.  Here we assume that each function is observed at the same time points.}
  \item{Y}{Response function. An n by T2 matrix, T2 is number of observed time points, which is the length of t2.  Here we assume that each function is observed at the same time points.}
  \item{ux}{Dimension of the predictor envelope. An integer between 0 and number of knots +2.}
  \item{uy}{Dimension of the response envelope. An integer between 0 and number of knots +2.}
  \item{t1}{The observed time points for the predictor functions.}
  \item{t2}{The observed time points for the response functions.}
  \item{knots}{The location of knots of the cubic splines used for estimation.  Locations should be positive. The default location of the knots are 0, 0.25, 0.5, 0.75, 1.}
}
\details{
 This function fits the envelope model to the function-on-function linear regression, \deqn{
 Y = \alpha + B X + \epsilon, 
 }
where X and Y are random functions in Hilbert spaces \eqn{H_X} and \eqn{H_Y}, \eqn{\alpha} is a fixed member in \eqn{H_Y}, \eqn{\epsilon} is a random member of \eqn{H_Y}, and B: \eqn{H_X -> H_Y} is a linear operator.  We use cubic splines as the basis for both \eqn{H_X} and \eqn{H_Y} in the estimation of the eigenfunctions of \eqn{Sigma_X} and \eqn{Sigma_\epsilon}. The coefficients \eqn{[X]} and \eqn{[Y]} with respect to the estimated eigenfunctions are computed.  The predictor and response envelope model is fitted on the linear regression model of \eqn{[Y]} on \eqn{[X]}.  Based on its result, the fitted value of \eqn{Y} is calculated.  The standard function-on-function regression model also works through the linear regression model of \eqn{[Y]} on \eqn{[X]}. But instead of fitting an envelope model, it fits a standard linear regression model, based on which the fitted value of \eqn{Y} is calculated.  The details are elaborated in Section 6, Karhunen-Lo`eve expansion based estimation, in the reference of Su et al. (2022).  
}
\value{The output is a list that contains the following components:
\item{beta}{The envelope estimator of the regression coefficients in the regression of \eqn{[Y]} on \eqn{[X]}.} 
\item{betafull}{The standard estimator, i.e., the OLS estimator of the regression coefficients in the regression of \eqn{[Y]} on \eqn{[X]}.}
\item{alpha}{The envelope estimator of the intercept in the regression of \eqn{[Y]} on \eqn{[X]}.}
\item{alphafull}{The standard estimator of the intercept in the regression of \eqn{[Y]} on \eqn{[X]}.}
\item{phihat.cord}{The estimated coordinates of eigenfunctions of \eqn{Sigma_\epsilon} with respect to the cubic splines.}
\item{psihat.cord}{The estimated coordinates of eigenfunctions of \eqn{Sigma_X} with respect to the cubic splines.}
\item{fitted.env}{The fitted value of Y computed from the functional envelope linear model.}
\item{fitted.full}{The fitted value of Y computed from the standard function-to-function linear model.}
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

m <- felmKL(X, Y, 4, 3, t1, t2)
head(m$fitted.env)
head(m$fitted.full)
}

