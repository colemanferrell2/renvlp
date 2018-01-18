\name{u.pois.env}
\alias{u.pois.env}
\title{Select the dimension of pois.env}
\description{
This function outputs dimensions selected by Akaike information criterion (AIC), Bayesian information criterion (BIC) and likelihood ratio testing with specified significance level for the envelope model in poisson regression.
}
\usage{
u.pois.env(X, Y, alpha = 0.01)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors and n is number of observations.  The predictors must be continuous variables.}
  \item{Y}{Response. An n by 1 matrix. The univariate response must be counts.}
  \item{alpha}{Significance level for testing.  The default is 0.01.}
}
\value{
\item{u.aic}{Dimension of the envelope subspace selected by AIC.}
\item{u.bic}{Dimension of the envelope subspace selected by BIC.}
\item{u.lrt}{Dimension of the envelope subspace selected by the likelihood ratio testing procedure.}
\item{loglik.seq}{Log likelihood for dimension from 0 to p.}
\item{aic.seq}{AIC value for dimension from 0 to p.}
\item{bic.seq}{BIC value for dimension from 0 to p.}
}

\examples{
data(horseshoecrab)
X1 <- as.numeric(horseshoecrab[ , 1] == 2)
X2 <- as.numeric(horseshoecrab[ , 1] == 3)
X3 <- as.numeric(horseshoecrab[ , 1] == 4)
X4 <- as.numeric(horseshoecrab[ , 2] == 2)
X5 <- as.numeric(horseshoecrab[ , 2] == 3)
X6 <- horseshoecrab[ , 3]
X7 <- horseshoecrab[ , 5]
X <- cbind(X1, X2, X3, X4, X5, X6, X7)
Y <- horseshoecrab[ , 4]

u <- u.pois.env(X, Y)
u
}