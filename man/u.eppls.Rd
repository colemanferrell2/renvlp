\name{u.eppls}
\alias{u.eppls}
\title{Select the dimension of eppls}
\description{
 This function outputs dimensions selected by Bayesian information criterion (BIC) for the Envelope-based Partial Partial Least Squares.
}


\usage{
u.eppls(X1, X2, Y)
}

\arguments{
  \item{X1}{
An \eqn{n} by \eqn{p1} matrix of continuous predictors, where \eqn{p1} is the  number of continuous predictors with \eqn{p1 < n}.
}
  \item{X2}{
An \eqn{n} by \eqn{p2} matrix of categorical predictors, where \eqn{p2} is the  number of categorical predictors with \eqn{p2 < n}.
}
  \item{Y}{
An \eqn{n} by \eqn{r} matrix of multivariate responses, where \eqn{r} is the number of responses.
}
}

\value{
  \item{u.bic}{
Dimension of the Envelope-based Partial Partial Least Squares.
  }
  \item{bic.seq}{
BIC value for dimension from \eqn{0} to \eqn{p1}.
  }
}

\references{
Park, Y., Su, Z. and Chung, D. (2022+)  Envelope-based Partial Partial Least Squares with Application to Cytokine-based Biomarker Analysis for COVID-19.
}

\examples{
  data(amitriptyline)
  
  Y <- amitriptyline[ , 1:2]
  X1 <- amitriptyline[ , 4:7]
  X2 <- amitriptyline[ , 3]
  u <- u.eppls(X1, X2, Y)
  
  u
}

