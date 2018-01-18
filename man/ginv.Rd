\name{ginv}
\alias{ginv}
\title{Generalized Inverse of a Matrix}
\description{
 Calculates the Moore-Penrose generalized inverse of a matrix X. 
}
\usage{
 ginv(X, tol = sqrt(.Machine$double.eps))
}
\arguments{
\item{X}{Matrix for which the Moore-Penrose inverse is required}
\item{tol}{A relative tolerance to detect zero sigular values}
}
\references{
 Venables, W.N. and Ripley, B.D. (1999) Modern Applied statistics with S-PLUS. Third Edition. Springer. p.100
}
