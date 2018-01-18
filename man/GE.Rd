\name{GE}
\alias{GE}
\title{Gaussian elimination}
\description{
 Gaussian elimination with partial pivoting. 
}
\usage{
GE(A)
}
\arguments{
\item{A}{An n by p matrix. n must be greater than or equal to p.}
}
\details{This function performs Gaussian elimination to the input matrix and returns the locations of pivoting elements.
}
\value{The output is a vector of length n.
\item{idx}{A vector of length n.  The first p elements are the indices of the pivoting elements, ordered accoridng to columns, and the rest n-p elements are the remaining indices from 1 to n.}
}
