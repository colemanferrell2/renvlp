\name{contr}
\alias{contr}
\title{Contraction matrix}
\description{
 Generate contraction matrix. 
}
\usage{
contr(d)
}
\arguments{
\item{d}{Dimension of the contraction matrix.  A positive integer.}
}
\details{The contraction and expansion matrices are links between the "vec" operator and "vech"operator: for an d by d symmetric matrix A, vech(A) = contr(d) * vec(A), and vec(A) = expan(d) * vech(A). The "vec" operator stacks the matrix A into an d ^ 2 dimensional vector columnwise.  The "vech" operator stacks the lower triangle or the upper triangle of a symmetric matrix into an d * (d + 1) / 2 vector. For more details of "vec", "vech", contraction and expansion matrix, refer to Henderson and Searle (1979).
}
\value{The output is a matrix.
\item{contrMatrix}{A contraction matrix that has dimension d * (d + 1) / 2 by d ^ 2.}
}
\references{
Henderson, H. V., and Searle, S. R. (1979). Vec and Vech operators for matrices, with some uses in Jacobians and multivariate statistics. Canadian J. Statist. 7, 65 - 81.
}
\examples{
contr(3)
}
