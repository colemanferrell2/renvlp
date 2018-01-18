\name{expan}
\alias{expan}
\title{Expansion matrix}
\description{
 Generate expansion matrix. 
}
\usage{
expan(d)
}
\arguments{
\item{d}{Dimension of the expansion matrix.  A positive integer.}
}
\details{The contraction and expansion matrices are links between the "vec" operator and "vech"operator: for an d by d symmetric matrix A, vech(A) = contr(d) * vec(A), and vec(A) = expan(d) * vech(A). The "vec" operator stacks the matrix A into an d ^ 2 dimensional vector columnwise.  The "vech" operator stacks the lower triangle or the upper triangle of a symmetric matrix into an d * (d + 1) / 2 vector. For more details of "vec", "vech", contraction and expansion matrix, refer to Henderson and Searle (1979).
}
\value{The output is a matrix.
\item{expanMatrix}{An expansion matrix that has dimension d ^ 2 by d * (d + 1) / 2.}
}
\references{
Henderson, H. V., and Searle, S. R. (1979). Vec and Vech operators for matrices, with some uses in Jacobians and multivariate statistics. Canadian J. Statist. 7, 65 - 81.
}
\examples{
expan(3)
}