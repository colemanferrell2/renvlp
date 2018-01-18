\name{cv.logit.env}
\alias{cv.logit.env}
\title{Cross validation for logit.env}
\description{
 Compute the prediction error for the envelope estimator in logistic regression using m-fold cross validation.
}
\usage{
cv.logit.env(X, Y, u, m, nperm)
}
\arguments{
  \item{X}{Predictors. An n by p matrix, p is the number of predictors and n is number of observations.  The predictors must be continuous variables.}
  \item{Y}{Response. An n by 1 matrix. The univariate response must be binary.}
  \item{u}{Dimension of the envelope. An integer between 0 and p.}
  \item{m}{A positive integer that is used to indicate m-fold cross validation.}
  \item{nperm}{A positive integer indicating number of permutations of the observations, m-fold cross validation is run on each permutation.}
}
\details{
This function computes prediction errors using m-fold cross validation.  For a fixed dimension u, the data is randomly partitioned into m parts, each part is in turn used for testing for the prediction performance while the rest m-1 parts are used for training.  This process is repeated for \code{nperm} times, and average prediction error is reported. 
}
\value{The output is a real nonnegative number.
\item{cvPE}{The prediction error estimated by m-fold cross validation.} 
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
Y <- as.numeric(ifelse(horseshoecrab[ , 4] > 0, 1, 0))

m <- 5
nperm <- 50
\dontrun{cvPE <- cv.logit.env(X, Y, 1, m, nperm)}
\dontrun{cvPE}

}

