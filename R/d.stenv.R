d.stenv <- function (X, Y, alpha = 0.01){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  if (a[1] != nrow(X)) 
    stop("X and Y should have the same number of observations.")
  if (sum(duplicated(cbind(X, Y), MARGIN = 2)) > 0) 
    stop("Some responses also appear in the predictors, or there maybe duplicated columns in X or Y.")
  d.test <- function (X, Y, d) {
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    a <- dim(Y)
    n <- a[1]
    r <- a[2]
    p <- ncol(X)
    m <- min(r, p)
    sigY <- stats::cov(Y) * (n - 1)/n
    sigYX <- stats::cov(Y, X) * (n - 1)/n
    sigX <- stats::cov(X) * (n - 1)/n
    invsigX <- chol2inv(chol(sigX))
    betaOLS <- tcrossprod(invsigX, sigYX)
    U <- sigYX %*% betaOLS
    M <- sigY - U
    E1 <- eigen(sigX)
    V1 <- E1$values
    Q1 <- E1$vectors
    D1 <- diag(sqrt(V1), ncol = length(V1))
    P1 <- Q1 %*% tcrossprod(D1, Q1)
    E2 <- eigen(M)
    V2 <- E2$values
    Q2 <- E2$vectors
    D2 <- diag(1 / sqrt(V2), ncol = length(V2))
    P2 <- Q2 %*% tcrossprod(D2, Q2)
    betastd <- P1 %*% betaOLS %*% P2 * sqrt((n - p - 1)/n)
    E3 <- svd(betastd)
    Rho <- E3$d
    R <- Rho^2
    b <- d + 1
    V <- R[b : m]
    test <- sum(V) * n
    return(list(lambda = test))
  }
  m <- min(r, p)
  m1 <- m - 1
  loglik.seq <- unlist(lapply(0 : m1, function(x) d.test(X, Y, x)$lambda))
  npara.seq <- (p - (0 : m1)) * (r - (0 : m1))
  lrt.test <- stats::pchisq(loglik.seq[1 : m], npara.seq[1 : m], lower.tail = F)
  if (any(lrt.test > alpha)) {
    u.lrt <- which(lrt.test > alpha)[1] - 1
  }
  else {
    u.lrt <- m
  }
  return(list(rank.beta = u.lrt))
}