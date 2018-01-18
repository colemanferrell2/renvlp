senvMU <- function(X, Y, u){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  if (a[1] != nrow(X)) 
    stop("X and Y should have the same number of observations.")
  if (u > r || u < 0) 
    stop("u must be an integer between 0 and r.")
  if (sum(duplicated(cbind(X, Y), MARGIN = 2)) > 0) 
    stop("Some responses also appear in the predictors, or there maybe duplicated columns in X or Y.")
  objfun <- function(d, Gamma, X, Y){
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    a <- dim(Y)
    n <- a[1]
    r <- a[2]
    p <- ncol(X)
    sigY <- stats::cov(Y) * (n - 1)/n
    sigYX <- stats::cov(Y, X) * (n - 1)/n
    sigX <- stats::cov(X) * (n - 1)/n
    invsigX <- chol2inv(chol(sigX))
    invsigY <- chol2inv(chol(sigY))
    betaOLS <- sigYX %*% invsigX
    U <- tcrossprod(betaOLS, sigYX)
    M <- sigY - U
    d1 <- c(1, d)
    Lambda <- diag(d1)
    invLambda <- diag(1 / d1)
    m1 <- crossprod(Gamma, Lambda)
    eigtem1 <- eigen(m1 %*% tcrossprod(invsigY, m1))
    m2 <- crossprod(Gamma, invLambda)
    eigtem2 <- eigen(m2 %*% tcrossprod(M, m2))
    temp1 <- sum(log(eigtem1$values))
    temp2 <- sum(log(eigtem2$values))
    objfun <- temp1 + temp2
    return(objfun)
  }
  
  sigY <- stats::cov(Y) * (n - 1)/n
  sigYX <- stats::cov(Y, X) * (n - 1)/n
  sigX <- stats::cov(X) * (n - 1)/n
  invsigX <- chol2inv(chol(sigX))
  invsigY <- chol2inv(chol(sigY))
  betaOLS <- sigYX %*% invsigX
  U <- tcrossprod(betaOLS, sigYX)
  M <- sigY - U
  tmp <- envMU(M, U, u)
  Gamma <- tmp$Gammahat
  
  d.init <- rep(1, (r - 1))
  k2 <- rep(0, (r - 1))
  
  tmp.init <- Rsolnp::solnp(pars = d.init, fun = objfun, LB = k2,
                    control = list(delta = 1e-8, tol = 1e-6, trace = 0),
                    Gamma = Gamma, X = X, Y = Y)
  d <- tmp.init$pars
  obj1 <- objfun(d, Gamma, X, Y)

  
  maxiter = 20000
  ftol = 1e-04
  i <- 1
  
  while (i < maxiter) {
    d1 <- c(1, d)
    L1 <- diag(d1)
    invL1 <- diag(1 / d1)
    M1 <- invL1 %*% M %*% invL1
    U1 <- invL1 %*% U %*% invL1
    tmp1 <- envMU(M1, U1, u)
    Gamma <- tmp1$Gammahat
    temp1 <- Rsolnp::solnp(pars = d, fun = objfun, LB = k2, 
                  control = list(delta = 1e-8, tol = 1e-6, trace = 0),
                  Gamma = Gamma, X = X, Y = Y)
    d <- temp1$pars
    obj2 <- objfun(d, Gamma, X, Y)
    if (abs(obj1 - obj2) < ftol) {
      break
    }
    else {
      obj1 <- obj2
      i <- i + 1
    }
  }
  Gamma0 <- tmp1$Gamma0hat
  d1 <- c(1, d)
  Lambda <- diag(d1)
  return(list(Gammahat = Gamma, Gamma0hat = Gamma0, Lambdahat = Lambda,
              objfun = obj2))
}