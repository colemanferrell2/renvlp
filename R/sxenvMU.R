sxenvMU <- function(X, Y, u, R){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  if (a[1] != nrow(X)) 
    stop("X and Y should have the same number of observations.")
  if (u > p || u < 0) 
    stop("u must be an integer between 0 and p.")
  if (sum(duplicated(cbind(X, Y), MARGIN = 2)) > 0) 
    stop("Some responses also appear in the predictors, or there maybe duplicated columns in X or Y.")
  sigY <- stats::cov(Y) * (n - 1)/n
  sigYX <- stats::cov(Y, X) * (n - 1)/n
  sigX <- stats::cov(X) * (n - 1)/n
  invsigX <- chol2inv(chol(sigX))
  invsigY <- chol2inv(chol(sigY))
  t1 <- crossprod(sigYX, invsigY)
  U <- t1 %*% sigYX
  M <- sigX - U
  tmp <- envMU(M, U, u)
  Gamma <- tmp$Gammahat
  if (length(R) == p) {
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
      t1 <- crossprod(sigYX, invsigY)
      U <- t1 %*% sigYX
      M <- sigX - U
      d1 <- c(1, d)
      Lambda <- diag(d1)
      invLambda <- diag(1 / d1)
      m1 <- crossprod(Gamma, Lambda)
      eigtem1 <- eigen(m1 %*% tcrossprod(invsigX, m1))
      m2 <- crossprod(Gamma, invLambda)
      eigtem2 <- eigen(m2 %*% tcrossprod(M, m2))
      temp1 <- sum(log(eigtem1$values))
      temp2 <- sum(log(eigtem2$values))
      objfun <- temp1 + temp2
      return(objfun)
    }
    
  d.init <- rep(1, (p - 1))
  
  k2 <- rep(0, (p - 1))
  tmp.init <- Rsolnp::solnp(pars = d.init, fun = objfun, LB = k2,
                  control = list(delta = 1e-8, tol = 1e-6, trace = 0),
                  Gamma = Gamma, X = X, Y = Y)
  d <- tmp.init$pars
  
  obj1 <- objfun(d, Gamma, X, Y)
  
  maxiter = 20000
  ftol = 1e-06
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
  } else {
    objfun1 <- function(d, Gamma, X, Y, R){
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
      t1 <- crossprod(sigYX, invsigY)
      U <- t1 %*% sigYX
      M <- sigX - U
      Lambda <- diag(d)
      invLambda <- diag(1 / d)
      m1 <- crossprod(Gamma, Lambda)
      eigtem1 <- eigen(m1 %*% tcrossprod(invsigX, m1))
      m2 <- crossprod(Gamma, invLambda)
      eigtem2 <- eigen(m2 %*% tcrossprod(M, m2))
      temp1 <- sum(log(eigtem1$values))
      temp2 <- sum(log(eigtem2$values))
      objfun <- temp1 + temp2
      return(objfun)
    }
    
    heq <- function (d, Gamma, X, Y, R){
      iter = sum(R)
      i <- 1
      cont <- NULL
      while (i < iter) {
        C <- matrix(0, (R[i] - 1), sum(R))
        for (j in 2 : (sum(R) - 1)){
          if (R[i] == j) {
            for (k in 1 : (j - 1)){
              s <- sum(R[0 : (i - 1)]) + 1
              C[k , s] <- 1
              C[k , (s + k)] <- -1
            }      
          }
        }
        cont <- rbind(cont, C)
        if (i == length(R)) {
          c1 <- c(1, rep(0, (sum(R) - 1)))
          cont <- rbind(c1, cont)
          break
        } else {
          i <- i + 1
        }
      }
      
      g <- rep(NA, nrow(cont))
      g[1] <- d[1]
      for (m in 2 : nrow(cont)){
        g[m] <- d[which(cont[m, ] == 1)] - d[which(cont[m, ] == -1)]  
      }
      g
    }
    d.init <- rep(1, p) 
    g1 <- heq(d.init, Gamma, X, Y, R)
    k1 <- c(1, rep(0, length(g1) - 1))
    k2 <- rep(0, p)
    
    tmp.init <- Rsolnp::solnp(pars = d.init, fun = objfun1, eqfun = heq, eqB = k1,
                   LB = k2, 
                   control = list(delta = 1e-8, tol = 1e-6, trace = 0),
                   R = R, Gamma = Gamma, X = X, Y = Y)
    d <- tmp.init$pars
    obj1 <- objfun1(d, Gamma, X, Y, R)
    
    maxiter = 20000
    ftol = 1e-05
    i <- 1
    
    while (i < maxiter) {
      L1 <- diag(d)
      invL1 <- diag(1 / d)
      M1 <- invL1 %*% M %*% invL1
      U1 <- invL1 %*% U %*% invL1
      tmp1 <- envMU(M1, U1, u)
      Gamma <- tmp1$Gammahat
      temp1 <- Rsolnp::solnp(pars = d, fun = objfun1, eqfun = heq, eqB = k1, 
              LB = k2, control = list(delta = 1e-8, tol = 1e-6, trace = 0),
              R = R, Gamma = Gamma, X = X, Y = Y)
      d <- temp1$pars
      obj2 <- objfun1(d, Gamma, X, Y, R)
      if (abs(obj1 - obj2) < ftol) {
        break
      }
      else {
        obj1 <- obj2
        i <- i + 1
      }
    }
    Gamma0 <- tmp1$Gamma0hat
    Lambda <- diag(d)
  }

  return(list(Gammahat = Gamma, Gamma0hat = Gamma0, Lambdahat = Lambda,
              objfun = obj2))
}