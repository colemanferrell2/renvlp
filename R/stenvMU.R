stenvMU <- function(X, Y, q, u){
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
  if (q > p || q < 0) 
    stop("q must be an integer between 0 and p.")
  if (sum(duplicated(cbind(X, Y), MARGIN = 2)) > 0) 
    stop("Some responses also appear in the predictors, or there maybe duplicated columns in X or Y.")
  sigY <- stats::cov(Y) * (n - 1)/n
  sigYX <- stats::cov(Y, X) * (n - 1)/n
  sigX <- stats::cov(X) * (n - 1)/n
  invsigX <- chol2inv(chol(sigX))
  invsigY <- chol2inv(chol(sigY))
  betaOLS <- sigYX %*% invsigX
  U1 <- crossprod(sigYX, invsigY) %*% sigYX
  M1 <- sigX - U1
  tmp1 <- envMU(M1, U1, q)
  Phi <- tmp1$Gammahat
  Phi0 <- tmp1$Gamma0hat
  E1 <- crossprod(Phi, sigX) %*% Phi
  invE1 <- chol2inv(chol(E1))
  E2 <- sigYX %*% Phi
  U2 <- E2 %*% tcrossprod(invE1, E2)
  M2 <- sigY - U2
  tmp2 <- envMU(M2, U2, u)
  Ga <- tmp2$Gammahat
  Ga0 <- tmp2$Gamma0hat
  m1 <- crossprod(Phi, sigX) %*% Phi
  m2 <- crossprod(Phi0, sigX) %*% Phi0
  m3 <- crossprod(Ga0, sigY) %*% Ga0
  m4 <- crossprod(Ga, M2) %*% Ga
  e1 <- eigen(m1)
  e2 <- eigen(m2)
  e3 <- eigen(m3)
  e4 <- eigen(m4)
  obj1 <- (sum(log(e1$values)) + sum(log(e2$values)) 
           + sum(log(e3$values)) + sum(log(e4$values))) 
  
  maxiter = 100
  ftol = 0.0001
  i <- 1
  while(i < maxiter){
    E3 <- crossprod(Ga, sigY) %*% Ga
    invE3 <- chol2inv(chol(E3))
    E4 <- crossprod(sigYX, Ga)
    U3 <- E4 %*% tcrossprod(invE3, E4)
    M3 <- sigX - U3
    tmp3 <- envMU(M3, U3, q)
    Phi <- tmp3$Gammahat
    Phi0 <- tmp3$Gamma0hat
    E5 <- crossprod(Phi, sigX) %*% Phi
    invE5 <- chol2inv(chol(E5))
    E6 <- sigYX %*% Phi
    U4 <- E6 %*% tcrossprod(invE5, E6)
    M4 <- sigY - U4
    tmp4 <- envMU(M4, U4, u)
    Ga <- tmp4$Gammahat
    Ga0 <- tmp4$Gamma0hat
    m5 <- crossprod(Phi, sigX) %*% Phi
    m6 <- crossprod(Phi0, sigX) %*% Phi0
    m7 <- crossprod(Ga0, sigY) %*% Ga0
    m8 <- crossprod(Ga, M4) %*% Ga
    e5 <- eigen(m5)
    e6 <- eigen(m6)
    e7 <- eigen(m7)
    e8 <- eigen(m8)
    obj2 <- (sum(log(e5$values)) + sum(log(e6$values)) 
             + sum(log(e7$values)) + sum(log(e8$values)))
    if (abs(obj1 - obj2) < ftol) {
      break
    }
    else {
      obj1 <- obj2
      i <- i + 1
    }
  }
    return(list(Gammahat = Ga, Gamma0hat = Ga0,
                Phihat = Phi, Phi0hat = Phi0, objfun = obj2))
}  
