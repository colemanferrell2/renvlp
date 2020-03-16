logit.env <- function(X, Y, u, asy = TRUE, init = NULL) {
  Y <- as.matrix(Y)
  X <- as.matrix(X)
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
  if (!is.null(init)) {
    if (nrow(init) != p || ncol(init) != u) stop("The dimension of init is wrong.")
  }
  
  sigY <- stats::cov(Y) * (n - 1)/n
  sigYX <- stats::cov(Y, X) * (n - 1)/n
  sigX <- stats::cov(X) * (n - 1)/n
  invsigY <- chol2inv(chol(sigY))
  tmp <- logit.envMU(X, Y, u, initial = init)
  
  
  Gammahat <- tmp$Gammahat
  Gamma0hat <- tmp$Gamma0hat
  objfun <- tmp$objfun
  muhat <- tmp$muhat
  etahat <- tmp$etahat
  weighthat <- tmp$weighthat
  Vhat <- tmp$Vhat
  avarhat <- tmp$avar
  
  covMatrix <- NULL
  asySE <- NULL
  ratio <- NULL
  if (u == 0) {
    Omegahat <- NULL
    Omega0hat <- sigX
    betahat <- matrix(0, p, 1)
    SigmaXhat <- sigX
    loglik <- objfun
    if (asy == T) 
      ratio <- matrix(1, p, r)
  }
  else if (u == p) {
    invsigX <- chol2inv(chol(sigX))
    Omegahat <- sigX
    Omega0hat <- NULL
    betahat <- etahat
    SigmaXhat <- sigX
    loglik <- objfun
    if (asy == T) {
      covMatrix <- avarhat
      asySE <- matrix(sqrt(diag(covMatrix)), nrow = p)
      ratio <- matrix(1, p, r)
    }
  }
  else {
    invsigX <- chol2inv(chol(sigX))
    Omegahat <- crossprod(Gammahat, sigX) %*% Gammahat
    Omega0hat <- crossprod(Gamma0hat, sigX) %*% Gamma0hat
    betahat <- Gammahat %*% etahat
    SigmaXhat <- Gammahat %*% tcrossprod(Omegahat, 
                 Gammahat) + Gamma0hat %*% tcrossprod(Omega0hat, Gamma0hat)
    betaOLS <- tcrossprod(invsigX, sigYX)
    PGamma <- tcrossprod(Gammahat)
    loglik <- objfun
    if (asy == T) {
      options(warn = -1) 
      fit <- stats::glm(Y ~ X, family = stats::binomial(link = "logit"))
      mu <- as.vector(fit$coefficients[1])
      beta <- as.vector(fit$coefficients[2 : (p + 1)])
      theta <- matrix(1, n, 1) %*% mu + X %*% beta
      c.theta <- - exp(theta) / (1 + exp(theta))^2
      c.theta.mean <- mean(c.theta)
      weight <- c.theta / c.theta.mean
      mu1 <- exp(theta) / (1 + exp(theta))
      V <- theta + ((Y - mu1) / weight)
      W <- diag(as.vector(weight))
      wx <- W %*% X
      mean.wx <- apply(wx, 2, mean)
      wxx <- X - matrix(1, nrow = n) %*% mean.wx
      sigwx <- crossprod(wxx, W) %*% wxx / n
      wv <- W %*% V
      mean.wv <- mean(wv)
      wvv <- V - matrix(1, nrow = n) %*% mean.wv
      sigwxv <- crossprod(wxx, W) %*% wvv / n
      inv.sigwx <- chol2inv(chol(sigwx))
      avar <- inv.sigwx / (- c.theta.mean)
      covMatrix <- avar
      asyFm <- matrix(sqrt(diag(covMatrix)), nrow = p)
      PGamma <- tcrossprod(Gammahat)
      t1 <- avarhat
      T1 <- PGamma %*% t1 %*% PGamma
      invt1 <- chol2inv(chol(t1))
      t2 <- kronecker(etahat, t(Gamma0hat))
      invOme <- chol2inv(chol(Omegahat))
      invOme0 <- chol2inv(chol(Omega0hat))
      t3 <- kronecker(Omegahat, invOme0) + kronecker(invOme, 
            Omega0hat) - 2 * kronecker(diag(u), diag(p - u))
      M <- t2 %*% tcrossprod(invt1, t2) + t3
      invM <- chol2inv(chol(M))
      T2 <- crossprod(t2, invM) %*% t2
      covMatrix <- T1 + T2
      asySE <- matrix(sqrt(diag(covMatrix)), nrow = p)
      ratio <- asyFm/asySE
    }
  }
  return(list(Gamma = Gammahat, Gamma0 = Gamma0hat, mu = muhat, 
              beta = as.matrix(betahat), SigmaX = SigmaXhat, eta = etahat, 
              Omega = Omegahat, Omega0 = Omega0hat, 
              loglik = loglik, n = n, covMatrix = covMatrix, asySE = asySE, 
              ratio = ratio))
}
