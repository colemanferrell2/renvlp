u.felmdir <- function(X, Y, t1, t2, knots = c(0, 0.25, 0.5, 0.75, 1)) {
  
  X <- as.matrix(X)
  Y <- as.matrix(Y) 
  dimX <- dim(X)
  n <- dimX[1]
  
  if (min(knots) < 0) stop("Location of the knots cannot be negative.")
  if (is.null(knots)) {
    knots <- orthogonalsplinebasis::expand.knots(c(0, 0.25, 0.5, 0.75, 1), order = 4)
  } else {
    knots <- knots / max(knots)
    knots <- orthogonalsplinebasis::expand.knots(knots, order = 4)
  }
  
  if (nrow(Y)!= nrow(X)) stop("X and Y should have the same number of observations.")
  if (min(t1) < 0) stop("Time points cannot be negative.")
  if (min(t2) < 0) stop("Time points cannot be negative.")
  if (max(t1) > 1) t1 <- t1 / max(t1)
  if (max(t2) > 1) t2 <- t2 / max(t2)
  
  spbasis <- orthogonalsplinebasis::SplineBasis(knots, order = 4)
  Gx <- orthogonalsplinebasis::GramMatrix(spbasis)
  eigX <- eigen(Gx)
  Gx12 <- eigX$vectors %*% diag(sqrt(eigX$values)) %*% t(eigX$vectors)
  Gx12inv <- eigX$vectors %*% diag(1 / sqrt(eigX$values)) %*% t(eigX$vectors)
  dx <- dim(Gx)[1]
  
  X.cord <- matrix(0, n, dx)
  Y.cord <- matrix(0, n, dx)
  basis.value.X <- orthogonalsplinebasis::evaluate(spbasis, t1) %*% Gx12inv
  basis.value.Y <- orthogonalsplinebasis::evaluate(spbasis, t2) %*% Gx12inv
  
  for (i in 1:n) {
      X.cord[i,] <- stats::lm(X[i, ] ~ basis.value.X -1)$coefficients
      Y.cord[i,] <- stats::lm(Y[i, ] ~ basis.value.Y -1)$coefficients
  }
  
  u.mat <- u.stenv(X.cord, Y.cord, alpha = 0.01)
  ux <- u.mat$u.bic[1]
  uy <- u.mat$u.bic[2]
  m <- stenv(X.cord, Y.cord, ux, uy)
  mfull <- stenv(X.cord, Y.cord, dx, dx)
  return(list(ux = ux, uy = uy, beta = m$beta, betafull = mfull$beta, alpha = m$mu, alphafull = mfull$mu))
}
