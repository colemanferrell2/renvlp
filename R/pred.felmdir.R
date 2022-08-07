pred.felmdir <- function(X, Y, ux, uy, t1, t2, Xnew, knots = c(0, 0.25, 0.5, 0.75, 1)) {

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
  
  m <- stenv(X.cord, Y.cord, ux, uy)
  
  X.cord.new <- stats::lm(Xnew ~ basis.value.X -1)$coefficients
  pred.env <- t(m$beta) %*% as.matrix(X.cord.new) + m$mu
  tmppred <- kronecker(diag(dx), t(X.cord.new))
  var.est.env.cord <- (tmppred %*% m$covMatrix %*% t(tmppred) + m$SigmaYcX) / n
  var.pred.env.cord <- (tmppred %*% m$covMatrix %*% t(tmppred) + m$SigmaYcX) / n + m$SigmaYcX
  
  value <- t(pred.env) %*% t(basis.value.Y)
  covMatrix.estm <- basis.value.Y %*% var.est.env.cord %*% t(basis.value.Y)
  covMatrix.pred <- basis.value.Y %*% var.pred.env.cord %*% t(basis.value.Y)
  SE.estm <- sqrt(diag(covMatrix.estm))
  SE.pred <- sqrt(diag(covMatrix.pred))
  
  return(
    list(
      value = value, covMatrix.estm = covMatrix.estm, SE.estm = SE.estm, covMatrix.pred = covMatrix.pred, SE.pred = SE.pred
    )
  )
  
}
