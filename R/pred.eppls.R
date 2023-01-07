pred.eppls <- function(m, X1new, X2new) {
    
    if(is.null(m$covMatrix2)) stop("The asymptotic distribution part of eppls is missing. Rerun eppls with asy = T")
    
    r <- ncol(m$SigmaYcX)
    n <- m$n
    if (is.null(m$Gamma)) {
        u <- 0
    } else {
        u <- ncol(m$Gamma)
    }
    
    if (u == 0) {
        value <- m$muY + m$beta2 %*% X2new
        temp <- kronecker(X2new, diag(r))
        temp2 <- crossprod(temp, m$covMatrix2) %*% temp / n
        covMatrix.estm <- m$SigmaYcX / n + temp2
        covMatrix.pred <- (1 + 1 / n) * m$Sigma + temp2
    } else {
        X1new <- as.matrix(X1new)
        X2new <- as.matrix(X2new)
        if (nrow(X1new) == 1) X1new <- t(X1new)
        if (nrow(X2new) == 1) X2new <- t(X2new)
       
        value <- m$muY + crossprod(m$beta1, X1new) + crossprod(m$beta2, X2new)
        temp1 <- kronecker(diag(r), X1new)
        temp2 <- kronecker(diag(r), X2new)
        covMatrix.estm <- m$SigmaYcX / n + crossprod(temp1, m$covMatrix1) %*% temp1 / n + crossprod(temp2, m$covMatrix2) %*% temp2 / n
        covMatrix.pred <- (1 + 1 / n) * m$SigmaYcX + crossprod(temp1, m$covMatrix1) %*% temp1 / n + crossprod(temp2, m$covMatrix2) %*% temp2 / n
    }
    
    SE.estm <- sqrt(diag(covMatrix.estm))
    SE.pred <- sqrt(diag(covMatrix.pred))
    return(
    list(value = value, covMatrix.estm = covMatrix.estm, SE.estm = SE.estm, covMatrix.pred = covMatrix.pred, SE.pred = SE.pred
    )
    )
    
    
}
