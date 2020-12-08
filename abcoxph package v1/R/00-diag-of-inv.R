# Get the diagonal of Q^-1 via inverting
diagOfInv <- function(x, verbose=FALSE,constrA = NULL,i = NULL) {
  if (is.null(i)){
    i <- 1:dim(x)[1]
  }
  result <- diag(solve(x))
  if (is.null(constrA)) {
    vars <- result
  } else {
    # Correct for linear constraints
    WW <- Matrix::solve(x,constrA)
    # VV <- Matrix::solve(t(constrA) %*% WW,WW)
    VV <- t(Matrix::solve(t(WW) %*% constrA,t(WW)))
    # Correct for the ones that are actually being returned
    correction <- rep(0,length(result))
    for (j in i) correction[j] <- VV[j, ] %*% WW[j, ]
    vars <- result - correction
  }
  vars[i]
}

