oneK <-
  function(ind,x) {
    n <- nrow(x) ; p <- ncol(x)
    mat1 <- matrix(0,p,p)
    ni <- length(ind)
    if (ni==1) {
      mat1[ind,ind] <- 1/((n-1)*var(x[,ind]))
    }else{
      mat1[ind,ind] <- solve((n-1)*cov(x[,ind]))
    }
    mat1
  }

