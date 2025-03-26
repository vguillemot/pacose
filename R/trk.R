trk <-
function(ind,x) {
  ni <- length(ind)
  n <- nrow(x)
  if (ni==1) {
    return(-(1/((n-1)*var(x[,ind])))**2)
  }else{
    mattemp <- solve((n-1)*cov(x[,ind]))
    return(-1/2*tr(mattemp**2) - 1/2*(tr(mattemp))**2 )
  }
}

