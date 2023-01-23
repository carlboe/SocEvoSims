

lifetable.5qx <- function(nqx){
            # Derive lifetable values from single year mortality rates
  n <- 5
  len.m <- length(nqx) - 1
  age <- seq(from=0,by=n,length= len.m+1)
  nax <- c(2.5, rep(2.5, len.m) )     # 5a0 is 1.15 approx; e.g. Swedish Female 1751, HMD
  nqx[nqx > 1] <- 1
  lx<- c(1,cumprod(1-nqx[1:len.m]))
  ndx <- -diff(c(lx,0))
  nLx <- (lx * n) - (ndx * (n - nax))
  Tx <- rev(cumsum(rev(nLx)))
  ex <- Tx/lx
  result <- data.frame(age = age, nqx = nqx, lx = lx, ndx = ndx,
                       nLx = nLx, Tx = Tx, ex = ex);
  return(result)
}

