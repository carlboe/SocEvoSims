
lifetable.1qx <- function(nqx){
            # Derive lifetable values from single year mortality rates
            n <- 1
                      len.m <- length(nqx) - 1
                      age <- seq(0, len.m)
                      nax <- c(0.5, rep(0.5, len.m) )
                      nqx[nqx > 1] <- 1
                      lx<- c(1,cumprod(1-nqx[1:len.m]))
                      ndx <- -diff(c(lx,0))
                      nLx <- (lx * n) - (ndx * (n - nax))
                      Tx <- rev(cumsum(rev(nLx)))
                      ex <- Tx/lx
                      result <- data.frame(age = age, nqx = nqx, lx = lx, ndx = ndx,
                                           nLx = nLx, Tx = Tx, ex = ex)
                      return(result)
          }

