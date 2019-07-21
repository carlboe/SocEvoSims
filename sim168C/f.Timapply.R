###################################################
# tapply function which includes all age categories
###################################################

timapply <- function (x,y,...,my.length=16){
  x <- c(x,rep(0,my.length))
  y <- c(y,seq(0,(my.length-1)))
  result <- tapply(x,y,...)
  return(result)}

# variant of timapply which does not add 0s to each age
timapply2 <- function(x,y,...,my.length=max.age+1){
  t.tmp=tapply(x,y,...);
  t.result=seq(0,my.length-1) + NA;
  t.matchups = seq(0,my.length-1) %in% as.numeric(labels(t.tmp)[[1]] );
  t.result[ t.matchups ] = t.tmp;
  t.result[ ! t.matchups ] = 0;
  return(t.result)
};


Timapply <- function (x,y,...,my.length=16){
  result <- tapply(x,y,...)
  mark <- 1+ as.numeric(names(result))
  final.result <- rep(NA,my.length)
  names(final.result) <- seq(0,,5,my.length)
  final.result[mark] <- result
  return(final.result)}

# variant which preserves (as 0 count) y entries which are not in INDEX x

timapply3 <- function(x,y,...,y1=NULL, y2=NULL){
  t.tmp=tapply(x,y,...);
  t.result=seq(0,my.length-1) + NA;
  t.matchups = seq(0,my.length-1) %in% as.numeric(labels(t.tmp)[[1]] );
  t.result[ t.matchups ] = t.tmp;
  t.result[ ! t.matchups ] = 0;
  return(t.result)
};

     n <- 17; fac <- factor(rep(1:3, len = n), levels = 1:5)
     table(fac)
     tapply(1:n, fac, sum)
     tapply(1:n, fac, sum, simplify = FALSE)
     tapply(1:n, fac, range)
     tapply(1:n, fac, quantile)

     ## example of ... argument: find quarterly means
     tapply(presidents, cycle(presidents), mean, na.rm = TRUE)

     ind <- list(c(1, 2, 2), c("A", "A", "B"))
     table(ind)
     tapply(1:3, ind) #-> the split vector
     tapply(1:3, ind, sum)

