## f.subsamptest.R

f.subsamptest <- function(x,sampsize){

## subsample the population down to a size sampsize.  The object x 
## contains many components which need to be subsampled, skipped, or
## specially handled

  n <- length(x$own.id.mat);
  mark <- sort(sample(seq(x$own.id.mat), sampsize)); #
   
  names.x <- names(x);
  
  for( i in seq(names.x)){
    if(length(x[[i]]) ==n){
      ## we have a component to subsample
      x[[i]] <- (x[[i]])[mark] ;
    };
  }

  # handle the gene case separately
  x$gene <- x$gene[mark,];

  return(x)

}
    
