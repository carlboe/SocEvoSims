f.disperse <- function(x,y,newval){
  # function to split the group identified by x using criterion y and, if satisfied, assign value newval
  # if  splitting is possible, return vector sizeof x with proper assignments;
  # if no splitting is possible, return FALSE

  # examine distinct values in y and split into groups
  y.table <- table(y);
  if(length(y.table) > 1){             # split possible
    y.table<-rev(sort(y.table));
    y.largest <- as.integer(names(y.table[1]));
    y.reassigns <- (y == y.largest);
    y[y.reassigns] <- newval;
    return(y);
  } else {
    return(FALSE)
  };

}
    
