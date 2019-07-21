# function to dump the storage modes of test

f.showmode <- function(x){
  names.x <- names(x);
  for(i in seq(names.x)){
    print(paste(names.x[i], storage.mode( x[[i]])));
  }

}

f.setmode <- function(x){
  names.x <- names(x);
  int.list <- c("own.id.mat","mom.id.mat","gmom.id.mat","ggmom.id.mat","mom.id.kin","gmom.id.kin",
 "ggmom.id.kin","gggmom.id.kin","ancestor.id","own.id.yr","mom.id.yr","gmom.id.yr",
 "ggmom.id.yr","gggmom.id.yr","group.id.mat","own.id.sg1","mom.id.sg1","gmom.id.sg1",
 "ggmom.id.sg1","group.id.sg1","ages","genes");
 
  for(i in seq(names.x)){
    if(names.x[i] %in% int.list){
      print(paste(i,names.x[i]));
      storage.mode(x[[i]]) <- "integer";
    }
    
  }
  return(x)
}
  

