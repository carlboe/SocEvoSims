
f.relatedness<-function(xtest){
# function to compute relatedness

vec.relatedness <- NA*unique(xtest$group.id.sg1) ;

t.uniqids <- sort(unique(xtest$group.id.sg1));
for(i in seq(t.uniqids) ) {
  isel = (xtest$group.id.sg1 == t.uniqids[i] ) ;
 
  t.own.id=xtest$own.id.sg1[isel]
  t.mom.id=xtest$mom.id.sg1[isel]
  t.gmom.id=xtest$gmom.id.sg1[isel]
  t.ggmom.id=xtest$ggmom.id.sg1[isel]


  ## count own.ids of current generation (c0) but exclude from the list those who are also living mothers, gm, or ggm
  isel.c0 <- (t.own.id %in% t.mom.id)  | (t.own.id %in% t.gmom.id) | (t.own.id %in% t.ggmom.id)

  t.c0  <- t.own.id[ ! isel.c0 ] ;


  ## count one generation back (c1) but exclude from the list those who are also living  gm, or ggm
  isel.c1 <-  (t.own.id %in% t.gmom.id) | (t.own.id %in% t.ggmom.id)
  t.c1  <- t.mom.id[ ! isel.c1 ] ;

  ## count two generation back (c2) but exclude from the list those who are also living  ggm
  isel.c2 <-   (t.own.id %in% t.ggmom.id)
  t.c2  <- t.gmom.id[ ! isel.c2 ] ;

  ## count three generation back (c3)
  t.c3  <- t.ggmom.id;

  t.nnodes <- length(table(t.c0)) + length(table(t.c1)) + length(table(t.c2))+ length(table(t.c3)) ;

  ## minimum possible nodes has all living group members related, possibly to the c4 level
  t.minnodes <- length(t.own.id) ;        # mom, gm, ggm present and living in group
  t.minnodes < - max(t.minnodes, 4);      # with c0 size of 1, c1,c2,c3 nodes give 4 total
  ## max has none related, separate lines at all generates
  t.maxnodes <- 4*length(t.own.id);
  if(t.maxnodes < t.nnodes | t.minnodes > t.nnodes) stop("Impossible computation ?")

  ## related goes from 0 to 1 with 0 being completely unrelated. 
  t.relatedness <- ifelse(t.maxnodes-t.minnodes==0,1, 1 - (t.nnodes - t.minnodes)/(t.maxnodes - t.minnodes) );

  vec.relatedness[i] <- t.relatedness;
};                                      #end of for()
  return(vec.relatedness)
}

# end of function

#t.groupsize<- table(test$group.id.sg1);
# vec.relatedness <- f.relatedness(test);
# histogram( ~ vec.relatedness | t.groupsize )
