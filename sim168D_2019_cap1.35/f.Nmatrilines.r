
f.Nmatrilines<-function(xtest){
## function to compute the number of matrilines present in each social
## group, as measured by the unique number of ggmom (after excluding
## the moms, gmom, and ggmoms fromt the group members

  vec.matrilines <- NA*unique(xtest$group.id.sg1) ;

  t.uniqids <- sort(unique(xtest$group.id.sg1));
  for(i in seq(t.uniqids) ) {
    isel = (xtest$group.id.sg1 == t.uniqids[i] ) ;
 
    t.own.id=xtest$own.id.sg1[isel]
    t.mom.id=xtest$mom.id.sg1[isel]
    t.gmom.id=xtest$gmom.id.sg1[isel]
    t.ggmom.id=xtest$ggmom.id.sg1[isel]


    ## count own.ids of current generation (c0) but exclude from the list those who are also living mothers, gm, or ggm
    isel.c0 <- (t.own.id %in% t.mom.id)  | (t.own.id %in% t.gmom.id) | (t.own.id %in% t.ggmom.id)

    t.nmatrilines  <- length(table(t.ggnom.id[ ! isel.c0 ])) ;


    vec.relatedness[i] <- t.nmatrilines;
  };                                      #end of for()
  return(vec.relatedness)
}

# end of function

#t.groupsize<- table(test$group.id.sg1);
# vec.relatedness <- f.relatedness(test);
# histogram( ~ vec.relatedness | t.groupsize )
