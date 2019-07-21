


f.newmomsurv <- function(test,test.old,usegroup=FALSE,subpops=seq(subPopTypes)){
##
## look at mothers of today's babies and determine whether or not they make it during the cycle
  isel.isbaby.t <- (test$age == 0 );
  t.momid.baby.t <- test$mom.id.sg1[isel.isbaby.t]; # all new mothers, including those who died
  
  t.newmom.lived <-    ( t.momid.baby.t %in% test$own.id.sg1 );
  t.newmom.lived <-  factor(t.newmom.lived, levels=c(T,F),labels=c('MomLives','MomDies') );

  t.isel <-  ( test.old$own.id.sg1 %in%  t.momid.baby.t  );

  t.ageofnewmom <- test.old$age[ t.isel  ];
  # classify new moms using pop membership of their babies.
  if(usegroup){
    t.subpop <- factor(test$subPop.id[isel.isbaby.t],levels=seq(subPopTypes)); 
    t.momcounts <- tapply(t.newmom.lived,list(t.newmom.lived,t.subpop),length);
  } else {
    t.momcounts <- tapply(t.newmom.lived,list(t.newmom.lived),length);
  };
  
  t.momcounts <- ifelse(is.na(t.momcounts),0,t.momcounts);  

  
  return(t.momcounts)
}
