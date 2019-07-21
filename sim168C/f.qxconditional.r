
f.qxconditional <- function(test,test.old){
##

## beginning of interval, number in each age class.

  t.motheralive.tm1 <- (test.old$mom.id.mat %in% test.old$own.id.mat);
  t.motheralive.tm1 <- factor(t.motheralive.tm1,levels=c(T,F),labels=c("MomAlive","MomDead"));
  t.age.factor.old <- factor(test.old$age,levels=0:15);
  t.subpop <-factor(test.old$subPop.id)
# 
  n.tm1 <- tapply(test.old$own.id.mat,
                list(t.age.factor.old,t.motheralive.tm1,t.subpop),
                length);
  n.tm1 <- ifelse(n.tm1==0,NA,n.tm1);  #replace with 0 NA

  
  d.tm1 <- tapply( !(test.old$own.id.mat  %in% test$own.id.mat),
                  list( t.age.factor.old, t.motheralive.tm1,t.subpop),
                  sum,na.rm=T);
                
  d.tm1 <- ifelse(is.na(d.tm1),0,d.tm1);  #replace NAs with 0

  q.tm1 <- d.tm1 / n.tm1 ;

  return(q.bymomalive = q.tm1)
}
