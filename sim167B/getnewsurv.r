## script to generate the ex-post survival info wanted

## first check the sim script and the looping script to make sure that
## the code is corrected for the variable own.id.sg1 This var was
## omitted previously because it was thought redundant.  That was an
## erroneous assumption, since the .mat variables get overwritten in
## part in NAs when the 'anc' style kinship measures are done.

load(".CheckPoint")
source("forward.5yr.altqx.Tcons.sociality.r");   
source("f.qxconditional.r");   
source("f.newmomsurv.r")
## use test.0, the launch state and iterate forward

Nnewcycles <- 20;
test.N <- test;
#test<- test.0;
test$own.id.sg1 <- test$own.id.mat;     #need this fix
qx.by.momalive <- NA + array(0,dim=c(Nnewcycles,16,2,length(subPopTypes) ) );
mom.survival   <- NA + array(0,dim=c(Nnewcycles,2,length(subPopTypes) ) );

pdf(file="expost.pdf",paper="letter",height=7.5,width=10);

for(tt in 1:Nnewcycles){
    test.old <- test;

  test <-
    forward.5yr.altqx(cycleno=as.integer(year),
                     beta.mat=beta.mat,
                     beta.kin=beta.kin,
                     beta.sg1=beta.sg1,
                     own.id.mat=test$own.id.mat,
                     mom.id.mat=test$mom.id.mat,
                     gmom.id.mat=test$gmom.id.mat,
                     ggmom.id.mat=test$ggmom.id.mat,
                     mom.id.kin= test$mom.id.kin,
                     gmom.id.kin=test$gmom.id.kin,
                     ggmom.id.kin=test$ggmom.id.kin,
                     gggmom.id.kin=test$gggmom.id.kin,
                     ggggmom.id.kin=test$ggggmom.id.kin,
                     gggggmom.id.kin=test$gggggmom.id.kin,
                     kinship3.id = test$kinship3.id,
                     kinship4.id = test$kinship4.id,
                     kinship5.id = test$kinship5.id,
                     group.id.anc=test$group.id.anc,
                     own.id.yr= test$own.id.yr,
                     mom.id.yr=  test$mom.id.yr,
                     gmom.id.yr =  test$gmom.id.yr,
                     ggmom.id.yr = test$ggmom.id.yr,
                     gggmom.id.yr=  test$gggmom.id.yr,
                     group.id.mat= test$group.id.mat,
                     gamma.mat= test$gamma.mat,
                     childhood.gamma.mat = test$childhood.gamma.mat,
                     gamma.kin = test$gamma.kin,
                     childhood.gamma.kin = test$childhood.gamma.kin,
                     gamma.anc=test$gamma.anc,
                     childhood.gamma.anc=test$childhood.gamma.anc,
                     own.id.sg1 = test$own.id.sg1, 
                     mom.id.sg1 = test$mom.id.sg1,
                     gmom.id.sg1 = test$gmom.id.sg1,
                     ggmom.id.sg1 = test$ggmom.id.sg1,
                     group.id.sg1 = test$group.id.sg1,
                     gamma.sg1 = test$gamma.sg1,
                     childhood.gamma.sg1 = test$childhood.gamma.sg1,
                     sharing.gamma = test$sharing.gamma,
                     sharing.childhood.gamma = test$sharing.childhood.gamma,             
                     subPop.id = test$subPop.id,
                     gene = test$gene,
                     age = test$age,
                     max.age = max.age,
                     e.fertility = e.fert,
                     e.mortality =  e.mort,
                     e.production =  e.prod,
                     e.cproduction = e.cprod,
                     e.density =  e.dens,
                     #pop.size= length(test$own.id.mat),
                     resource.size = Cresources
                     );

    qx.by.momalive[tt,,,] <- f.qxconditional(test,test.old,usegroup=T);
    
    mom.survival[tt,,] <- f.newmomsurv(test,test.old,usegroup=T)

    print(paste("cycle ",tt));
  }


t.desc<-"Initial" ;
#t.desc<-"Final";
matplot(1:20,qx.by.momalive[,1,2,],type="l",col=subPopTypes,
        xlab="Cycle",ylab="");
title(paste(sim.num,t.desc," 5q0|Mom Dies, by subpopulation "))
legend(x="bottomright",fill=subPopTypes,legend=subPopLabels)

matplot(1:20,qx.by.momalive[,1,1,],type="l",col=subPopTypes,
        xlab="Cycle",ylab="");
title(paste(sim.num,t.desc," 5q0|Mom Alive, by subpopulation "))
legend(x="bottomright",fill=subPopTypes,legend=subPopLabels)


t.denom <- apply(mom.survival, c(1,3),sum); #total women, both categories
t.numer <- mom.survival[,2,];
t.ans  <- t.numer / t.denom;            # proportion of new mothers dying in interval

matplot(1:20,t.ans,type="l",col=subPopTypes,
        xlab="Cycle",ylab="");
title(paste(sim.num,"Proportion of New Mothers dying in next interval\n",t.desc))
legend(x="bottomright",fill=subPopTypes,legend=subPopLabels)
