## sim167B

## investigation of the 1st 1000 cycles of sim16x
load(".CheckPoint");

pdf(file="history.pdf",paper="letter",height=7.5,width=10);

# total population and proportion population by subpopulation
t.pop.series <- apply(pop.series,c(1,3),sum,na.rm=T);
# proportion in age, by subpopulation
t.pop.prop.series <- sweep(pop.series,c(1,3),t.pop.series,"/");

fyr=1; lyr=300;

par(mfrow=c(1,1));
# TOTPOP.bysub
matplot(fyr:lyr,t.pop.series[fyr:lyr,],type="l",xlab="Cycle",ylab="",col=subPopTypes )
title(paste(sim.num,"Total population size, by subpopulation\n TOTPOP.bysub  "))
legend(x="bottomright",fill=subPopTypes,legend=subPopLabels)
# TOTL10POP.bysub
matplot(fyr:lyr,log10(t.pop.series[fyr:lyr,]),type="l",xlab="Cycle",ylab="log10(Pop)",col=subPopTypes );
title(paste(sim.num,"log10 Total population size, by subpopulation\n  TOTL10POP.bysub "))
legend(x="bottomright",fill=subPopTypes,legend=subPopLabels)


# population in age 0 (proxy for births)
i.ageclass=1;                           # age 0-5

lyr=400;
matplot(fyr:lyr,pop.series[fyr:lyr,i.ageclass,],type="l",xlab="Cycle",ylab="",col=subPopTypes )
title(paste(sim.num,"Total population, ages 0-4, by subpopulation\n TOTPOP.X.bysub"))
legend(x="bottomright",fill=subPopTypes,legend=subPopLabels)


matplot(fyr:lyr,t.pop.prop.series[fyr:lyr,i.ageclass,],type="l",ylim=c(0.05,0.2),xlab="Cycle",ylab="",col=subPopTypes )
title(" Proportion ages 0-4, by subpopulation\n  PPOP.X.bysub")
legend(x="bottomright",fill=subPopTypes,legend=subPopLabels)

## Age Structure
f.plot<-function(tt){
matplot(1:16,(t.pop.prop.series[tt+1,,]),type="l",col=subPopTypes,xlab="Age Class",ylab="")
title(paste(sim.num,"Age Structure,t=",tt,"\n ASPOP.t.bysub "));
legend(x="bottomright",fill=subPopTypes,legend=subPopLabels)

};


par(mfrow=c(2,2))
f.plot(0); f.plot(20); f.plot(150); f.plot(1000)

# Sharing gamma profiles GAMMA.t.bysub
f.plot<-function(tt){
matplot(1:16,(gamma.series[tt+1,,]),type="l",col=subPopTypes,xlab="Age Class",ylab="")
title(paste(sim.num,"Gamma Profile,t=",tt," \n GAMMA.t.bysub"))
legend(x="bottomright",fill=subPopTypes,legend=subPopLabels)
};
par(mfrow=c(2,2))
f.plot(0);abline(h=1.2); f.plot(20);abline(h=1.2);
f.plot(100);abline(h=1.2); f.plot(1000);abline(h=1.2)

# production is captured in yz series  YX.t.bysub
f.plot<-function(tt){
matplot(1:16,(yx.series[tt+1,,]),type="l",col=subPopTypes,xlab="Age Class",ylab="")
title(paste(sim.num,"Production Profile,t=",tt,"\n YX.t.bysub "));
legend(x="bottomright",fill=subPopTypes,legend=subPopLabels)
};
par(mfrow=c(2,2))
f.plot(0); f.plot(5); f.plot(20); f.plot(100)

#matplot(1:16,qx.by.momalive[tt+1, ,1,],type="l",col=subPopTypes,xlab="Age Class",ylab="")

par(mfrow=c(1,1))

lyr=5000
matplot(fyr:lyr,density.series[fyr:lyr,],type="l",xlab="Cycle",col=subPopTypes)
title("Subpopulation Density measure \n(proportional to  total consumption w/in subpopulation) ")

lyr=5000;
par(mfrow=c(2,1))
matplot(fyr:lyr,kinship5.series[fyr:lyr,],type="l",xlab="Cycle",ylab="",col=subPopTypes)
title("Mean kin group (kin5) size within subpopulation ")

matplot(fyr:lyr,sg.series[fyr:lyr,],type="l",xlab="Cycle",ylab="",col=subPopTypes)
title("Mean social group  size within subpopulation ")


par(mfrow=c(1,1));
SkipSmooth<- FALSE;
if(! SkipSmooth){

lyr=5000;  
t.t <- kinship5.series[fyr:lyr,]

c.list<- as.list(NA+subPopTypes);
t.ra <- 999; t.rb <- 0;
#c.range <- NA+subPopTypes;
for(i in seq(subPopTypes)){
  t.i <-  supsmu(fyr:lyr,t.t[,i]);
  c.list[[i]] <- t.i
  #lines(c.list[[i]], col=subPopTypes[i])
  t.ra <- min(t.ra, min(t.i$y));
  t.rb <- max(t.rb, max(t.i$y));
  
}

matplot(fyr:lyr,t.t ,type="n",xlab="Cycle",ylab="",col=subPopTypes,ylim=c(t.ra,t.rb) )
title("Mean kin5 group size, by subpopulation\n bl=M, red=k5,gr=5050");
for(i in seq(subPopTypes)){
  lines(c.list[[i]], col=subPopTypes[i]);
}


## social group sizes
t.t <- sg.series[fyr:lyr,]

c.list<- as.list(NA+subPopTypes);
t.ra <- 999; t.rb <- 0;
#c.range <- NA+subPopTypes;
for(i in seq(subPopTypes)){
  t.i <-  supsmu(fyr:lyr,t.t[,i]);
  c.list[[i]] <- t.i
  #lines(c.list[[i]], col=subPopTypes[i])
  t.ra <- min(t.ra, min(t.i$y));
  t.rb <- max(t.rb, max(t.i$y));
  
}

matplot(fyr:lyr,t.t ,type="n",xlab="Cycle",ylab="",col=subPopTypes,ylim=c(t.ra,t.rb) )
title("Mean social group size, by subpopulation\n bl=M, red=k5,gr=5050");
for(i in seq(subPopTypes)){
  lines(c.list[[i]], col=subPopTypes[i]);
}

}



par(mfrow=c(2,2))
for(i in c(0,2,4,6,8,10,12,14)){
  matplot(fyr:lyr,gamma.series[fyr:lyr,i+1,],type="l",col=subPopTypes,xlab="Cycle")
  title(paste("Age",5*i," Gamma (50% kin5, 50% sg)\n  within subpopulation "))

}

lyr=5000
for(i in c(0,2,4,6,8,10,12,14)){
  matplot(((lyr-100):lyr),gamma.series[((lyr-100):lyr),i+1,],type="l",col=subPopTypes,xlab="Cycle")
  title(paste("Age",5*i," Gamma (50% kin5, 50% sg)\n  within subpopulation "))

}



for(i in c(0,2,4,6,8,10,12,14)){
  t.t <- qx.series[fyr:lyr, i+1,]
  c.list<- as.list(NA+subPopTypes);
  t.ra <- 1; t.rb <- 0;
  
  for(j in seq(c.list)){
    ci <- supsmu(fyr:lyr,qx.series[fyr:lyr,i+1,j])
    t.ra <- min(t.ra, min(ci$y));
    t.rb <- max(t.rb, max(ci$y));
    c.list[[j]] <- ci
  }

  matplot(fyr:lyr,t.t ,type="n",xlab="Cycle",ylab="",col=subPopTypes,ylim=c(t.ra,t.rb) )
  title(paste("qx age class",i,"by subpopulation\n bl=M, red=k5,gr=5050"));
  for(k in seq(subPopTypes)){
    lines(c.list[[k]], col=subPopTypes[k]);
  }
  
}
      



par(mfrow=c(1,1))

lyr=5000;
t.t <- fission.series[fyr:lyr,];

c.list<- as.list(NA+subPopTypes);
t.ra <- 1; t.rb <- 0;
#c.range <- NA+subPopTypes;
for(i in seq(subPopTypes)){
  t.i <-  supsmu(fyr:lyr,t.t[,i]);
  c.list[[i]] <- t.i
  #lines(c.list[[i]], col=subPopTypes[i])
  t.ra <- min(t.ra, min(t.i$y));
  t.rb <- max(t.rb, max(t.i$y));
  
}
matplot(fyr:lyr,t.t ,type="n",xlab="Cycle",ylab="",col=subPopTypes,ylim=c(t.ra,t.rb) )
title("Proportion Fissioning, by subpopulation\n bl=M, red=k5,gr=5050");
for(i in seq(subPopTypes)){
  lines(c.list[[i]], col=subPopTypes[i]);
}



lyr=5000;
t.t <- fusion.series[fyr:lyr,];

c.list<- as.list(NA+subPopTypes);
t.ra <- 1; t.rb <- 0;
#c.range <- NA+subPopTypes;
for(i in seq(subPopTypes)){
  t.i <-  supsmu(fyr:lyr,t.t[,i]);
  c.list[[i]] <- t.i
  #lines(c.list[[i]], col=subPopTypes[i])
  t.ra <- min(t.ra, min(t.i$y));
  t.rb <- max(t.rb, max(t.i$y));
  
}
matplot(fyr:lyr,t.t ,type="n",xlab="Cycle",ylab="",col=subPopTypes,ylim=c(t.ra,.04) )
title("Proportion Fusioning, by subpopulation\n bl=M, red=k5,gr=5050");
for(i in seq(subPopTypes)){
  lines(c.list[[i]], col=subPopTypes[i]);
}
