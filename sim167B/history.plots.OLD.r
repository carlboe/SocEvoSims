

## investigation of the 1st 1000 cycles of sim16x
load(".CheckPoint");

pdf(file="history.pdf",paper="letter",height=7.5,width=10);

# total population and proportion population by subpopulation
t.pop.series <- apply(pop.series,c(1,3),sum,na.rm=T);
# proportion in age, by subpopulation
t.pop.prop.series <- sweep(pop.series,c(1,3),t.pop.series,"/");

fyr=1; lyr=1000;

matplot(fyr:lyr,t.pop.series[fyr:lyr,],type="l",xlab="Cycle",ylab="",col=subPopTypes )
title(paste(sim.num,"Total population size, by subpopulation\n bl=M, red=k5,gr=5050"))


matplot(fyr:lyr,log(t.pop.series[fyr:lyr,]),type="l",xlab="Cycle",ylab="log(Pop)",col=subPopTypes );
title(paste(sim.num,"log Total population size, by subpopulation\n bl=M, red=k5,gr=5050"))

# population in age 0 (proxy for births)
i.ageclass=1;                           # age 0-5

lyr=400;
matplot(fyr:lyr,pop.series[fyr:lyr,i.ageclass,],type="l",xlab="Cycle",ylab="",col=subPopTypes )
title(paste(sim.num,"Total population, ages 0-4, by subpopulation\n bl=M, red=k5,gr=5050"))


matplot(fyr:lyr,t.pop.prop.series[fyr:lyr,i.ageclass,],lty=1,type="l",xlab="Cycle",ylab="",col=subPopTypes )
title(" Proportion ages 0-4, by subpopulation\n bl=M, red=k5,gr=5050")

# Age Structure
f.plot<-function(tt){
matplot(1:16,(t.pop.prop.series[tt+1,,]),type="l",col=subPopTypes,xlab="Age Class",ylab="")
title(paste(sim.num,"Age Structure,t=",tt," by subpopulation\n bl=M, red=k5,gr=5050"))
};

par(mfrow=c(2,2))
f.plot(0); f.plot(50); f.plot(200); f.plot(1000)

# Sharing gamma profiles
f.plot<-function(tt,mycol=subPopTypes){
matplot(1:16,(gamma.series[tt+1,,]),type="l",col=mycol,lty=1,xlab="Age Class",ylab="")
title(paste(sim.num,"Gamma Profile,t=",tt," by subpopulation\n bl=M, red=k5,gr=5050"))
};
par(mfrow=c(2,2))
f.plot(0); f.plot(50); f.plot(200); f.plot(1000)

t.t <- t.pop0.series[fyr:lyr,]/t.pop.series[fyr:lyr,];


matplot(fyr:lyr,t.t ,type="n",xlab="Cycle",ylab="",col=subPopTypes,ylim=c(t.ra,t.rb) )
title("Proportion ages 25-30, by subpopulation\n bl=M, red=k5,gr=5050");
for(i in seq(subPopTypes)){
  lines(c.list[[i]], col=subPopTypes[i]);
}

t.ageclasses=0:15;
par(mfrow=c(2,2))
for(iyr in (0:7)*1000){
  t.t<- pop.series[iyr+1,,];
  t.sum<-apply(t.t,2,sum);
  t.t<- sweep(t.t,2,t.sum,"/");           #proportions within subpops
  matplot(t.ageclasses,t.t,col=subPopTypes,type="l",xlab="Cycle");
  title(paste("Age structure, cycle=",iyr,"\n bl=M, red=k5,gr=5050"));
}  


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

lyr=11000;  
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

lyr=15000
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

lyr=9000;
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



lyr=9000;
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
