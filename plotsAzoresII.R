## pulled from history.plots.r

# original plot for Fig 3a
#Figure 3. Other outcomes of the simulations by social arrangement: consumption, production, density,
#dependency ratio, and fertility.

## sim168B
#setwd("/home/boe/Dropbox/LeeBoe2018/sim168B")
setwd("~boe/Dropbox/LeeBoe2018")
setwd("./sim168B")
## investigation of the 1st 1000 cycles of sim168
load(".CheckPoint");


f.plot<-function(tt){
  matplot(1:16,(gamma.series[tt+1,,]),type="l",col=subPopTypes,xlab="Age Class",ylab="")
  title(paste(sim.num,"Gamma Profile,t=",tt," \n GAMMA.t.bysub"))
  legend(x="bottomright",fill=subPopTypes,legend=subPopLabels)
};
#par(mfrow=c(1,1))
#  f.plot(100);abline(h=1.2)  # original plot


## modified plot
library(devEMF)
emf(file="Fig3a.emf")
#png(file="Fig3a.png", width=7, height=7, units="in", res=300)
og<- par(cex.lab=1.5, cex.axis=1.25)
age.labels<- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", 
               "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", 
               "70-74", "75-79")

tt=100

matplot(1:16,(gamma.series[tt+1,,]),type="l",col=subPopTypes,xlab="age group",ylab="consumption ratio",
        axes=FALSE)
odd.ticks <- 1:16 %% 2 == 1
axis(1,at=(1:16)[ odd.ticks ], labels=age.labels[odd.ticks])
axis(1, at=(1:16), tick=TRUE, labels=FALSE)

 # replace by labels right next to each set of three lines: M.100, K5.100, K5.50.
text(10, 1.57,"M.100", cex=1.2)
text(12.49,1.36,"K5.100", cex=1.2)
text(12.8, 1.18,"SG.K5.50", cex=1.2)
axis(2)
dev.off()


## Fig3b 

# production is captured in yz series
f.plot<-function(tt){
  matplot(1:16,(yx.series[tt+1,,]),type="l",col=subPopTypes,xlab="Age Class",ylab="")
  title(paste(sim.num,"Production Profile,t=",tt,"\n YX.t.bysub "));
  legend(x="bottomright",fill=subPopTypes,legend=subPopLabels)
};

# f.plot(100) # original plot

emf(file="Fig3b.emf")
#png(filename="Fig3b.png", width=7, height=7, units="in", res=300)
og<- par(cex.lab=1.5, cex.axis=1.25)
matplot(1:16,(yx.series[tt+1,,]),type="l",col=subPopTypes,xlab="age group",ylab="production (in Kcal)",
        axes=FALSE)
odd.ticks <- 1:16 %% 2 == 1
axis(1,at=(1:16)[ odd.ticks ], labels=age.labels[odd.ticks])
axis(1, at=(1:16), tick=TRUE, labels=FALSE)
text(6.13, 6851.240,"M.100", cex=1.2)
text(7.14,5614.966,"K5.100", cex=1.2)
text(8.52, 3974,"SG.K5.50", cex=1.2)
axis(2)
dev.off()

## Fig4a is from sim167B  (Fig6a in new document)
setwd("../sim167B")
## investigation of the 1st 1000 cycles of sim168
load(".CheckPoint");
# total population and proportion population by subpopulation
t.pop.series <- apply(pop.series,c(1,3),sum,na.rm=T);
# proportion in age, by subpopulation
t.pop.prop.series <- sweep(pop.series,c(1,3),t.pop.series,"/");
# TOTPOP.bysub -- original plot
# matplot(fyr:lyr,t.pop.series[fyr:lyr,],type="l",xlab="Cycle",ylab="",col=subPopTypes )
# title(paste(sim.num,"Total population size, by subpopulation\n TOTPOP.bysub  "))
# legend(x="bottomright",fill=subPopTypes,legend=subPopLabels)

# log(Pop).bysub
fyr=1
lyr=1000
emf(file="Fig4a.emf")
og<- par(cex.lab=1.5, cex.axis=1.25)
matplot(fyr:lyr,log(t.pop.series[fyr:lyr,]),type="l",xlab="Cycle",ylab="log(Pop)",col=subPopTypes, ylim=c(4,10.5) );
#title(paste(sim.num,"log Total population size, by subpopulation\n  TOTL10POP.bysub "))
#legend(x="bottomright",fill=subPopTypes,legend=subPopLabels)
text(250,9.96,"M.100", cex=1.2)
text(60, 4.92,"SG.K5.50", cex=1.2)
text(300,5.33,"K5.100", cex=1.2)

dev.off()

## Fig4b is from sim168D not sim168C (sim.num was not set correctly) [Fig6b in revised document]
## here childhood sharing contribution is erased, using e.cprod=0 instead of e.cprod=e.prod=0.5

setwd("../sim168D")
## investigation of the 1st 1000 cycles of sim168
load(".CheckPoint");
# total population and proportion population by subpopulation
t.pop.series <- apply(pop.series,c(1,3),sum,na.rm=T);
fyr=1
lyr=3000
emf(file="Fig4b.emf")
og<- par(cex.lab=1.5, cex.axis=1.25)
matplot(fyr:lyr,log(t.pop.series[fyr:lyr,]),type="l",xlab="Cycle",ylab="log(Pop)",col=subPopTypes, ylim=c(4,10.5) );
text(550, 6.46,"M.100", cex=1.2)
text(1560, 7.95,"K5.100", cex=1.2)
text(1846, 9.66,"SG.K5.50", cex=1.2)
dev.off()


## alternative tests with different caps on the effect of childhood sharing gamma

## alternative, capping e.cprod effect at 1.5
## Fig4b is from sim168D not sim168C (sim.num was not set correctly)
setwd("../sim168D_2019")
## investigation of the 1st 1000 cycles of sim168
load(".CheckPoint");
# total population and proportion population by subpopulation
t.pop.series <- apply(pop.series,c(1,3),sum,na.rm=T);
fyr=1
lyr=500
emf(file="Fig4b_scg.cap1.5.emf")
matplot(fyr:lyr,log(t.pop.series[fyr:lyr,]),type="l",xlab="Cycle",ylab="log(Pop)",col=subPopTypes, ylim=c(4,10.5) );
text(366.1062, 9.595500,"M.100")
text(374., 7.513230,"K5.100")
text(296., 4.802351,"SG.K5.50")
title("Sim 168D_2019 Log(Population) with e.cprod = e.prod = 0.5, \n with 1.5 cap on scg")
dev.off()

## boxplot showing distribution of scg in the final population (N=Ncycles=500)
t.labels <- c("M.100\nSim1", "M.100\nSim2", "M.100\nSim3",
              "K5.100\nSim1", "K5.100\nSim2", "K5.100\nSim3",
              "K5.50.50\nSim1", "K5.50.50\nSim2", "K5.50.50\nSim3")

emf(file="newFig4_scg.cap1.5.emf")
#pdf(file="newFig4_scg.cap1.5.pdf")
op<- par(mar=c(5+2,4+2,4,2)+0.1, cex=1.1)
boxplot(test$sharing.childhood.gamma ~ test$subPop.id, axes=FALSE)
axis(2)
axis(1, labels=FALSE)
box()
# multiplier on scg is min( 1.5, (sharing.childhood.gamma)^0.5 )  so 1.5^2 is the cut point
# above which the multiplier is not increasing
#abline(h= 1.5^2, col="red")
mtext("Consumption ratio", side=2, line=3, cex=1.5)
mtext(t.labels[1:6], side=1, at=1:6, line=3)
par(op)
dev.off()
## alternative, capping e.cprod effect at 1.25
## Fig4b is from sim168D not sim168C (sim.num was not set correctly)
setwd("../sim168D_2019_cap1.25")
## investigation of the 1st 1000 cycles of sim168
load(".CheckPoint");
# total population and proportion population by subpopulation
t.pop.series <- apply(pop.series,c(1,3),sum,na.rm=T);
fyr=1
lyr=500
emf(file="Fig4b_scg.cap1.25.emf")
matplot(fyr:lyr,log(t.pop.series[fyr:lyr,]),type="l",xlab="Cycle",ylab="log(Pop)",col=subPopTypes, ylim=c(4,10.5) );
text(313.9978, 8.139354,"M.100")
text(384.7164, 9.083422,"K5.100")
text(453.9461, 10.081130,"SG.K5.50")
title("Sim 168D_2019 Log(Population) with e.cprod = e.prod = 0.5, \n with 1.25 cap on scg")
dev.off()



## alternative, capping e.cprod effect at 1.35
## Fig4b is from sim168D not sim168C (sim.num was not set correctly)
setwd("../sim168D_2019_cap1.35")
## investigation of the 1st 1000 cycles of sim168
load(".CheckPoint");
# total population and proportion population by subpopulation
t.pop.series <- apply(pop.series,c(1,3),sum,na.rm=T);
fyr=1
lyr=500
emf(file="Fig4b_scg.cap1.35.emf")
matplot(fyr:lyr,log(t.pop.series[fyr:lyr,]),type="l",xlab="Cycle",ylab="log(Pop)",col=subPopTypes, ylim=c(4,10.5) );
text(430.1251, 10.252778,"M.100")
text(401., 8.729397,"K5.100")
text(200.8481, 7.731689,"SG.K5.50")
title("Sim 168D_2019_cap1.35 Log(Population) with e.cprod = e.prod = 0.5, \n with 1.35 cap on scg")
dev.off()



## alternative, capping e.cprod effect at 1.30
## Fig4b is from sim168D not sim168C (sim.num was not set correctly)
setwd("../sim168D_2019_cap1.30")
## investigation of the 1st 1000 cycles of sim168
load(".CheckPoint");
# total population and proportion population by subpopulation
t.pop.series <- apply(pop.series,c(1,3),sum,na.rm=T);
fyr=1
lyr=1500
#emf(file="Fig4b_scg.cap1.30.emf")
pdf(file="Fig4b_scg.cap1.30.pdf")
png(file="Fig4b_scg.cap1.30.png")

matplot(fyr:lyr,log(t.pop.series[fyr:lyr,]),type="l",xlab="Cycle",ylab="log(Pop)",col=subPopTypes, ylim=c(8,10.5) );
text(458.56, 9.19,"M.100")
text(282.85, 9.89,"K5.100")
text(434.56, 8.41,"SG.K5.50")
title("Sim 168D_2019_cap1.30 Log(Population) with e.cprod = e.prod = 0.5, \n with 1.30 cap on scg")
dev.off()
