## candidate plots for fig3 of Ron's NRC submission
## 'Intergenerational Transfers, Life Histories, and the Elderly'



## investigation of the 1st 1000 cycles of sim168
load(".CheckPoint");

pdf(file="NRC_A.pdf", colormodel="gray");

# total population and proportion population by subpopulation
t.pop.series <- apply(pop.series,c(1,3),sum);
# proportion in age, by subpopulation
t.pop.prop.series <- sweep(pop.series,c(1,3),t.pop.series,"/");

fyr=1; lyr=600;

par(mfrow=c(1,1));
# TOTPOP.bysub
matplot(fyr:lyr,t.pop.series[fyr:lyr,],type="l",xlab="Cycle",ylab="Population",col=subPopTypes, lwd=1.8 )
#title(paste(sim.num,"Total population size, by subpopulation\n TOTPOP.bysub  "))
#legend(x="bottomright",fill=subPopTypes,legend=subPopLabels)

text(210,3517.67,label="Matriarchy")
text(460,10200, label="3rd cousin")
text(265,17000,label="Compound 8-25")
dev.off()
embedFonts("NRC_A.pdf")

matplot(fyr:lyr,log(t.pop.series[fyr:lyr,]),type="l",xlab="Cycle",ylab="log(Pop)",col=subPopTypes );
title(paste(sim.num,"log Total population size, by subpopulation\n  TOTLPOP.bysub "))
legend(x="bottomright",fill=subPopTypes,legend=subPopLabels)
