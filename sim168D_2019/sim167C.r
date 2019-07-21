## sim167. 15 independent groups compete with weak
## density dependence.  There are 3 types of groups, 5 maternal, 5 kin5, and 5 with 50/50 (kin5) sharing

## new density measurement, based on consumption weighted population.

## New mortality variant  h(x) = *K(x)*delta + xi)*gamma^(e.mortality)




options(save.defaults=list(compress=TRUE)); #avoid unnecessarily huge files


start.time <- date()
  
source("forward.5yr.altqx.Tcons.sociality.r");         #main simulation routine
source("lifetable.1qx.r")
#source("f.disperse.r");
source("f.showmode.r");
source("f.kinshipidX.R");                #define CB kinship algorithms.
source("f.Timapply.R");
#############################################################################################
#############################################################################################
# Fertility experiments that use different scalings of the baseline ferfility
# This is controlled with parameter f.adjust.

# 
############################################################################################


#############################
# Control parameters
############################
start.year = 1;                         # launch time
Ncycles =  5000;            # number of iterations in this simulation
Ncheckpoint = 500;                      # write out checkpoints at this interval
#subPopTypes <- c( 1,1,1,1,1,   2,2,2,2,2,   3,3,3,3,3);
subPopTypes <- c(3,3 );

subPopLabels<-   c("M","kin5","50/50")[subPopTypes];

max.age <- 15;                          # 

###############################
# Input Parameters
##############################

sim.num <- c("Sim 167B")
mute.rate <- .01
additive.gene.risk <- .1/5;		# effect of a mutation on mortality
background.mort <- 0.005;               # extrinsic mortality, affected by gamma, not genes;
e.mort <- -1
e.fert <- 1
e.prod <- 0.5
e.dens <- -1

adjust.increment <- 0.1
# beta are sharing mixture weights
beta.kin <- 0.50     ;                        # kinship.id indexed kin shares
beta.mat <- 0      ;                        # older matriline sharing
# group sharing beta is 1-beta.kin - beta.mat
beta.sg1 <-  (1- beta.kin - beta.mat);

i.popsize <- 10000;                      # initial and target population size; resources will adjust round this level
f.adjust  <- 1.0;                       #for fertility experiments, 1=baseline
dispersal.rate <- 0.5;                  # for now, use tim's dispersal
resources <- 69361 * (i.popsize/100000);# resources was set empirically, scale it for smaller popsizes

## algorithm choices
setFusion<-'useTMfusion';               #Fusion of small groups based on odd/even merges; no relatedness preference in merge
setFission<-'useTMpk5fission';           #Fission by splitting large social groups; kin groups are kept together in the splits
setDispersal<- 'useNone'


################################
# Initialize Population
################################

###################
# Initial Mortality
###################
qx.baseline <- rep(0,16)
qx.baseline[16] <- 1.0
px <- 1-qx.baseline
lx <- c(1,cumprod(px)[1:15])

###################
# Initial Fertility
###################
fx <- c(0,0,0.0088,0.1536,0.275,
        0.298,0.318,0.279,0.219,0.0622,rep(0,6))
adjust <- sum(fx*lx[1:16])
fx <- fx/adjust
fx <- fx* f.adjust;

####################
# Initial Population
####################

i.gamma.sg1 <- rep(1,i.popsize)
i.childhood.gamma.sg1 <- rep(1,i.popsize)
i.gamma.mat <- rep(1,i.popsize)
i.childhood.gamma.mat <- rep(1,i.popsize)
i.genes <- matrix(as.integer(0),i.popsize,16)
i.ages <- as.integer( sample(seq(0,15),i.popsize,T,prob=lx[1:16]) )
#i.ages <- factor(i.ages, levels=0:15);

# Caloric Production among 3 Kaplan populations; see M,F combined C,P schedules
# in kapconv.xls
yxz<- # 0:5:90
      c(0.956462585,158.6018677,640.812987,2065.538596,4807.555662,6180.768293,
        6180.768293,6180.768293,7093.225806,7093.225806,7378,7378,5035.5,
        5035.5,3000,3000,1500,1500,0);


#yxz <- yxz*1

# Caloric Consumption among 3 hunter gatherer pops; from Kaplan schedules;
#  see kapconv.xls
cxz<-# 0:5:90
  c(1394.606803,2191.669202,2726.972078,3594.029825,3771.813397,3558.792683,
    3558.792683,3558.792683,3596.451613,3596.451613,3383.083333,3383.083333,
    2710.7, 2710.7,2500,2500,2000,2000,0);




adjust <- sum(yxz[i.ages+1])/sum(cxz[i.ages+1])
i.gamma.sg1 <- rep(adjust,length(i.ages))


#initialize id variables.
# assume one big family to start

i.own.id.mat <- seq (1,i.popsize)
i.mom.id.mat <- i.own.id.mat - 1* i.popsize; 
i.gmom.id.mat <- i.own.id.mat - 2* i.popsize;
i.ggmom.id.mat <- i.own.id.mat - 3* i.popsize;
i.gggmom.id.mat <-i.own.id.mat - 3* i.popsize;
i.group.id.mat <- i.own.id.mat

i.own.id.sg1 <- i.own.id.mat
i.mom.id.sg1 <- as.integer(rep(NA,i.popsize)); # curiously, NAs default to reals
i.gmom.id.sg1 <- as.integer(rep (NA,i.popsize))
i.ggmom.id.sg1 <- as.integer(rep(NA,i.popsize))
i.group.id.sg1 <- i.own.id.sg1


# extra stuff needed for 18c
  i.gggmom.id.mat <- NA+ i.ggmom.id.mat;
  i.kinship.id <-  i.ggmom.id.mat;
  i.own.id.yr <- i.mom.id.yr <-  i.gmom.id.yr <-  i.ggmom.id.yr <-  i.gggmom.id.yr <- NA+ i.ggmom.id.mat;                    
  i.mom.id.kin<- NA+i.own.id.mat;
  i.gmom.id.kin<- NA+ i.gmom.id.mat;
  i.ggmom.id.kin<-NA+ i.ggmom.id.mat;
  i.gggmom.id.kin<- NA+ i.gggmom.id.mat;

## with consumption-weighted population size, adjust 'resources' so
## that the level is comparable with this new measure

#resources <- 69361 * (i.popsize/100000);#  population based resource size
t.TotC <- sum( cxz[i.ages + 1]);           #total consumption across all population
t.adjust <- t.TotC/i.popsize;
Cresources <- t.adjust * resources;


# crudely clump starting population into groups of approx target size to avoid
# very slow startup; this is not necessary if starting from a prior run

i.group.id.sg1 <- as.integer( 1+ floor( i.own.id.sg1 / 8) );

# group population into almost equal starting pops, keeping starting clumps together
t.nclumps <- length(unique(i.group.id.sg1)); # number of clumps to reassign to subPop types
t.isubPop <- 1+ (seq(t.nclumps) %% length(subPopTypes) ); #split up clumps into subPop groupings
i.subPop.id <- t.isubPop[match( i.group.id.sg1, seq(t.isubPop) ) ];

#i.subPop.id <- factor(subPopTypes[i.subPop.id],levels=subPopTypes);


###################################
### Forecast 1 cycle forward in time
###################################


test<-
  forward.5yr.altqx(cycleno=as.integer(0),
                     beta.mat=beta.mat,
                     beta.kin=beta.kin,
                     beta.sg1=beta.sg1,                     
                     own.id.mat=i.own.id.mat,
                     mom.id.mat= NA+as.integer(i.own.id.mat - 1*i.popsize),
                     gmom.id.mat=NA+as.integer(i.own.id.mat - 2*i.popsize),
                     ggmom.id.mat=NA+as.integer(i.own.id.mat - 3*i.popsize),
                     mom.id.kin = as.integer(i.own.id.mat - 1*i.popsize),
                     gmom.id.kin= as.integer(i.own.id.mat - 2*i.popsize),
                     ggmom.id.kin=as.integer(i.own.id.mat - 3*i.popsize),
                     gggmom.id.kin=as.integer(i.own.id.mat - 4*i.popsize),
                     ggggmom.id.kin=as.integer(i.own.id.mat - 5*i.popsize),
                     gggggmom.id.kin= as.integer(i.own.id.mat - 6*i.popsize),                    
                     kinship3.id = i.group.id.sg1, #should not matter
                     kinship4.id = i.group.id.sg1,
                     kinship5.id = i.group.id.sg1,                   
                     group.id.anc= as.integer(i.own.id.mat - 4*i.popsize) , #Tim's ancestory grouping variable
                     own.id.yr= as.integer(0*i.own.id.mat),
                     mom.id.yr= NA+  i.mom.id.mat,
                     gmom.id.yr =  NA+i.gmom.id.mat,
                     ggmom.id.yr = NA+ i.ggmom.id.mat,
                     gggmom.id.yr= NA+ i.gggmom.id.mat,
                     group.id.mat= i.group.id.mat,
                     gamma.mat= i.gamma.mat,
                     childhood.gamma.mat = i.childhood.gamma.mat,
                     gamma.kin = i.gamma.mat,
                     childhood.gamma.kin = i.childhood.gamma.mat,
                     gamma.anc = i.gamma.mat, #new
                     childhood.gamma.anc = i.childhood.gamma.mat, #new
                     ##own.id.sg1 = i.own.id.sg1,  # DELETE
                     mom.id.sg1 = i.mom.id.sg1,
                     gmom.id.sg1 = i.gmom.id.sg1,
                     ggmom.id.sg1 = i.ggmom.id.sg1,
                     group.id.sg1 = i.group.id.sg1,
                     gamma.sg1 = i.gamma.sg1,
                     childhood.gamma.sg1 = i.childhood.gamma.sg1,
                     subPop.id = i.subPop.id,       #new vector to track subpopulations
                     gene = i.genes,
                     age = i.ages,      
                     max.age = as.integer(max.age) ,
                     e.fertility = e.fert,
                     e.mortality =  e.mort,
                     e.production =  e.prod,
                     e.density =  e.dens,
                     #pop.size= length(i.own.id.mat),
                     resource.size = Cresources
                     );




rm(i.own.id.mat,i.mom.id.mat, i.gmom.id.mat, i.ggmom.id.mat,
   i.gggmom.id.mat, i.kinship.id, i.group.id.mat, i.gamma.mat,
   i.gggmom.id.kin,i.gggmom.id.yr,i.ggmom.id.kin,
   i.ggmom.id.yr,i.gmom.id.kin ,i.gmom.id.yr,
   i.mom.id.kin , i.mom.id.yr ,i.own.id.yr,
   i.childhood.gamma.mat, i.own.id.sg1, i.mom.id.sg1, i.gmom.id.sg1,
   i.ggmom.id.sg1, i.group.id.sg1, i.gamma.sg1, i.childhood.gamma.sg1,
   i.genes, i.ages)


###################################
## Reinitialize test using previous simulation, if appropriate
###################################
#load("sim107.test.Rbin");               #reloads test, final state from sim107

## reset mutation load to zero
#test$genes <- ( 0* test$genes);

## convert to integer to save memory and increase speed
#test<- f.setmode(test);

########################
## Forecast Ncycles years
########################
yx.series <- array(NA, dim=c(Ncycles,16,length(subPopTypes) ) );
sim.results <-  array(NA, dim=c(Ncycles,16,length(subPopTypes) ) );
pop.series <-  array(NA, dim=c(Ncycles,16,length(subPopTypes) ) );


gamma.series <- array(NA, dim=c(Ncycles,16,length(subPopTypes) ) );
density.series <- matrix(NA,Ncycles,length(subPopTypes));
fusion.series <-  matrix(NA,Ncycles,length(subPopTypes));
fission.series <- matrix(NA,Ncycles,length(subPopTypes));
mat.series <-  matrix(NA,Ncycles,length(subPopTypes));
anc.series <-  matrix(NA,Ncycles,length(subPopTypes));
kinship3.series <-  matrix(NA,Ncycles,length(subPopTypes));
kinship4.series <-  matrix(NA,Ncycles,length(subPopTypes));
kinship5.series <-  matrix(NA,Ncycles,length(subPopTypes));
sg.series       <-  matrix(NA,Ncycles,length(subPopTypes));

qx.series <-  NA+array(0, dim=c(Ncycles,16,length(subPopTypes) ) );

# 100K size test
too.high <- ifelse(length(test$own.id.mat)> i.popsize,"yes","no");

# homeostasis correction flags
do.homeo.resource <- FALSE;                #toggle resource homeostasis
do.2Xpopsize.check     <- TRUE ;                #toggle 2Xpopsize  chop down

# interval at which to apply homeostatic checks
homeo.N <- 100;

# window to average over to determine if pop is growing or shrinking
pop.avgwindow <- 20;

set.seed(999);

for (year in seq(Ncycles)){
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
                     ##own.id.sg1 = test$own.id.sg1,  # DELETE
                     mom.id.sg1 = test$mom.id.sg1,
                     gmom.id.sg1 = test$gmom.id.sg1,
                     ggmom.id.sg1 = test$ggmom.id.sg1,
                     group.id.sg1 = test$group.id.sg1,
                     gamma.sg1 = test$gamma.sg1,
                     childhood.gamma.sg1 = test$childhood.gamma.sg1,
                     subPop.id = test$subPop.id,
                     gene = test$gene,
                     age = test$age,
                     max.age = max.age,
                     e.fertility = e.fert,
                     e.mortality =  e.mort,
                     e.production =  e.prod,
                     e.density =  e.dens,
                     #pop.size= length(test$own.id.mat),
                     resource.size = Cresources
                     );


 
  #sim.results[year,] <- test$mean.harm
  gamma.series[year,,] <-  tapply(test$sharing.gamma,list(factor(test$age,levels=0:15),
                                  factor(test$subPop.id,levels=seq(subPopTypes) )),mean);
  t.age.factor <- factor(test$age,levels=0:15);
  t.subpop.factor <- factor(test$subPop.id,levels=seq(subPopTypes));
for(i in seq(subPopTypes)){
  # avg mutation load, by age class, for each subpop
  t.subpop<- test$subPop.id==i;
  if(any(t.subpop)){                    #if there is a 0 subpop, output remains NA
    sim.results[year,,i] <-apply(test$gene[t.subpop,, drop=F],2,mean);
     anc.series[year,i] <-      mean( table(test$group.id.anc[t.subpop]));
    kinship3.series[year,i] <- mean(table(test$kinship3.id[t.subpop]));
    kinship4.series[year,i] <- mean(table(test$kinship4.id[t.subpop]));
    kinship5.series[year,i] <- mean(table(test$kinship5.id[t.subpop]));
    mat.series[year,i]    <-   mean(table(test$group.id.mat[t.subpop]));
    sg.series[year,i]      <-  mean(table(test$group.id.sg1[t.subpop])) ;
  }
}

  pop.series[year,,] <- tapply(test$own.id.mat,list(t.age.factor,t.subpop.factor),length)
  yx.series[year,,]  <-  tapply(test$yx,list(t.age.factor,t.subpop.factor),mean);

  t.totc <- tapply( cxz[ test$age + 1 ], t.subpop.factor, sum );
  density.series[year,] <- t.totc/Cresources
  fusion.series[year,] <- test$Pfusions;
  fission.series[year,] <- test$Pfissions;
 
  qx.series[year,,] <-  test$qx.by.age;

 print(paste("cycle ",year));
  
 if( (year %% Ncheckpoint) == 0){
    print(paste(date(), "CHECKPOINT:  Completed ", start.year + year," of ", start.year + Ncycles));
    # save the results to this point
    save.image(file=".CheckPointC",compress=TRUE);
  }

  # every homeo.N years, adjust the resource.target, if adj. enabled
  if ( do.homeo.resource && (year %% homeo.N == 0) ) {
    print("...Resource adjustment...");
    # large and growing pop -- then shrink it.
    t.totP <- apply(pop.series[(year-pop.avgwindow+1):year,],1,sum);
    if ( (t.totP[pop.avgwindow] > i.popsize) && ( mean(exp(diff(log(t.totP)))) > 1 ) ) {
       #Pop is above target size and is growing, so decrease available resources
      if (too.high=="no") {adjust.increment <- ifelse( adjust.increment>.001, adjust.increment*0.9, .001) }
      too.high <- c("yes")
      resources <- resources - (adjust.increment*resources);
      Cresources <- Cresources - (adjust.increment*Cresources);
    }
    # small and shrinking pop -- then expand it
    if ( (t.totP[pop.avgwindow] < i.popsize)  && ( mean(exp(diff(log(t.totP)))) < 1 )   ) {
                                        #Pop is below target size and is declining, so increase resources
      if (too.high=="yes") {adjust.increment <- ifelse(adjust.increment>.001, adjust.increment*0.9, .001) }
      too.high <- c("no")
      resources <- resources + (adjust.increment*resources);
      Cresources <- Cresources + (adjust.increment*Cresources);
     
    }
  }

  density.series[year,] <-  t.totc/Cresources

#if exceeds twice the normative population size, chop it down
  if (do.2Xpopsize.check && length(test$own.id.mat)> 2*i.popsize){
   # Select i.popsize folks at random
    print(".....Popsize forced reduction....");
    test<- f.subsamptest(test,i.popsize); 
  }

}

end.time <- date()
#rbind(test$own.id.mat,test$mom.id.mat,test$gmom.id.mat,test$ggmom.id.mat, test$gggmom.id.mat,test$kinship.id)[,1:10]


