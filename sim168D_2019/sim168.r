## sim168D_2019.  Redo of sim168B with childhood sharing gamma effect limited to a factor of 1.5

## to see if the winning strategy changes

## 9 independent groups compete with weak
## density dependence.  There are 3 types of groups, 3 maternal, 3 kin5, and 3 with 50/50 (kin5) sharing


## new density measurement, based on consumption weighted population.

## New mortality variant  h(x) = *K(x)*delta + xi)*gamma^(e.mortality)

options(save.defaults=list(compress=TRUE)); #avoid unnecessarily huge files
library(compiler)

## load in and save the results of the individual type simulations
load(".CheckPointA");                   #Maternal type
testA <- test;
load(".CheckPointB");                   #two kin5 type
testB <- test;
load(".CheckPointC");                   #two 50/50 types
testC <- test;


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
Ncycles =  500;            # number of iterations in this simulation
Ncheckpoint = 200;                      # write out checkpoints at this interval
subPopTypes <- c( 1,1,1,  2,2,2,  3,3,3);

subPopLabels<-   c("M","kin5","50/50")[subPopTypes];

max.age <- 15;                          # 

###############################
# Input Parameters
##############################

sim.num <- c("Sim 168D_2019")
mute.rate <- .01
additive.gene.risk <- .1/5;		# effect of a mutation on mortality
background.mort <- 0.005;               # extrinsic mortality, affected by gamma, not genes;
e.mort <- -1
e.fert <- 1
e.prod <- 0.5;
e.cprod <- 0.5;                           # e.cprod = e.prod normally,
cprod.cap = 1.5                       #
e.dens <- -1

adjust.increment <- 0.1
# beta are sharing mixture weights
beta.kin <- 0.50     ;                        # kinship.id indexed kin shares
beta.mat <- 0      ;                        # older matriline sharing
# group sharing beta is 1-beta.kin - beta.mat
beta.sg1 <-  (1- beta.kin - beta.mat);

i.popsize <- 100000;                      # initial and target population size; resources will adjust round this level
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
f.setinitpop <- function(cname, convert=TRUE){
  # cname is the component name
  res<- NULL;
  t.A <- testA[[cname]];
  t.B <- testB[[cname]];
  t.C <- testC[[cname]];
  
  if(is.null(dim(t.A))){    
    res <- c(t.A,t.A,t.A,  t.B,t.B,t.B,  t.C,t.C,t.C);
  } else {
    
    if(cname=="gene"){
      res <- rbind(t.A,t.A,t.A,  t.B,t.B,t.B,  t.C,t.C,t.C);
    } else {
      error("Houston...we have a problem")
    }
  }
  if(convert)
    storage.mode(res) <- "integer";
  
  return(res)
}
  
i.gamma.sg1 <-f.setinitpop("gamma.sg1",F);
i.childhood.gamma.sg1 <- f.setinitpop("childhood.gamma.sg1",F);
i.gamma.mat <-f.setinitpop("gamma.mat",F);
i.childhood.gamma.mat <-f.setinitpop("childhood.gamma.mat",F);

i.gamma.anc <-f.setinitpop("gamma.anc",F);
i.childhood.gamma.anc <-f.setinitpop("childhood.gamma.anc",F);

i.gamma.kin <-f.setinitpop("gamma.kin",F);
i.childhood.gamma.kin <-f.setinitpop("childhood.gamma.kin",F);


i.genes <-f.setinitpop("gene");
i.ages <-f.setinitpop("age");

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




#adjust <- sum(yxz[i.ages+1])/sum(cxz[i.ages+1])
#i.gamma.sg1 <- rep(adjust,length(i.ages))


#initialize id variables.
# all the individual pops have ids less than 5000000, so we build up unique ids

f.setinitids <- function(cname,convert=TRUE){
  # cname is the component name
  res<- NULL;
  t.A <- testA[[cname]];
  t.B <- testB[[cname]];
  t.C <- testC[[cname]];

  t.s <- 5000000;
  t.d <- 4000000;

  # make sure the separate groups share no ids
  res <- c(t.A + t.s,   t.A + t.s + t.d,    t.A + t.s + 2*t.d,
           t.B + 3*t.s, t.B + t.s + 4*t.d,  t.B + t.s + 5*t.d,
           t.C + 6*t.s, t.C + t.s + 7*t.d,  t.C + t.s + 8*t.d);

  if(convert)
    storage.mode(res)<- "integer";
 
  return(res)
}
  
i.own.id.mat <- f.setinitids("own.id.mat");
i.mom.id.mat <- f.setinitids("mom.id.mat");
i.gmom.id.mat <- f.setinitids("gmom.id.mat");
i.ggmom.id.mat <- f.setinitids("ggmom.id.mat");
i.gggmom.id.mat <-f.setinitids("gggmom.id.mat");
                  
i.group.id.mat <- f.setinitids("group.id.mat");

i.own.id.sg1 <- f.setinitids("own.id.sg1");
i.mom.id.sg1 <- f.setinitids("mom.id.sg1");
i.gmom.id.sg1 <- f.setinitids("gmom.id.sg1");
i.ggmom.id.sg1 <-  f.setinitids("ggmom.id.sg1");
i.group.id.sg1 <-  f.setinitids("group.id.sg1");

i.group.id.anc <-  f.setinitids("group.id.anc");

# extra stuff needed for 18c

  zero <- as.integer( 0 );
  i.kinship3.id <-   f.setinitids("kinship3.id");
  i.kinship4.id <-   f.setinitids("kinship4.id");
  i.kinship5.id <-   f.setinitids("kinship5.id");

  i.own.id.yr <-  zero * f.setinitids("own.id.yr");
  i.mom.id.yr <- zero * f.setinitids("mom.id.yr");

  i.gmom.id.yr <- zero * f.setinitids("gmom.id.yr");
  i.ggmom.id.yr <-zero * f.setinitids("ggmom.id.yr");
  i.gggmom.id.yr <- zero * f.setinitids("gggmom.id.yr");

  i.mom.id.kin <-   f.setinitids("mom.id.kin");
  i.gmom.id.kin <-  f.setinitids("gmom.id.kin");

  i.ggmom.id.kin <- f.setinitids("ggmom.id.kin");

  i.gggmom.id.kin<- f.setinitids("gggmom.id.kin");
  i.ggggmom.id.kin<- f.setinitids("ggggmom.id.kin");
  i.gggggmom.id.kin<- f.setinitids("gggggmom.id.kin");


  i.group.id.sg1 <- f.setinitids("group.id.sg1");


## with consumption-weighted population size, adjust 'resources' so
## that the level is comparable with this new measure

i.popsize <- length(i.own.id.mat);      #true init popsize
resources <- 69361 * (i.popsize/100000);#  population based resource size
t.TotC <- sum( cxz[i.ages + 1]);           #total consumption across all population
t.adjust <- t.TotC/i.popsize;
Cresources <- t.adjust * resources;


# crudely clump starting population into groups of approx target size to avoid
# very slow startup; this is not necessary if starting from a prior run

#i.group.id.sg1 <- as.integer( 1+ floor( i.own.id.sg1 / 8) );

# group population into almost equal starting pops, keeping starting clumps together
#t.nclumps <- length(unique(i.group.id.sg1)); # number of clumps to reassign to subPop types
#t.isubPop <- 1+ (seq(t.nclumps) %% length(subPopTypes) ); #split up clumps into subPop groupings
#i.subPop.id <- t.isubPop[match( i.group.id.sg1, seq(t.isubPop) ) ];

#i.subPop.id <- factor(subPopTypes[i.subPop.id],levels=subPopTypes);

i.subPop.id <- c( 1+0*testA$subPop.id, 2+ 0*testA$subPop.id, 3 + 0*testA$subPop.id,
                  4+0*testB$subPop.id, 5+ 0*testB$subPop.id, 6 + 0*testB$subPop.id,
                  7 +0*testC$subPop.id, 8 + 0*testC$subPop.id, 9 + 0*testC$subPop.id );

i.subPop.id <- as.integer(i.subPop.id);

# compute initial shared gammas
i.sharing.gamma <- NA+ i.gamma.mat;
i.sharing.childhood.gamma <- NA + i.gamma.mat;

t.subPops <- seq(subPopTypes)
for( iSub in (t.subPops) ){
  isel <- (t.subPops[iSub] == i.subPop.id ) ;
  t.thistype <- subPopTypes[iSub]; # regime type of this subpop
  if(t.thistype == 1){                    # Matriarchal
    t.beta.sg1 = 0; t.beta.kin=0; t.beta.mat = 1;
  }
  if(t.thistype == 2){                    # 100% kinN
    t.beta.sg1 = 0; t.beta.kin=1; t.beta.mat = 0;
  }
  if(t.thistype == 3){                    # mix of social and kin, no maternal
    t.beta.sg1 = .5; t.beta.kin=.5; t.beta.mat = 0;
  }

  i.sharing.gamma[isel] <- ( t.beta.sg1 * i.gamma.sg1[isel] ) +
    (t.beta.kin * i.gamma.kin[isel] ) + ( t.beta.mat* i.gamma.mat[isel] );
  
  i.sharing.childhood.gamma[isel] <- ( t.beta.sg1 * i.childhood.gamma.sg1[isel] ) +
    (t.beta.kin * i.childhood.gamma.kin[isel]) + ( t.beta.mat * i.childhood.gamma.mat[isel] );
};
 
      
###################################
### Forecast 1 cycle forward in time
###################################


test<-
  forward.5yr.altqx(cycleno=as.integer(0),
                     beta.mat=beta.mat,
                     beta.kin=beta.kin,
                     beta.sg1=beta.sg1,                     
                     own.id.mat=i.own.id.mat,
                     mom.id.mat= i.mom.id.mat,
                     gmom.id.mat= i.gmom.id.mat,
                     ggmom.id.mat= i.ggmom.id.mat,
                     mom.id.kin = i.mom.id.kin,
                     gmom.id.kin= i.gmom.id.kin,
                     ggmom.id.kin= i.ggmom.id.kin,
                     gggmom.id.kin= i.gggmom.id.kin,
                     ggggmom.id.kin= i.ggggmom.id.kin,
                     gggggmom.id.kin=  i.gggggmom.id.kin,
                     kinship3.id = i.kinship3.id,
                     kinship4.id = i.kinship4.id,
                     kinship5.id = i.kinship5.id,
                     group.id.anc= i.group.id.anc,
                     own.id.yr= i.own.id.yr,
                     mom.id.yr=  i.mom.id.yr,
                     gmom.id.yr =  i.gmom.id.yr,
                     ggmom.id.yr = i.ggmom.id.yr,
                     gggmom.id.yr= i.gggmom.id.yr,
                     group.id.mat= i.group.id.mat,
                     gamma.mat= i.gamma.mat,
                     childhood.gamma.mat = i.childhood.gamma.mat,
                     gamma.kin = i.gamma.kin,
                     childhood.gamma.kin = i.childhood.gamma.mat,
                     gamma.anc = i.gamma.anc, #new
                     childhood.gamma.anc = i.childhood.gamma.anc, #new
                     ##own.id.sg1 = i.own.id.sg1,  # DELETE
                     mom.id.sg1 = i.mom.id.sg1,
                     gmom.id.sg1 = i.gmom.id.sg1,
                     ggmom.id.sg1 = i.ggmom.id.sg1,
                     group.id.sg1 = i.group.id.sg1,
                     gamma.sg1 = i.gamma.sg1,
                     childhood.gamma.sg1 = i.childhood.gamma.sg1,
                     sharing.gamma = i.sharing.gamma,
                     sharing.childhood.gamma = i.sharing.childhood.gamma,
                     subPop.id = i.subPop.id,       #new vector to track subpopulations
                     gene = i.genes,
                     age = i.ages,      
                     max.age = as.integer(max.age) ,
                     e.fertility = e.fert,
                     e.mortality =  e.mort,
                     e.production =  e.prod,
                     e.cproduction=  e.cprod,
                     cproduction.cap = cprod.cap,
                     e.density =  e.dens,
                     #pop.size= length(i.own.id.mat),
                     resource.size = Cresources
                     );




rm(i.own.id.mat,i.mom.id.mat, i.gmom.id.mat, i.ggmom.id.mat,
   i.gggmom.id.mat, i.kinship3.id, i.kinship4.id, i.kinship5.id, i.group.id.mat, i.gamma.mat,
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
test$gene <- ( zero * test$gene);

## convert to integer to save memory and increase speed
#test<- f.setmode(test);

########################
## Forecast Ncycles years
########################
yx.series <- array(NA, dim=c(Ncycles,16,length(subPopTypes) ) );
sim.results <-  array(NA, dim=c(Ncycles,16,length(subPopTypes) ) );
pop.series <-  array(NA, dim=c(Ncycles,16,length(subPopTypes) ) );


gamma.series <- array(NA, dim=c(Ncycles,16,length(subPopTypes) ) );
childhood.gamma.series <- array(NA, dim=c(Ncycles,16,length(subPopTypes) ) );
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
                     cproduction.cap = cprod.cap,
                     e.density =  e.dens,
                     #pop.size= length(test$own.id.mat),
                     resource.size = Cresources
                     );


 
  #sim.results[year,] <- test$mean.harm
  gamma.series[year,,] <-  tapply(test$sharing.gamma,list(factor(test$age,levels=0:15),
                                  factor(test$subPop.id,levels=seq(subPopTypes) )),mean);
  childhood.gamma.series[year,,] <-  tapply(test$sharing.childhood.gamma,list(factor(test$age,levels=0:15),
                                                          factor(test$subPop.id,levels=seq(subPopTypes) )),mean);
  t.age.factor <- factor(test$age,levels=0:15);
  t.subpop.factor <- factor(test$subPop.id,levels=seq(subPopTypes));
for(i in seq(subPopTypes)){
  # avg mutation load, by age class, for each subpop
  t.subpop<- test$subPop.id==i;
  if(any(t.subpop)){                    #if there is a 0 subpop, output remains NA
    sim.results[year,,i] <-apply(test$gene[t.subpop,, drop=F],2,mean);
     anc.series[year,i] <-      mean( table(test$group.id.anc[t.subpop]));
    #kinship3.series[year,i] <- mean(table(test$kinship3.id[t.subpop]));
    #kinship4.series[year,i] <- mean(table(test$kinship4.id[t.subpop]));
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
 

 print(paste("cycle ",year));
  
 if( (year %% Ncheckpoint) == 0){
    print(paste(date(), "CHECKPOINT:  Completed ", start.year + year," of ", start.year + Ncycles));
    # save the results to this point
    save.image(file=".CheckPoint",compress=TRUE);
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
# final results saved
save.image(file=".CheckPoint",compress=TRUE);

