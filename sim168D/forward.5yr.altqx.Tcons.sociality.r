
## Sim166 Evolution of sociality -- density based on TotConsumption, N
## subgroups grow independently, except coupled through total density;
## The  pops follow one of Maternal, kin5,  kin5 50/50 kin/sg sharing,
## and so the only difference between the populations is from random
## initial configuration and demographic randomness of events.  ##

# forward.5yr.altqx.Tcons.sociality

## No dispersion
source("f.Timapply.R")


forward.5yr.altqx <-
  function (cycleno,
            beta.mat,
            beta.kin,
            beta.sg1,
            own.id.mat,
            mom.id.mat,
            gmom.id.mat,
            ggmom.id.mat,
            mom.id.kin,
            gmom.id.kin,
            ggmom.id.kin,
            gggmom.id.kin,
            ggggmom.id.kin,
            gggggmom.id.kin,
            kinship3.id,
            kinship4.id,
            kinship5.id,            
            group.id.anc,
            own.id.yr,
            mom.id.yr,
            gmom.id.yr,
            ggmom.id.yr,
            gggmom.id.yr,
            group.id.mat,
            gamma.mat,
            childhood.gamma.mat,
            gamma.kin,
            childhood.gamma.kin,
            gamma.anc,
            childhood.gamma.anc,
            #own.id.sg1,
            mom.id.sg1,
            gmom.id.sg1,
            ggmom.id.sg1,
            group.id.sg1,
            gamma.sg1,
            childhood.gamma.sg1,
            sharing.gamma,
            sharing.childhood.gamma,
            subPop.id,       #new vector to track subpopulations
            gene,
            age,
            max.age,
            e.fertility,
            e.mortality,
            e.production,
            e.cproduction,
            e.density = -0.2,
            #pop.size= length(own.id.mat);
            resource.size=1e5
            )
{
  ## Carl Boe
  ## Tim Miller
  ## July 23, 2005
  ## This function projects population one year forward in time.
   
   
  
    ############################################################################
    # Step 1. Assign degree of food sharing from Matriarchy versus Social Group.
    ############################################################################



  ## the sharing parameter gamma is a weighted average of the sharing
  ## within matriarchies, sharing within kinship lines, and sharing
  ## within social groups.  The betas add to 1.  When beta.sg1=1, then
  ## sharing only takes place within social groups.

  ## Data-check: does sharing.gamma*consumption sum to production?

  ## Gamma is measure of actual consumption relative to baseline
  ## consumption needs (age pattern taken from Ache).

  ## NB: gamma values persevere from the end of the previous cycle and
  ## are recalculated as the last stage of this program
  

      #####################################################################
      # Step 2. Define density as *effective* population size divided by resource size.
      #####################################################################

      # Resource size is an arbitary constant which sets the equilibrium
      # population size.  A larger value of resources leads to lower density and
      # increased productivity and a larger equlibrium population size.
      # Resource size can be set at a fixed value for the duration of the simulation
      # or alternatively, can be updated periodically (not too frequently!)
      # during the simulation to insure a final equilibrium
      # population size of about N individuals.

      # redefine 'population' to be consumption-weighted population across *all* subgroups
      totC <- sum(cxz[age+1]);          #sum of each i's age-indexed consumption
      density <- sum(totC) /resource.size;

      #############################################
      # Step 3.  Births in the year to those age x.
      #############################################

      #  gave.birth == 1, indicating which moms give birth
      adjust.fx <-  sharing.gamma^(e.fertility)
      gave.birth <- rep(0,length(age))
      gave.birth[ runif(length(age)) < (fx[age+1])*adjust.fx ] <- 1

      #Newborns inherit mom's genes.
      newborn.gene <-  gene[gave.birth==1,]

      # Mutations in genes
      mutants <- rep(0,length(newborn.gene))
      mutants[runif(length(newborn.gene))<mute.rate] <- 1 # spontaneous mutation
      newborn.gene <- newborn.gene + as.integer(mutants)
      ## storage.mode(newborn.gene) <- "integer"; # ?? needed ?

      ## Newborns are age 0
      newborn.age <- rep(as.integer(0),sum(gave.birth))
      newborn.qx  <- newborn.age;       #newborns survive the first interval

      # Newborns inherit mom's ids and gamma values, and association (subPop) values
      gave.birth.isel <- (gave.birth==1); # used repeatedly, so calc once
      newborn.ggmom.id.sg1 <- gmom.id.sg1[gave.birth.isel]
      newborn.gmom.id.sg1 <-  mom.id.sg1[gave.birth.isel]
      newborn.mom.id.sg1 <-   own.id.mat[gave.birth.isel]; # own.id for sg1 is same as for .mat
      newborn.group.id.sg1 <-   group.id.sg1[gave.birth.isel]
      newborn.subPop.id <- subPop.id[gave.birth.isel]; # congenital association meme
      
      #newborn.own.id.sg1   <- as.integer(max(own.id.sg1) + seq(1,length(newborn.age)) )
      newborn.gamma.sg1 <-    gamma.sg1[gave.birth.isel]
      newborn.ggmom.id.mat <- gmom.id.mat[gave.birth.isel]
 ##     newborn.gggmom.id.mat <- ggmom.id.mat[gave.birth.isel]
      newborn.gmom.id.mat <-  mom.id.mat[gave.birth.isel]
      newborn.mom.id.mat <-   own.id.mat[gave.birth.isel]
      newborn.group.id.mat <-   group.id.mat[gave.birth.isel]
      newborn.own.id.mat   <- as.integer(max(own.id.mat) + seq(newborn.age))      
      newborn.gamma.mat <-    gamma.mat[gave.birth.isel]
      newborn.gamma.kin <-    gamma.kin[gave.birth.isel]
      # Tim's .anc measure
      newborn.gamma.anc <-    gamma.anc[gave.birth.isel]
      newborn.group.id.anc <-   group.id.anc[gave.birth.isel]
      newborn.sharing.gamma <- sharing.gamma[gave.birth.isel]
      newborn.sharing.childhood.gamma <- sharing.childhood.gamma[gave.birth.isel]
      

 ## new stuff for kin-based measurements      
      newborn.mom.id.kin <-    newborn.mom.id.mat; # newborns get mom ids regardless of .mat, .kin, .sg1
      newborn.gmom.id.kin <-  mom.id.kin[gave.birth.isel]
      newborn.ggmom.id.kin <-  gmom.id.kin[gave.birth.isel]
      newborn.gggmom.id.kin <-  ggmom.id.kin[gave.birth.isel]
      newborn.ggggmom.id.kin <-  gggmom.id.kin[gave.birth.isel]
      newborn.gggggmom.id.kin <-  ggggmom.id.kin[gave.birth.isel]
      
      #newborn.kinship.id  <-  kinship.id[gave.birth.isel] # assignment not needed, recomputed below
      
      # newborn birth year info
      newborn.own.id.yr   <-   as.integer(rep(cycleno,length(newborn.age)) );
      newborn.mom.id.yr   <-   own.id.yr[gave.birth.isel]
      newborn.gmom.id.yr   <-  mom.id.yr[gave.birth.isel]
      newborn.ggmom.id.yr   <-  gmom.id.yr[gave.birth.isel]
      newborn.gggmom.id.yr   <-  ggmom.id.yr[gave.birth.isel] 

      
      ###################################################
      # STEP 4.  Survive the population from age x to x+1.
      ####################################################

      # Some people die
      # effect of phenotype
      # mortality consists of two separable components, working on the hazard.  There is
      # (1) a component indexing the lethality of genetic mutations;
      # (2) a gamma-dependent component that indicates nutrition or other effect of resources
      
      lethality <- additive.gene.risk ; # CAB: use more appropriate name for this setting under


                                   # new mortality formulation
      gene.risk <- rep( as.integer(0),length(age));
      for (cnt in seq(gene.risk)){
        # get mutation load for each ind. at that individual's current age.
        gene.risk[cnt] <-  ( gene[cnt,][[(age[cnt]+1)]] ); # CAB: just the mutation count
      };

      ## new mortality formulation, see RLee memo 10/2/06 and subsequent email.  New
      ## hazard formulation is h(x) = K(x)*delta + background.mort * gamma^alpha (K is mutation count, delta is lethality,
      ##  epsilon is the level of background mortality, gamma is a resource, e.g. food consumption; epsilon
      ##  set to .01 corresponds to a background 1qx of .01 
        
      ## the baseline hazard is the -log(1-q_base)/5 for 5yr age groups.  But, we keep it in the qx form
      ## to avoid problems with log(0) or an infinity hazard rate.
      adjust.hx <- (gene.risk*lethality  + background.mort ) * sharing.gamma^(e.mortality) ;
       # incorporate baseline with gene-dependent mortality
      qx.risk <- 1 - (1 - qx.baseline[age+1])*( exp( -5*adjust.hx)  ) ;

      died <- (runif(length(age)) <= qx.risk ); #logical
      died.id <- own.id.mat[died] # ID of those who died

## new stuff: we keep a history of the death probabilities so we can look at a distribution over time
## and get some feel for the variation

      died.by.age.subpop <- tapply(died,
                                   list(factor(age,levels=0:15),factor(subPop.id,levels=seq(subPopTypes))),
                                   sum);
      pop.by.age.subpop <- tapply(rep(1,length(died)),
                                   list(factor(age,levels=0:15),factor(subPop.id,levels=seq(subPopTypes))),
                                   sum);

      qx.by.age.subpop <- ifelse(pop.by.age.subpop==0, NA, died.by.age.subpop/pop.by.age.subpop);
      
      #died.by.age <- timapply2(died, age,sum);
      #pop.by.age  <- timapply2(rep(1,length(died)),age,sum);
      #qx.by.age   <- died.by.age / pop.by.age; # cohort mortality, a 16-vector
      #qx.by.age <- ifelse( is.nan(qx.by.age),1,qx.by.age);
      #On death, remove rows of those who died

      age <- age[!died]
      gene <- gene[!died,]
      gamma.sg1 <- gamma.sg1[!died]
      childhood.gamma.sg1 <- childhood.gamma.sg1[!died]
      gamma.mat <- gamma.mat[!died]
      childhood.gamma.mat <- childhood.gamma.mat[!died]
      gamma.kin <- gamma.kin[!died]
      childhood.gamma.kin <- childhood.gamma.kin[!died]
      childhood.gamma.anc <- childhood.gamma.anc[!died]
      gamma.anc <- gamma.anc[!died]
      sharing.gamma <- sharing.gamma[!died]
      sharing.childhood.gamma <- sharing.childhood.gamma[!died]

      age <- age+ as.integer( 1 ) # age population by 1 year
      age[age>max.age] <- max.age

      qx.risk <- qx.risk[!died];
      
      ##own.id.sg1 <- own.id.sg1[!died]   ## DELETE THIS
      mom.id.sg1 <- mom.id.sg1[!died]
      gmom.id.sg1 <- gmom.id.sg1[!died]
      ggmom.id.sg1 <- ggmom.id.sg1[!died]
      group.id.sg1 <- group.id.sg1[!died]
      own.id.mat <- own.id.mat[!died]
      mom.id.mat <- mom.id.mat[!died]
      gmom.id.mat <- gmom.id.mat[!died]
      ggmom.id.mat <- ggmom.id.mat[!died]
      ##gggmom.id.mat <- gggmom.id.mat[!died]
      group.id.mat <- group.id.mat[!died]
      subPop.id <- subPop.id[!died];
      
      ## new stuff for kin-based measurements      
      mom.id.kin <-  mom.id.kin[!died]      
      gmom.id.kin <-  gmom.id.kin[!died]   
      ggmom.id.kin <- ggmom.id.kin[!died]
      gggmom.id.kin <- gggmom.id.kin[!died]
      ggggmom.id.kin <- ggggmom.id.kin[!died]
      gggggmom.id.kin <- gggggmom.id.kin[!died]
      
      ##kinship.id recomputed later

      ## Tim's ancestor id
      group.id.anc <- group.id.anc[!died]

      
      own.id.yr <- own.id.yr[!died]
      mom.id.yr <- mom.id.yr[!died]
      gmom.id.yr <- gmom.id.yr[!died]
      ggmom.id.yr <- ggmom.id.yr[!died]
      gggmom.id.yr <- gggmom.id.yr[!died] 

      ##ADD in NEWBORNS
      age <- c(age,newborn.age)
      gene <- rbind(gene,newborn.gene)
      ##own.id.sg1 <- c(own.id.sg1,newborn.own.id.sg1)  # DELETE 
      mom.id.sg1 <- c(mom.id.sg1,newborn.mom.id.sg1)
      gmom.id.sg1 <- c(gmom.id.sg1,newborn.gmom.id.sg1)
      ggmom.id.sg1 <- c(ggmom.id.sg1,newborn.ggmom.id.sg1)
      group.id.sg1 <- c(group.id.sg1,newborn.group.id.sg1)
      gamma.sg1 <- c(gamma.sg1,newborn.gamma.sg1)
      childhood.gamma.sg1 <- c(childhood.gamma.sg1,newborn.gamma.sg1)
      own.id.mat <- c(own.id.mat,newborn.own.id.mat)
      mom.id.mat <- c(mom.id.mat,newborn.mom.id.mat)
      gmom.id.mat <- c(gmom.id.mat,newborn.gmom.id.mat)
      ggmom.id.mat <- c(ggmom.id.mat,newborn.ggmom.id.mat)
      ##gggmom.id.mat <- c(gggmom.id.mat,newborn.gggmom.id.mat)
      group.id.mat <- c(group.id.mat,newborn.group.id.mat)
      subPop.id  <- c(subPop.id, newborn.subPop.id);
      gamma.mat <- c(gamma.mat,newborn.gamma.mat)
      childhood.gamma.mat <- c(childhood.gamma.mat,newborn.gamma.mat)
      
      ##new stuff for kin-based measurements
      mom.id.kin <- c(mom.id.kin,newborn.mom.id.kin)
      gmom.id.kin <- c(gmom.id.kin,newborn.gmom.id.kin)
      ggmom.id.kin <- c(ggmom.id.kin,newborn.ggmom.id.kin)
      gggmom.id.kin <- c(gggmom.id.kin,newborn.gggmom.id.kin)
      ggggmom.id.kin <- c(ggggmom.id.kin,newborn.ggggmom.id.kin)
      gggggmom.id.kin <- c(gggggmom.id.kin,newborn.gggggmom.id.kin)
      
      gamma.kin <- c(gamma.kin,newborn.gamma.kin)
      childhood.gamma.kin <- c(childhood.gamma.kin,newborn.gamma.kin)

      # Tim's ancestor measure
      gamma.anc <- c(gamma.anc,newborn.gamma.anc)
      childhood.gamma.anc <- c(childhood.gamma.anc,newborn.gamma.anc)
      sharing.gamma <- c(sharing.gamma, newborn.sharing.gamma)
      sharing.childhood.gamma <- c(sharing.childhood.gamma, newborn.sharing.childhood.gamma);
      
      group.id.anc <- c(group.id.anc,newborn.group.id.anc)

      
      own.id.yr <- c(own.id.yr,newborn.own.id.yr)
      mom.id.yr <- c(mom.id.yr,newborn.mom.id.yr)
      gmom.id.yr <- c(gmom.id.yr,newborn.gmom.id.yr)
      ggmom.id.yr <- c(ggmom.id.yr,newborn.ggmom.id.yr)
      gggmom.id.yr <- c(gggmom.id.yr,newborn.gggmom.id.yr)
      
      #########################################################################
      # STEP 5.  Update kin groups and matriarchies 
      #########################################################################



      ##
      kinship3.id <- NULL; kinship4.id <- NULL; kinship5.id <- NULL;
      #kinship3.id <-  f.kinshipid3(gggmom.id.kin, ggmom.id.kin, gmom.id.kin, mom.id.kin,own.id.mat)
      #kinship4.id <- f.kinshipid4(ggggmom.id.kin, gggmom.id.kin, ggmom.id.kin, gmom.id.kin, mom.id.kin,own.id.mat)
      kinship5.id <- f.kinshipid5(gggggmom.id.kin,ggggmom.id.kin, gggmom.id.kin, ggmom.id.kin, gmom.id.kin, mom.id.kin,own.id.mat)

      

      # Older measure of matriarchy; for reference
      #On death, replace ids used to identify matriarchies with NAs

      mom.id.mat[mom.id.mat %in% died.id] <- NA
      gmom.id.mat[gmom.id.mat %in% died.id] <- NA
      ggmom.id.mat[ggmom.id.mat %in% died.id] <- NA

      group.id.mat <- ggmom.id.mat
      group.id.mat[is.na(group.id.mat)] <- gmom.id.mat[is.na(group.id.mat)]
      group.id.mat[is.na(group.id.mat)] <- mom.id.mat[is.na(group.id.mat)]
      group.id.mat[is.na(group.id.mat)] <- own.id.mat[is.na(group.id.mat)]

      #########################################################################
      # STEP 5b.  Calculate ancestory ids
      #########################################################################
      # ancestor id is the grand grandmom of the matriarch

      ## Need a sort of look-up table to match matriarchy id with anc id (which is ggmom of matriarchs)
      is.matriarch <- group.id.mat == own.id.mat # T = matriarch
      matriarchies <- group.id.mat[is.matriarch==TRUE]
      ancestories <-  ggmom.id.kin[is.matriarch==TRUE]; # ggmom.id.kin is the same as tim's ggmom.id.anc

      locate <- match(group.id.mat, matriarchies)
      group.id.anc <- ancestories[locate]



      ######################################################################
      # STEP 5.5 Randomly reassign kinship members into different groups
      ######################################################################

      ## Because kin groups are always in in the same subpopulation, we can operate at the population level
      
      # Every n cycles, randomly put all members of a kin group into a different social group. 
      ## no dispersion, set the freq to 0 ##
    if(setDispersal !=  'useNone'){        #skip this section
      disperse.freq <- 0;               # freq in cycles of  the dispersion event; 0 disables 
      if(disperse.freq != 0 &&  (cycleno %% disperse.freq) == 0 ){
      # List of social group foreach kin line
        weighted.list.of.sg.ids <- group.id.sg1[!duplicated(kinship5.id)]

      #Randomly re-assign this list
        random.list.of.sg.ids <- sample(weighted.list.of.sg.ids,
                                        length(weighted.list.of.sg.ids),
                                        replace=T)
      # Create a look-up table to match kinship id with new
      # randomly assign sg1
        random.id.sg1 <- random.list.of.sg.ids[match(kinship5.id,unique(kinship5.id))]

      # Reassign to this new social group
        group.id.sg1 <- as.integer(random.id.sg1);
      }
    };
      
      ################################################################
      # STEP 6  Fission (splitting along maternal ancestor lines)
      ################################################################
      
      Pfissions=NA;
   
      ## Splitting happens when kinship lines (as measured by
      ## kinship.id) get too big. Unless splitting or dispersion happens, the
      ## social group is mostly the maternal group.
      
      #print("...Step 5.5");
      fission.size <- 25   

    if(setFission=='useTMpk5fission'){
      
      ## kinship5 groups are kept together in the splits, and kin5
      ## groups do not overlap subpop boundaries. So, there is no need
      ## to iterate over subpops, but we do anyway, because we want the subpop detail

      Pfissions <- NA + seq(subPopTypes);# container of dim=#subpops
      t.subPops <- seq(subPopTypes);
      for( iSub in (t.subPops) ){
        isel <- (t.subPops[iSub] == subPop.id ) ;
        t.thistype <- subPopTypes[iSub]; # regime type of this subpop
        t.group.id.sg1 <- group.id.sg1[isel]; #subpop specific group

        ## Find Big Groups
        group.size <- table(t.group.id.sg1)
        fission.size <- 25
        big.group.size <- group.size[group.size > fission.size]
        big.groups <- as.numeric(names(big.group.size))

        ##Create new group IDs for big groups.
        new.big.group.id.sg1 <- max(group.id.sg1)+seq(big.groups)

        ##Fission the groups within subpops
        new.group.id.sg1 <- t.group.id.sg1
        locate <- match(new.group.id.sg1,big.groups)
        new.group.id.sg1[!is.na(locate)] <- new.big.group.id.sg1[locate[!is.na(locate)]]
        ## Re-assign about half of population to original groups if kinship id is an odd integer.
        ## This keeps kinships from being split when the social group fissions.
        t.kinship <- kinship5.id[isel];
        odd.integer <- (t.kinship  > (2*floor( t.kinship /2)) )
        new.group.id.sg1[odd.integer] <- t.group.id.sg1[odd.integer]
        group.id.sg1[isel] <- as.integer(new.group.id.sg1);
        Pfissions[iSub] <- length(big.group.size)/length(group.size);
      }
      
    }; #end if setFission
 

      ##########################################################
      # STEP 7.  Fuse social groups which have too few members.
      ##########################################################

      ## The strategy here is similar to the fission computation.  We
      ## order groups defined by group.id.sg1 according to size, with
      ## the smallest first.  Next process all the groups with size
      ## falling below the fusion.size threshold. For each small
      ## group, fuse them to the group which has the highest
      ## relatedness index (based on kinship.id, gggmom.id, etc.) and
      ## which can accomodate the small group.  
      
      ##print("...Step 6 fusion");     
      # Find Small Groups

      Pfusions = NA;
      

    if(setFusion=='useTMfusion'){
      Pfusions <- NA+ seq(subPopTypes);# container of dim=#subpops
      
      ## Find Small Groups

      
      t.subPops <- seq(subPopTypes);
      for( iSub in (t.subPops) ){
        isel <- (t.subPops[iSub] == subPop.id ) ;
        t.group.id.sg1 <- group.id.sg1[isel]; #subpop specific group
        group.size <- table(t.group.id.sg1)
        fusion.size <- 8
        small.group.size <- group.size[group.size<fusion.size]
        small.groups <- as.numeric(names(small.group.size))
      
        if(length(small.groups)>1){ # Must have at least 2 groups
          ## Combine small groups two at a time.
          ## Take half of the ids of small groups.
          ## These groups will be eliminated.
          small.groups.eliminated <- small.groups[seq(2,length(small.groups),2)]
          ## Take the other half of the id of small groups.
          ## These groups will receive new members.
          small.groups.receive <- small.groups[seq(1,length(small.groups),2)]

          ## Assign new group ids to each individual
          new.group.id.sg1 <- t.group.id.sg1;
          ## Find those in small groups which are to be eliminated
          locate <- match(t.group.id.sg1,small.groups.eliminated)
          ## Assign them to other small groups
          new.group.id.sg1[!is.na(locate)] <- small.groups.receive[locate[!is.na(locate)]]
          ## we have reassigned groups within the subpopulation; write out the results
          ## but just for the subpopulation
          group.id.sg1[isel] <- as.integer(new.group.id.sg1);
          Pfusions[iSub] <- length(small.group.size)/length(group.size);
        }
      } # end for(iSub
      
    } #end setFusion 

      ################################################################################################
      # STEP 8.  Update economic variables based on new composition of matriarchies and social groups.
      ################################################################################################
      ##print("...Step 8");

      ### all of the gamma computations involve restrictions to the subpopulation level.
      ### Therefore, we loop through the entire gamma calculation set at one time

      ## production is a function of prior sharing gammas which are passed into this routine


      
## begin loop      
      t.subPops <- seq(subPopTypes);
      yx <- NA + (own.id.mat);
      
      for( iSub in (t.subPops) ){
        isel <- (t.subPops[iSub] == subPop.id ) ;
        t.group.id.sg1 <- group.id.sg1[isel]; #subpop specific group
        t.group.id.mat <- group.id.mat[isel]; #subpop specific group
        t.age <- age[isel];             #subpop specific ages
        t.scg <-  sharing.childhood.gamma[isel]
        t.sg <-  sharing.gamma[isel]
        t.kinship5.id <- kinship5.id[isel];
        t.group.id.anc <- group.id.anc[isel];
      
      ##################################################################################
      # STAGE 8a. Find gamma (relative consumption) based on membership in social group.
      ##################################################################################

     # Sum of consumption weights within each social group, specific within each subpopulation
        c.group <- tapply(cxz[t.age+1],t.group.id.sg1,sum)
        g.group.ids <- as.integer(names(c.group))
        c.group <- as.numeric(c.group)

        ## Sum of production weights within each social group
        t.yx <- yxz[t.age+1]*(t.scg^(e.cproduction))*(t.sg^(e.production))
        t.yx <- t.yx * (density^e.density); ## !! GLOBAL density
        p.group <- as.numeric(tapply( t.yx,t.group.id.sg1,sum))

        ## Define gamma parameter for each subpopulation
        g.group <- p.group/c.group

        ## Assign gamma to members of group
        gamma.sg1[isel] <- g.group[match(t.group.id.sg1,g.group.ids)]
      

      ################################################################################
      # Step 8b.  Find gamma (relative consumption) based on membership in matriarchy.
      #################################################################################

      ## Sum of consumption weights within each matriarchy
        c.group <- tapply(cxz[t.age+1],t.group.id.mat,sum)
        g.group.ids <- as.integer(names(c.group))
        c.group <- as.numeric(c.group)

      ## Sum of production weights within each social group -- already computed
        #yx <- yxz[t.age+1]*(t.scg ^(e.production))*(t.sg^(e.production))
        #yx <- yx * (density^e.density)
        p.group <- as.numeric(tapply(t.yx,t.group.id.mat,sum))

      ## Define gamma parameter for each social group
        g.group <- p.group/c.group

      ## Assign gamma to members of group
        gamma.mat[isel] <- g.group[match(t.group.id.mat,g.group.ids)]

      ################################################################################
      # Step 8c.  Find gamma (relative consumption) based on membership in kinship group
      #################################################################################

        ## Sum of consumption weights within each kinship group
        c.group <- tapply(cxz[t.age+1],t.kinship5.id,sum)
        g.group.ids <- as.integer(names(c.group))
        c.group <- as.numeric(c.group)

        ## Sum of production weights within each kinship group -- already computed above
        #yx <- yxz[t.age+1]*(t.scg^(e.production))*(t.sg ^(e.production))
        #yx <- yx * (density^e.density)
        p.group <- as.numeric(tapply( t.yx,t.kinship5.id,sum))

        ## Define gamma parameter for each kinship group
        g.group <- p.group/c.group

        ## Assign gamma to members of kinship group
        gamma.kin[isel] <- g.group[match(t.kinship5.id,g.group.ids)]

     #################################################################################
     # Step 8c.  Find gamma (relative consumption) based on membership in ancestory group.
     #################################################################################


        ## Sum of consumption weights within each related group
        c.group <- tapply(cxz[t.age+1],t.group.id.anc,sum)
        g.group.ids <- as.numeric(names(c.group))
        c.group <- as.numeric(c.group)

        ## Sum of production weights within each  group  -- already computed
        #yx <- yxz[t.age+1]*(t.scg^(e.production))*(t.sg^(e.production))
        #yx <- yx * (density^e.density)
        p.group <- as.numeric(tapply( t.yx,t.group.id.anc,sum))

        ## Define gamma parameter for each  group
        g.group <- p.group/c.group

        ## Assign gamma to members of group
        gamma.anc[isel] <- g.group[match(t.group.id.anc,g.group.ids)]

        ## save the production information for the subgroup into the whole container
        yx[isel] <- t.yx;

      } # for(iSub ...

     #################################################################################
     ## Step 8d.  Compute sharing gamma based on component sg,kin,mat gammas
     #################################################################################
      ##
      ## each subpopulation has its own type of sharing gamma calculation;
      ## so do each subpop separately
      
      t.subPops <- seq(subPopTypes)
      for( iSub in (t.subPops) ){
        isel <- (t.subPops[iSub] == subPop.id ) ;
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

        sharing.gamma[isel] <- ( t.beta.sg1 *gamma.sg1[isel] ) + (t.beta.kin * gamma.kin[isel] ) +
          ( t.beta.mat*gamma.mat[isel] );
        sharing.childhood.gamma[isel] <- ( t.beta.sg1*childhood.gamma.sg1[isel] ) +
          (t.beta.kin * childhood.gamma.kin[isel]) + ( t.beta.mat *childhood.gamma.mat[isel] );
      };
 
      
      
     ###############################
     # Step 9.  Summarize the data.
     ###############################

     # average number of mutant genes carried by newborns
      mean.harm <- apply(newborn.gene,2,mean)

 
      result <- list(own.id.mat=own.id.mat,
                     mom.id.mat=mom.id.mat,
                     gmom.id.mat=gmom.id.mat,
                     ggmom.id.mat=ggmom.id.mat,
                     ##gggmom.id.mat=gggmom.id.mat,
                       
                     mom.id.kin = mom.id.kin,
                     gmom.id.kin = gmom.id.kin,
                     ggmom.id.kin = ggmom.id.kin,
                     gggmom.id.kin = gggmom.id.kin ,
                     ggggmom.id.kin = ggggmom.id.kin ,
                     gggggmom.id.kin = gggggmom.id.kin ,
                     kinship3.id = kinship3.id,
                     kinship4.id = kinship4.id,
                     kinship5.id = kinship5.id,                     
                     group.id.anc=group.id.anc, #Tim's ancestory grouping variable
                     own.id.yr = own.id.yr,
                     mom.id.yr = mom.id.yr,
                     gmom.id.yr = gmom.id.yr,
                     ggmom.id.yr=ggmom.id.yr,
                     gggmom.id.yr = gggmom.id.yr,                                       
                     group.id.mat = group.id.mat,
                     gamma.mat = gamma.mat,
                     childhood.gamma.mat = childhood.gamma.mat,
                     gamma.kin = gamma.kin,
                     childhood.gamma.kin = childhood.gamma.kin, 
                     gamma.anc = gamma.anc, #new
                     childhood.gamma.anc = childhood.gamma.anc, #new
                     ##own.id.sg1 = own.id.sg1,  # DELETE
                     mom.id.sg1 = mom.id.sg1,
                     gmom.id.sg1 = gmom.id.sg1,
                     ggmom.id.sg1 = ggmom.id.sg1,                                  
                     group.id.sg1 = group.id.sg1,
                     gamma.sg1 = gamma.sg1,
                     childhood.gamma.sg1 = childhood.gamma.sg1,
                     sharing.gamma = sharing.gamma,
                     sharing.childhood.gamma = sharing.childhood.gamma,     
                     subPop.id=subPop.id,
                     gene=gene,
                     age=age,                   
                     mean.harm=mean.harm,
                     yx=yx,
                     c.group=c.group,
                     p.group=p.group,
                     g.group=g.group,
                     Pfissions=Pfissions,
                     Pfusions=Pfusions,
                     qx.by.age=qx.by.age.subpop
                     )
      return(result)
  }           
                                
         
          
