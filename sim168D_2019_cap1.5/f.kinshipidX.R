## Define the kinship relatedness measures
      ## Determine the kinship id.  The concept is pretty simple. The
      ## kinship.id is the id of the oldest matriarch back to the
      ## ggmom generation.  To determine this use the gggmom ids and
      ## assume that no gggmom is currently alive.  This assumption
      ## lets us avoid counting too far  back in time, double
      ## counting, etc. We want to follow the line of an individual
      ## back 3 generations, but if that individual is a mom, gmom, or
      ## ggmom, we do not go back so far, e.g. the kinship.id of
      ## living ggmom is herself. kinship.id equates 3rd cousins.

      ## Unlike the matriarchal measure, this measures across both
      ## living and dead ancestors. This is a measure over generations
      ## and not over time,

      ## lesson learned: keep separate track of mom,gmom,ggmom,gggmom
      ## arrays for kin,mat cases since the .mat history contains
      ## censoring of IDs from deaths while the .kin does not. It
      ## should only be necessary to track gggmom for the .kin case;
      ## it is not used otherwise

 
      ## Note that kinship.id is not necessarily a strictly heritable
      ## trait.  progeny of two sisters will eventually have distinct
      ## kinship.ids

f.kinshipid3 <- function(gggmom.id.kin, ggmom.id.kin, gmom.id.kin, mom.id.kin,own.id.kin){
  t.lgggm <- length(gggmom.id.kin);
  t.lggm  <- length(ggmom.id.kin);
  t.lgm   <- length(gmom.id.kin);
  t.lm    <- length(mom.id.kin);
  t.o     <- length(own.id.kin);
  if( ! all( c( t.lgggm == t.lggm, t.lgggm==t.lgm,t.lgggm==t.lm, t.lgggm==t.o ))){
    stop("arguments must be of equal length")
  }
      ## go three generations back to look for ancestor
  gggmoms<-sort(unique(gggmom.id.kin));
  t.maxid<- max(own.id.kin);
  isel <- (own.id.kin %in% gggmoms ); # should be empty.  Test
  if(any(isel)) warning("living gggmom!", immediate=TRUE)
  
  t.bumpids <- as.integer(100);
  isel<- mom.id.kin %in% gggmoms ;
  t.mom.id <- mom.id.kin;
  t.mom.id[isel]<- t.maxid + t.bumpids;   # replace any moms who are also gggmoms with out of range values
  
  t.gmom.id <- gmom.id.kin;
  isel<- isel | (t.gmom.id %in% gggmoms);
  t.gmom.id[isel] <- t.maxid + t.bumpids;   #replace gmoms in ggg moms or gmoms whose daughters were gggmoms


  t.ggmom.id <- ggmom.id.kin;
  isel<- isel | (t.ggmom.id %in% gggmoms);
  t.ggmom.id[isel] <- t.maxid+t.bumpids;

  t.tree<- rbind(own.id.kin,t.mom.id, t.gmom.id, t.ggmom.id);
  storage.mode(t.tree)<- "integer";
  ## we have replaced gggmom ids and those ids that are even older with high values. Look back
  ## across each line and choose the minimum value. That is the kinship.id
  
                                        # kinship.id is oldest matriarch in common, 3 gens back, suitably adjusted
                                        # For those with new history (e.g. from dispersal), the past is ignored
  kinship.id <- apply(t.tree, 2, min, na.rm=TRUE);
  return(kinship.id)
}

f.kinshipid4 <- function(ggggmom.id.kin, gggmom.id.kin, ggmom.id.kin, gmom.id.kin, mom.id.kin,own.id.kin){

  ## go 4 generations back to look for the ancestor
  ggggmoms<-sort(unique(ggggmom.id.kin));
  t.maxid<- max(own.id.kin);
  isel <- (own.id.kin %in% ggggmoms ); # should be empty.  Test
  if(any(isel)) warning("living ggggmom!", immediate=TRUE)

  t.bumpids <- as.integer(100);
  isel<- mom.id.kin %in% ggggmoms ;
  t.mom.id <- mom.id.kin;
  t.mom.id[isel]<- t.maxid + t.bumpids;   # replace any moms who are also ggggmoms with out of range values
      
  t.gmom.id <- gmom.id.kin;
  isel<- isel | (t.gmom.id %in% ggggmoms);
  t.gmom.id[isel] <- t.maxid + t.bumpids;   #replace gmoms in gggg moms or gmoms whose daughters were ggggmoms


  t.ggmom.id <- ggmom.id.kin;
  isel<- isel | (t.ggmom.id %in% ggggmoms);
  t.ggmom.id[isel] <- t.maxid+t.bumpids;


  t.gggmom.id <- gggmom.id.kin;
  isel<- isel | (t.gggmom.id %in% ggggmoms);
  t.gggmom.id[isel] <- t.maxid+t.bumpids;

  t.tree<- rbind(own.id.kin,t.mom.id, t.gmom.id, t.ggmom.id,t.gggmom.id);
  storage.mode(t.tree)<- "integer";
  ## we have replaced gggmom ids and those ids that are even older with high values. Look back
  ## across each line and choose the minimum value. That is the kinship.id
  
                                        # kinship.id is oldest matriarch in common, 3 gens back, suitably adjusted
                                        # For those with new history (e.g. from dispersal), the past is ignored
  kinship2.id <- apply(t.tree, 2, min, na.rm=TRUE);
  return(kinship2.id)
}

f.kinshipid5 <- function(gggggmom.id.kin,ggggmom.id.kin, gggmom.id.kin, ggmom.id.kin, gmom.id.kin, mom.id.kin,own.id.kin){
## go 5 generations back to look for the ancestor
  gggggmoms<-sort(unique(gggggmom.id.kin));
  t.maxid<- max(own.id.kin);
  isel <- (own.id.kin %in% gggggmoms ); # should be empty.  Test
  if(any(isel)) warning("living gggggmom!", immediate=TRUE)
  
  t.bumpids <- as.integer(100);
  isel<- mom.id.kin %in% gggggmoms ;
  t.mom.id <- mom.id.kin;
  t.mom.id[isel]<- t.maxid + t.bumpids;   # replace any moms who are also gggggmoms with out of range values
  
  t.gmom.id <- gmom.id.kin;
  isel<- isel | (t.gmom.id %in% gggggmoms);
  t.gmom.id[isel] <- t.maxid + t.bumpids;   #replace gmoms in ggggg moms or gmoms whose daughters were gggggmoms


  t.ggmom.id <- ggmom.id.kin;
  isel<- isel | (t.ggmom.id %in% gggggmoms);
  t.ggmom.id[isel] <- t.maxid+t.bumpids;


  t.gggmom.id <- gggmom.id.kin;
  isel<- isel | (t.gggmom.id %in% gggggmoms);
  t.gggmom.id[isel] <- t.maxid+t.bumpids;

  t.ggggmom.id <- ggggmom.id.kin;
  isel<- isel | (t.ggggmom.id %in% gggggmoms);
  t.ggggmom.id[isel] <- t.maxid+t.bumpids;



  t.tree<- rbind(own.id.kin,t.mom.id, t.gmom.id, t.ggmom.id,t.gggmom.id,t.ggggmom.id);
  storage.mode(t.tree)<- "integer";
  ## we have replaced gggmom ids and those ids that are even older with high values. Look back
  ## across each line and choose the minimum value. That is the kinship.id
  
                                        # kinship.id is oldest matriarch in common, 3 gens back, suitably adjusted
      # For those with new history (e.g. from dispersal), the past is ignored
  kinship3.id <- apply(t.tree, 2, min, na.rm=TRUE);

  return(kinship3.id)
}
