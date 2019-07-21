# Tim Miller
# Oct 27, 2005
# I added Carl's code for calculating means to the graphs.

# Count Number of Matriarchies

mat.count <- table(table(test$group.id.mat))
ind.mat.count <- mat.count*as.numeric(names(mat.count))
sum(mat.count) # number of matrirachies
sum(ind.mat.count) # number of individuals

# Count Number of Social Groups

sg1.count <- table(table(test$group.id.sg1))
ind.sg1.count <- sg1.count*as.numeric(names(sg1.count))
sum(sg1.count) # number of social groups
sum(ind.sg1.count) # number of individuals 

## Count Number of TM Ancestor Groups

anc.count <- table(table(test$group.id.anc))
ind.anc.count <- anc.count*as.numeric(names(anc.count))
sum(anc.count) # number of anc
sum(ind.anc.count) # number of individuals

## Count Number of kin3 groups
kin3.count <- table(table(test$kinship3.id))
ind.kin3.count <- kin3.count*as.numeric(names(kin3.count))
sum(kin3.count) # number of anc
sum(ind.kin3.count) # number of individuals

## Count Number of kin4 groups
kin4.count <- table(table(test$kinship4.id))
ind.kin4.count <- kin4.count*as.numeric(names(kin4.count))
sum(kin4.count) # number of anc
sum(ind.kin4.count) # number of individuals
## Count Number of kin5 groups
kin5.count <- table(table(test$kinship5.id))
ind.kin5.count <- kin5.count*as.numeric(names(kin5.count))
sum(kin5.count) # number of anc
sum(ind.kin5.count) # number of individuals




tim.barplot<- function(input,...){
  bp.max <- max(input)
  offset <- .015*bp.max
  bp <- barplot(input,space=0,ylim=c(0,1.03*bp.max),...)
  text(bp,input+offset,input,col="blue")
  tim.stamp()
}

#ps.options(paper="letter")
#postscript("counts.ps")
t.mean =
sum(mat.count*as.numeric(dimnames(mat.count)[[1]]))/sum(mat.count);

tim.barplot(mat.count,
            xlab="Size of Matriarchy",
            ylab="Number of Matriarchies",
            main=paste(sum(mat.count),"Matriarchies by Size of Matriarchy\n Mean=",round(t.mean,3)))

t.mean =
sum(ind.mat.count*as.numeric(dimnames(ind.mat.count)[[1]]))/sum(ind.mat.count);
tim.barplot(ind.mat.count,
            xlab="Size of Matriarchy",
            ylab="Number of Individuals",
            main=paste(sum(ind.mat.count),"Individuals by Size of Matriarchy in Which They Live\n Mean=",round(t.mean,3)));

tim.barplot(round(100*mat.count/sum(mat.count),2),
            xlab="Size of Matriarchy",
            ylab="Percent of Matriarchies",
            main="Percent of Matriarchies by Size of Matriarchy")
tim.barplot(round(100*ind.mat.count/sum(ind.mat.count),2),
            xlab="Size of Matriarchy",
            ylab="Percent of Individuals",
            main="Percent of Individuals by Size of Matriarchy in Which They Live")

t.mean =
sum(sg1.count*as.numeric(dimnames(sg1.count)[[1]]))/sum(sg1.count);
tim.barplot(sg1.count,
            xlab="Size of Social Group",
            ylab="Number of Social Groups",
            main=paste(sum(sg1.count),"Social Groups by Size of Social Group\n Mean=",round(t.mean,3)))

t.mean =
sum(ind.sg1.count*as.numeric(dimnames(ind.sg1.count)[[1]]))/sum(ind.sg1.count);

tim.barplot(ind.sg1.count,
            xlab="Size of Social Group",
            ylab="Number of Individuals",
            main=paste(sum(ind.sg1.count),"Individuals by Size of Social Group in Which They Live\n Mean=",round(t.mean,3)))
tim.barplot(round(100*sg1.count/sum(sg1.count),2),
            xlab="Size of Social Group",
            ylab="Percent of Social Groups",
            main="Percent of Social Groups by Size of Social Group")
tim.barplot(round(100*ind.sg1.count/sum(ind.sg1.count),2),
            xlab="Size of Social Group",
            ylab="Percent of Individuals",
            main="Percent of Individuals by Size of Social Group in Which They Live")


par.old<-par(no.readonly=T);
par(mfrow=c(2,2));

t.mean =
sum(anc.count*as.numeric(dimnames(anc.count)[[1]]))/sum(anc.count);
tim.barplot(anc.count,
            xlab="Size of TM Ancestry",
            ylab="Number ",
            main=paste(sum(anc.count),"TM Ancestry Groups by Size of Group\n Mean=",round(t.mean,3)))

t.mean =
sum(kin3.count*as.numeric(dimnames(kin3.count)[[1]]))/sum(kin3.count);
tim.barplot(kin3.count,
            xlab="Size of Kinship3 Group",
            ylab="Number ",
            main=paste(sum(kin3.count),"CB Kin3 Groups by Size of Group\n Mean=",round(t.mean,3)))

t.mean =
sum(kin4.count*as.numeric(dimnames(kin4.count)[[1]]))/sum(kin4.count);
tim.barplot(kin4.count,
            xlab="Size of Kinship4 Group",
            ylab="Number ",
            main=paste(sum(kin4.count),"CB Kin4 Groups by Size of Group\n Mean=",round(t.mean,3)))

t.mean =
sum(kin5.count*as.numeric(dimnames(kin5.count)[[1]]))/sum(kin5.count);
tim.barplot(kin5.count,
            xlab="Size of Kinship5 Group",
            ylab="Number ",
            main=paste(sum(kin5.count),"CB Kin5 Groups by Size of Group\n Mean=",round(t.mean,3)))



t.mean =
sum(ind.anc.count*as.numeric(dimnames(ind.anc.count)[[1]]))/sum(ind.anc.count);
tim.barplot(ind.anc.count,
            xlab="Size of TM Ancestry",
            ylab="Number of Individuals",
            main=paste(sum(ind.anc.count),"Individuals by Size of Group in Which They Live\n Mean=",round(t.mean,3)));

t.mean =
sum(ind.kin3.count*as.numeric(dimnames(ind.kin3.count)[[1]]))/sum(ind.kin3.count);
tim.barplot(ind.kin3.count,
            xlab="Size of Kinship3 Group",
            ylab="Number of Individuals",
            main=paste(sum(ind.kin3.count),"Individuals by Size of Group in Which They Live\n Mean=",round(t.mean,3)));

t.mean =
sum(ind.kin4.count*as.numeric(dimnames(ind.kin4.count)[[1]]))/sum(ind.kin4.count);
tim.barplot(ind.kin4.count,
            xlab="Size of Kinship4 Group",
            ylab="Number of Individuals",
            main=paste(sum(ind.kin4.count),"Individuals by Size of Group in Which They Live\n Mean=",round(t.mean,3)));

t.mean =
sum(ind.kin5.count*as.numeric(dimnames(ind.kin5.count)[[1]]))/sum(ind.kin5.count);
tim.barplot(ind.kin5.count,
            xlab="Size of Kinship5 Group",
            ylab="Number of Individuals",
            main=paste(sum(ind.kin5.count),"Individuals by Size of Group in Which They Live\n Mean=",round(t.mean,3)));
                                        #dev.off()
par(par.old);

