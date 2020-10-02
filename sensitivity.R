library(data.table)
library(lsr)
library(parallel)

maxgroups=104

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Usage: script binSize ncores",call.=FALSE)
} else {
  k=as.numeric(args[1])
  ncores=as.numeric(args[2])
}

load(file="sensitivityData.RData")

##Wrapper for ciMean so that it returns a real NA instead of a logical one. Needed for data.table
myCiMean=function(x,...){
    val=ciMean(x,...)[1]
    if(is.na(val)) {
        return(NA_real_)
    }else{
        return(val)
    }
}

##Done previously to prepare the data
#####################################
#samples=data[,unique(Sample)]
#nsamples=length(samples)
#sampleIs=seq(1,length(samples))
#dataSens=data[conditionPAF>0.05,.(Sample,condition,scorefiltNABcovB_PAFFPAF0.05)]
#toCompareWith=data[conditionPAF>0.05,.(score=first(scorefiltNABcovB_PAFFPAF0.05_ciMean_bycond)),by=condition][order(-score)]

groups=combn(sampleIs,k,simplify=FALSE)

if (length(groups)>=maxgroups) {
    groups=sample(groups,size=maxgroups,replace=FALSE)
}

resSens=mclapply(groups,function(thisgroup){
  thisres=dataSens[Sample%in%samples[thisgroup],.(score=myCiMean(scorefiltNABcovB_PAFFPAF0.05,conf=.9,na.rm=TRUE)[1]),keyby=condition][order(-score)]
  nvalue=which(toCompareWith[,condition]==thisres[1,condition])
  return(c(k,nvalue,ecdf(toCompareWith[,score])(toCompareWith[nvalue,score])))
},mc.cores=ncores) ##resSens will be a list now I may need to so something to it before exporting

resSens=do.call(rbind.data.frame,resSens)
colnames(resSens)=c("gsize","ncond","rAc")
write.csv(resSens,file=paste0("resultsSensitivity_",k,".csv"),quote=FALSE,row.names=FALSE)

