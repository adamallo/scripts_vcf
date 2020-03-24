library(lsr)
library(methods)
library(plotly)
library(cowplot)
library(xtable)
library(dunn.test)
library(boot)
library(dplyr)
library(doParallel)
library(Hmisc)
library(tidyr)
library(data.table)

n_cores=4 ##WARNING:Only use more than one if the sytem has a lot of RAM!

# Data parsing #HARDCODED
############################################

###Conf variables
#dir="~/dcis/new_replicates_withpaf_snv_indel_newcovB_newparseTSV"
n_sim_param_simple=6###WARNING: CHANGE THIS ACCORDINGLY WITH YOUR DATA. TOTAL NUMBER OF PARAMETERS non-NAB and non-covB
n_sim_param_covB=3###WARNING: CHANGE THIS ACCORDINGLY WITH YOUR DATA. PARAMETERS RELATIVE TO covB
n_sim_param_NAB=3###WARNING: CHANGE THIS ACCORDINGLY WITH YOUR DATA. PARAMETERS RELATIVE TO NAB
n_sim_param_PAF=1###WARNING: CHANGE THIS ACCORDINGLY WITH YOUR DATA. PARAMETERS RELATIVE TO popAF
n_sim_params=n_sim_param_simple+n_sim_param_covB+n_sim_param_NAB+n_sim_param_PAF
#setwd(dir)
#################


###My functions
#################################################

makekgroups=function(thisdata,thisformula="Sample ~ DNAcon",k=5){
  groupsizes=aggregate(as.formula(thisformula),thisdata,length)
  ntimes=min(groupsizes[,2]%/%k)
  if (ntimes == 0) {
    stop("Group subsampling has not been implemented yet")
  }
  
  ksamples=list()
  for(i in 1:k) {
    thisthisdata=groupsizes
    thisthisdata[,2]=rep(ntimes,n=nrow(groupsizes))    
    ksamples[[i]]=thisthisdata
  }
  
  remainders=groupsizes[,2]%%(ntimes*k)
  groupextras=sapply(remainders,function(x){sample(seq(1,k,1),size = x,replace = FALSE)})
  
  for (dnaconi in 1:length(groupextras)){
    for (reptoadd in groupextras[[dnaconi]]){
      ksamples[[reptoadd]][dnaconi,2]=ksamples[[reptoadd]][dnaconi,2]+1
    }
  }
  groups=aggregate(as.formula(thisformula),thisdata,identity)
  finalgroups=list()
  for (kgroup in 1:k) {
    outdata=vector()
    for (igroup in 1:nrow(ksamples[[kgroup]])){
      namegroup=ksamples[[kgroup]][igroup,1]
      sampledi=sample(seq(1,length(groups[[igroup,2]]),1),size=ksamples[[kgroup]][igroup,2],replace=FALSE) #Not one
      sampled=groups[[igroup,2]][sampledi]
      groups[[igroup,2]]=groups[[igroup,2]][-sampledi]
      outdata=c(outdata,as.character(sampled))
    }
    finalgroups[[kgroup]]=outdata
  }
  return(finalgroups)
}

##Function that collapses options
collapsecondst=function(x){
	ux=unique(x)
	paste(ux,collapse="|")
}

##Wrapper for ciMean so that it returns a real NA instead of a logical one. Needed for data.table
myCiMean=function(x,...){
	val=ciMean(x,...)[1]
	if(is.na(val)) {
		return(NA_real_)
	}else{
		return(val)
	}
}
  

##################################################



if (file.exists("analysis.RData")) {
  load(file="analysis.RData")
} else {
  dnacon=read.csv("dnacon.csv",header = FALSE)
  #Dict-like array. "Key" == sample, value = dnacon
  dnacons=as.vector(dnacon[,2])
  names(dnacons)=as.vector(dnacon[,1])
  datanofilt=read.csv("results.noscript.csv")
  datanofilt$dnacon=dnacons[as.character(datanofilt$sample)]
  
  data=fread(input='gunzip -c summary_results.csv.gz',header=TRUE,sep=",")

  #Making condition information
  data[,original_condition:=do.call(paste,c(.SD,sep=",")),.SDcols=c(2:(n_sim_params+1))]
  data[,conditionFilt:=do.call(paste,c(.SD,sep=",")),.SDcols=c(2:(n_sim_param_simple+1))]
  data[,conditioncovB:=do.call(paste,c(.SD,sep=",")),.SDcols=c((n_sim_param_simple+2):(n_sim_param_simple+n_sim_param_covB+1))]
  data[,conditionNAB:=do.call(paste,c(.SD,sep=",")),.SDcols=c((n_sim_param_simple+n_sim_param_covB+2):(n_sim_param_simple+n_sim_param_covB+n_sim_param_NAB+1))]
  ##We are using just one PAF condition, so I do not use the generalized option like in the previous ones
  #data[,conditionPAF:=do.call(paste,c(.SD,sep=",")),.SDcols=c((n_sim_param_simple+n_sim_param_covB+n_sim_param_NAB+2):(n_sim_params+1))]
  data[,conditionPAF:=PAFmax_pAF]

  #Adding info on original amount of DNA
  data[,dnacon:=dnacons[Sample]]
  data[,dnagroup:=as.character(dnacon)]
  data[dnacon>=100,dnagroup:=">=100"]

  #Making NAB conditions in numeric
  tomod=colnames(data)[(n_sim_param_simple+n_sim_param_covB+2):(n_sim_param_simple+n_sim_param_covB+n_sim_param_NAB+1)]
  data[,(tomod):=lapply(.SD,sub,pattern="0_",replacement=""),.SDcols=tomod]
  data[,(tomod):=lapply(.SD,as.numeric),.SDcols=tomod]
  data[,(tomod):=lapply(.SD,function(x){ifelse(x==-1,1/0,x)}),.SDcols=tomod]
  
  #New condition with NAB fixed
  data[,condition:=do.call(paste,c(.SD,sep=",")),.SDcols=c(2:(n_sim_params+1))]


  #Selection of filtering options attending to propU taking into consideration all paremeters except PAF and all samples
  ###########################################################################################################
  toworkwith=colnames(data)[2:(n_sim_params+1)]
  data[,`:=`(filtNABcovB_propU_ciMean_bycond=myCiMean(filtNABcovB_propU,conf=.9,na.rm=TRUE)[1]),keyby=condition]
  lowciMean_cond_uniq=data[,.(condition=paste(lapply(.SD,collapsecondst),collapse=","),repcond=paste(lapply(.SD,first),collapse=",")),by=filtNABcovB_propU_ciMean_bycond,.SDcols=toworkwith]

  #data[,`:=`(filtNABcovB_propU_stquantile_bycond=quantile(filtNABcovB_propU,p=.25)),by=condition]
  #stquantile_cond_uniq=data[,.(condition=paste(lapply(.SD,collapsecondst),collapse=","),repcond=paste(lapply(.SD,first),collapse=",")),by=filtNABcovB_propU_stquantile_bycond,.SDcols=toworkwith]
  
  #Selection of filtering options attending to the proportion of common variants that have a PAF<0.05 taking into consideration all paremeters except PAF and all samples
  ###########################################################################################################
  data[,`:=`(filtNABcovB_U_PAF_0.05_ciMean_bycond=myCiMean(filtNABcovB_U_PAF_0.05,conf=.9,na.rm=TRUE)[1]),keyby=condition]
  lowciMean_FPAF0.05_cond_uniq=data[,.(condition=paste(lapply(.SD,collapsecondst),collapse=","),repcond=paste(lapply(.SD,first),collapse=",")),by=filtNABcovB_U_PAF_0.05_ciMean_bycond,.SDcols=toworkwith]
  
  #Max Euclidean distance
  maxE=sqrt(2)
  
  #Selection of filtering options attending to the euclidean distance to the best score 1 1 attending to similarity without taking into consideration PAF and the FPAF at 0.05 (two dimensions)
  ###########################################################################################################
  data[,`:=`(scorefiltNABcovB_FPAF0.05=1-sqrt((1-filtNABcovB_U_PAF_0.05)**2+(1-filtNABcovB_propU)**2)/maxE)]
  data[,`:=`(scorefiltNABcovB_FPAF0.05_ciMean_bycond=myCiMean(scorefiltNABcovB_FPAF0.05,conf=.9,na.rm=TRUE)[1]),keyby=condition]
  #data[,`:=`(scorefiltNABcovB_FPAF0.05=myCiMean(sqrt((1-filtNABcovB_U_PAF_0.05)**2+(1-filtNABcovB_propU)**2),conf=.9,na.rm=TRUE)[1]),keyby=condition]
  lowciMean_Euclidean_cond_uniq=data[,.(condition=paste(lapply(.SD,collapsecondst),collapse=","),repcond=paste(lapply(.SD,first),collapse=",")),by=scorefiltNABcovB_FPAF0.05_ciMean_bycond,.SDcols=toworkwith]
  #data[,`:=`(filtNABcovB_U_PAF_0.05_stquantile_bycond=quantile(filtNABcovB_U_PAF_0.05,p=.25)),by=condition]
  #stquantile_FPAF0.05_cond_uniq=data[,.(condition=paste(lapply(.SD,collapsecondst),collapse=","),repcond=paste(lapply(.SD,first),collapse=",")),by=filtNABcovB_U_PAF_0.05_stquantile_bycond,.SDcols=toworkwith]
  
  
  #Selection of filtering options attending to the euclidean distance to the best score 1 1 attending to similarity eliminating variants at the PAF 0.25 level and co-opimizing the FPAF at 0.05 (two dimensions)
  ###########################################################################################################
  data[,`:=`(scorefiltNABcovB_PAFFPAF0.05=1-sqrt((1-filtNABcovBPAF_U_PAF_0.05)**2+(1-filtNABcovBPAF_propU)**2)/maxE)]
  data[,`:=`(scorefiltNABcovB_PAFFPAF0.05_ciMean_bycond=myCiMean(scorefiltNABcovB_PAFFPAF0.05,conf=.9,na.rm=TRUE)[1]),keyby=condition]
  #data[,`:=`(scorefiltNABcovB_FPAF0.05=myCiMean(sqrt((1-filtNABcovB_U_PAF_0.05)**2+(1-filtNABcovB_propU)**2),conf=.9,na.rm=TRUE)[1]),keyby=condition]
  lowciMean_EuclideanPAF_cond_uniq=data[conditionPAF==0.25,.(condition=paste(lapply(.SD,collapsecondst),collapse=","),repcond=paste(lapply(.SD,first),collapse=",")),by=scorefiltNABcovB_PAFFPAF0.05_ciMean_bycond,.SDcols=toworkwith]
  #data[,`:=`(filtNABcovB_U_PAF_0.05_stquantile_bycond=quantile(filtNABcovB_U_PAF_0.05,p=.25)),by=condition]
  #stquantile_FPAF0.05_cond_uniq=data[,.(condition=paste(lapply(.SD,collapsecondst),collapse=","),repcond=paste(lapply(.SD,first),collapse=",")),by=filtNABcovB_U_PAF_0.05_stquantile_bycond,.SDcols=toworkwith]
  
  #Selection of filtering options attending to the euclidean distance to the best score 1 1 attending to similarity eliminating variants at the PAF 0.25 level and co-opimizing the MeanPAF (two dimensions)
  ###########################################################################################################
  data[,`:=`(scorefiltNABcovBPAF_MPAF=1-sqrt((filtNABcovBPAF_U_PAF_Mean)**2+(1-filtNABcovBPAF_propU)**2)/maxE)]
  data[,`:=`(scorefiltNABcovBPAF_MPAF_ciMean_bycond=myCiMean(scorefiltNABcovBPAF_MPAF,conf=.9,na.rm=TRUE)[1]),keyby=condition]
  lowciMean_EuclideanMPAF_cond_uniq=data[conditionPAF==0.25,.(condition=paste(lapply(.SD,collapsecondst),collapse=","),repcond=paste(lapply(.SD,first),collapse=",")),by=scorefiltNABcovBPAF_MPAF_ciMean_bycond,.SDcols=toworkwith]
  #data[,`:=`(filtNABcovB_U_PAF_0.05_stquantile_bycond=quantile(filtNABcovB_U_PAF_0.05,p=.25)),by=condition]
  #stquantile_FPAF0.05_cond_uniq=data[,.(condition=paste(lapply(.SD,collapsecondst),collapse=","),repcond=paste(lapply(.SD,first),collapse=",")),by=filtNABcovB_U_PAF_0.05_stquantile_bycond,.SDcols=toworkwith]

  #Selection of NABfiltering options attending to the proportion of common variants that have a PAF<0.05 
  ###########################################################################################################
  data[,`:=`(filtNABcovB_U_PAF_0.05_ciMean_bycondNAB=myCiMean(filtNABcovB_U_PAF_0.05,conf=.9,na.rm=TRUE)[1]),by=conditionNAB]

#  data[,`:=`(filtNABcovB_U_PAF_0.05_stquantile_bycond=quantile(filtNABcovB_U_PAF_0.05,p=.25)),by=condition]
#  stquantile_FPAF0.05_cond_uniq=data[,.(condition=paste(lapply(.SD,collapsecondst),collapse=","),repcond=paste(lapply(.SD,first),collapse=",")),by=filtNABcovB_U_PAF_0.05_stquantile_bycond,.SDcols=toworkwith]

  #Selection of filtering options attending to propU taking into consideration all paremeters and all samples including PAF filter
  ###########################################################################################################
  data[,`:=`(filtNABcovBPAF_propU_ciMean_bycond=myCiMean(filtNABcovBPAF_propU,conf=.9,na.rm=TRUE)[1]),by=condition]
  lowciMeanPAF_cond_uniq=data[,.(condition=paste(lapply(.SD,collapsecondst),collapse=","),repcond=paste(lapply(.SD,first),collapse=",")),by=filtNABcovBPAF_propU_ciMean_bycond,.SDcols=toworkwith]

  #data[,`:=`(filtNABcovBPAF_propU_stquantile_bycond=quantile(filtNABcovBPAF_propU,p=.25)),by=condition]
  #stquantilePAF_cond_uniq=data[,.(condition=paste(lapply(.SD,collapsecondst),collapse=","),repcond=paste(lapply(.SD,first),collapse=",")),by=filtNABcovBPAF_propU_stquantile_bycond,.SDcols=toworkwith]
  
  #Selection of filtering options attending to propU taking into consideration all paremeters and all samples, selecting only filters with PAF>0.05
  #####################################################################################################################################
  lowciMeanPAF0.05_cond_uniq=data[conditionPAF==0.05,.(condition=paste(lapply(.SD,collapsecondst),collapse=","),repcond=paste(lapply(.SD,first),collapse=",")),by=filtNABcovBPAF_propU_ciMean_bycond,.SDcols=toworkwith]

#  #k-fold cross-validation
#  ########################
#  datakgroups=as.data.frame(dnacon)
#  colnames(datakgroups)=c("Sample","DNAcon")
#  datakgroups$DNAcon=as.factor(sapply(datakgroups$DNAcon,function(x){if(x>=100) ">=100" else x}))
#  listks=makekgroups(datakgroups)
#  setindex(data,Sample) 
#
#  #NoPAF
#  resultscrossvalidation=mclapply(listks,function (mydataks){
#	resultski=data[(mydataks),.(lowciMean_test=ciMean(filtNABcovB_propU,conf=.9)[1]),keyby=condition,on="Sample"]
#	resultski[,lowciMean_training:=data[!(mydataks),ciMean(filtNABcovB_propU,conf=.9)[1],by=condition,on="Sample"][[2]]]
#	resultski[,stquantile_test:=data[(mydataks),quantile(filtNABcovB_propU,p=.25)[1],by=condition,on="Sample"][[2]]]
#	resultski[,stquantile_training:=data[!(mydataks),quantile(filtNABcovB_propU,p=.25)[1],by=condition,on="Sample"][[2]]]
#	
#    return(resultski)
#  }, mc.cores = n_cores)
#  
#  resultscrossvalidation=bind_rows(resultscrossvalidation,.id="k")
#  setDT(resultscrossvalidation)
#  resultscrossvalidation=melt(resultscrossvalidation,measure=patterns("^lowciMean","^stquantile"),value.name=c("lowciMean","stquantile"),variable.name=c("type"))  
#  levels(resultscrossvalidation$type)=list(test="1",training="2")
#
#  resultscrossvalidation=resultscrossvalidation[,c(lapply(.SD,mean),lapply(.SD,sd)),keyby=c("condition","type"),.SDcols=c("lowciMean","stquantile")]
#  setnames(resultscrossvalidation,3,"mean_lowciMean")
#  setnames(resultscrossvalidation,4,"mean_stquantile")
#  setnames(resultscrossvalidation,5,"sd_lowciMean")
#  setnames(resultscrossvalidation,6,"sd_stquantile")
#  resultscrossvalidation=merge(resultscrossvalidation,data[,lapply(.SD,first),keyby=condition,.SDcols=toworkwith])
#
#  #PAF
#  resultscrossvalidationPAF=mclapply(listks,function (mydataks){
#	resultski=data[(mydataks),.(lowciMeanPAF_test=ciMean(filtNABcovBPAF_propU,conf=.9)[1]),keyby=condition,on="Sample"]
#	resultski[,lowciMeanPAF_test:=data[mydataks,ciMean(filtNABcovBPAF_propU,conf=.9)[1],by=condition,on="Sample"][[2]]]
#	resultski[,lowciMeanPAF_training:=data[!(mydataks),ciMean(filtNABcovBPAF_propU,conf=.9)[1],by=condition,on="Sample"][[2]]]
#	resultski[,stquantilePAF_test:=data[(mydataks),quantile(filtNABcovBPAF_propU,p=.25)[1],by=condition,on="Sample"][[2]]]
#	resultski[,stquantilePAF_training:=data[!(mydataks),quantile(filtNABcovBPAF_propU,p=.25)[1],by=condition,on="Sample"][[2]]]
#	
#    return(resultski)
#  }, mc.cores = n_cores)
#  
#  resultscrossvalidationPAF=bind_rows(resultscrossvalidationPAF,.id="k")
#  setDT(resultscrossvalidationPAF)
#  resultscrossvalidationPAF=melt(resultscrossvalidationPAF,measure=patterns("^lowciMean","^stquantile"),value.name=c("lowciMeanPAF","stquantilePAF"),variable.name=c("type"))  
#  levels(resultscrossvalidationPAF$type)=list(test="1",training="2")
#
#  resultscrossvalidationPAF=resultscrossvalidationPAF[,c(lapply(.SD,mean),lapply(.SD,sd)),keyby=c("condition","type"),.SDcols=c("lowciMeanPAF","stquantilePAF")]
#  setnames(resultscrossvalidationPAF,3,"mean_lowciMeanPAF")
#  setnames(resultscrossvalidationPAF,4,"mean_stquantilePAF")
#  setnames(resultscrossvalidationPAF,5,"sd_lowciMeanPAF")
#  setnames(resultscrossvalidationPAF,6,"sd_stquantilePAF")
#  resultscrossvalidationPAF=merge(resultscrossvalidationPAF,data[,lapply(.SD,first),keyby=condition,.SDcols=toworkwith])
 
  save.image("analysis.RData")
  
}

#Fixes
ordered_samples=data[,.(median_filtNABcovB_propU=median(filtNABcovB_propU)),by=Sample][order(median_filtNABcovB_propU),Sample]
ordered_dnas=c("20","40","60","80",">=100")
data[,`:=`(dnagroupF=factor(dnagroup,ordered=TRUE,levels=ordered_dnas))]
names(data)[which(names(data)=="filtNABcovBPAF_#U")]="filtNABcovBPAF_TU"

##Best condition PAF<0.05
#########################

best_cond_PAF0.05=lowciMeanPAF0.05_cond_uniq[order(-filtNABcovBPAF_propU_ciMean_bycond)][1,.(repcond)]

##Best condition Euclidean co-optimization
##########################################

best_cond_euclidean=lowciMean_Euclidean_cond_uniq[order(-scorefiltNABcovB_FPAF0.05_ciMean_bycond)][1,.(repcond)]

##Best condition EuclideanPAF co-optimization
#############################################
best_cond_euclideanPAF=lowciMean_EuclideanPAF_cond_uniq[order(-scorefiltNABcovB_PAFFPAF0.05_ciMean_bycond)][1,.(repcond)]

##Best condition EuclideanMeanPAF co-optimization
#############################################
best_cond_euclideanMPAF=lowciMean_EuclideanMPAF_cond_uniq[order(-scorefiltNABcovBPAF_MPAF_ciMean_bycond)][1,.(repcond)]

##Small summary tables
write.table(file="EPAF_summary.csv",data[best_cond_euclideanPAF,.(Sample,filtNABcovBPAF_propU,filtNABcovBPAF_TU,filtNABcovBPAF_U_PAF_Mean)],quote=FALSE,sep="\t")
write.table(file="PAF0.05_summary.csv",data[best_cond_PAF0.05,.(Sample,filtNABcovBPAF_propU,filtNABcovBPAF_TU,filtNABcovBPAF_U_PAF_Mean)],quote=FALSE,sep="\t")

##Plots
##########################################


##Best condition PAF<0.05
#########################
data[,`:=`(SampleF=factor(Sample,ordered=TRUE,levels=ordered_samples))]


plot_sim_bestPAF0.05_text=ggplot(data=data,aes(y=filtNABcovBPAF_propU*100,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_PAF0.05],aes(ymin=filtNABcovBPAF_propU*100,ymax=filtNABcovBPAF_propU*100,linetype="solid"),color="#c51b7d",size=0.35)+geom_point(data=datanofilt,aes(x=sample,y=sim*100,shape=""),size=2.5,stroke=1)+scale_y_continuous(name="Similarity (%)")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_shape_manual(values=c(4),labels="Default\npipeline",drop=FALSE,name="")+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))

if (!file.exists("plot_sim_bestPAF0.05_text.pdf"))
{
  save_plot("plot_sim_bestPAF0.05_text.pdf",plot_sim_bestPAF0.05_text,base_height=6,base_aspect_ratio=1.6)
}


plot_vars_bestPAF0.05_text=ggplot(data=data,aes(y=filtNABcovBPAF_TU,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_PAF0.05],aes(ymin=filtNABcovBPAF_TU,ymax=filtNABcovBPAF_TU,linetype="solid"),color="#c51b7d",size=0.35)+geom_point(data=datanofilt,aes(x=sample,y=X.AN+X.BN+X.comun,shape=""),size=2.5,stroke=1)+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_shape_manual(values=c(4),labels="Default\npipeline",drop=FALSE,name="")+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))


if (!file.exists("plot_vars_bestPAF0.05_text.pdf"))
{
  save_plot("plot_vars_bestPAF0.05_text.pdf",plot_vars_bestPAF0.05_text,base_height=6,base_aspect_ratio=1.6)
}

plot_vars_log_bestPAF0.05_text=ggplot(data=data,aes(y=filtNABcovBPAF_TU,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_PAF0.05],aes(ymin=filtNABcovBPAF_TU,ymax=filtNABcovBPAF_TU,linetype="solid"),color="#c51b7d",size=0.35)+geom_point(data=datanofilt,aes(x=sample,y=X.AN+X.BN+X.comun,shape=""),size=2.5,stroke=1)+scale_y_log10(name="Log10 number of variants")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_shape_manual(values=c(4),labels="Default\npipeline",drop=FALSE,name="")+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))


if (!file.exists("plot_vars_log_bestPAF0.05_text.pdf"))
{
  save_plot("plot_vars_log_bestPAF0.05_text.pdf",plot_vars_log_bestPAF0.05_text,base_height=6,base_aspect_ratio=1.6)
}

#plot_fpaf_bestPAF0.05_text=ggplot(data=data,aes(y=filtNABcovB_U_PAF_0.05*100,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_PAF0.05],aes(ymin=filtNABcovB_U_PAF_0.05*100,ymax=filtNABcovB_U_PAF_0.05*100,linetype="solid"),color="#c51b7d",size=0.35)+geom_point(data=datanofilt,aes(x=sample,y=fpaf*100,shape=""),size=2.5,stroke=1)+scale_y_continuous(name="Common samples with PAF<=0.05 (%)")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_shape_manual(values=c(4),labels="Default\npipeline",drop=FALSE,name="")+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))

##I do not have the default yet
plot_fpaf_bestPAF0.05_text=ggplot(data=data,aes(y=filtNABcovB_U_PAF_0.05*100,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_PAF0.05],aes(ymin=filtNABcovB_U_PAF_0.05*100,ymax=filtNABcovB_U_PAF_0.05*100,linetype="solid"),color="#c51b7d",size=0.35)+scale_y_continuous(name="Common samples with PAF<=0.05 (%)")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))

if (!file.exists("plot_fpaf_bestPAF0.05_text.pdf"))
{
  save_plot("plot_fpaf_bestPAF0.05_text.pdf",plot_fpaf_bestPAF0.05_text,base_height=6,base_aspect_ratio=1.6)
}

#MeanPAF
plot_mpaf_bestPAF0.05_text=ggplot(data=data,aes(y=filtNABcovBPAF_U_PAF_Mean,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_PAF0.05],aes(ymin=filtNABcovBPAF_U_PAF_Mean,ymax=filtNABcovBPAF_U_PAF_Mean,linetype="solid"),color="#c51b7d",size=0.35)+scale_y_continuous(name="Mean PAF")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))

if (!file.exists("plot_mpaf_bestPAF0.05_text.pdf"))
{
  save_plot("plot_mpaf_bestPAF0.05_text.pdf",plot_mpaf_bestPAF0.05_text,base_height=6,base_aspect_ratio=1.6)
}


##Best co-optimization
######################

ordered_samples=data[,.(median_scorefiltNABcovB_FPAF0.05=median(scorefiltNABcovB_FPAF0.05)),by=Sample][order(median_scorefiltNABcovB_FPAF0.05),Sample]
data[,`:=`(SampleF=factor(Sample,ordered=TRUE,levels=ordered_samples))]
names(data)[which(names(data)=="filtNABcovB_#U")]="filtNABcovB_TU"

##I do not have the default for the score yet
plot_score_bestEuclidean_text=ggplot(data=data,aes(y=scorefiltNABcovB_FPAF0.05,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_euclidean],aes(ymin=scorefiltNABcovB_FPAF0.05,ymax=scorefiltNABcovB_FPAF0.05,linetype="solid"),color="#c51b7d",size=0.35)+scale_y_continuous(name="Score (complement of prop. max. theoretical distance)")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))


##plot_score_bestEuclidean_text=ggplot(data=data,aes(y=scorefiltNABcovB_FPAF0.05,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_euclidean],aes(ymin=scorefiltNABcovB_FPAF0.05,ymax=scorefiltNABcovB_FPAF0.05,linetype="solid"),color="#c51b7d",size=0.35)+geom_point(data=datanofilt,aes(x=sample,y=score*100,shape=""),size=2.5,stroke=1)+scale_y_continuous(name="Similarity (%)")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_shape_manual(values=c(4),labels="Default\npipeline",drop=FALSE,name="")+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))

if (!file.exists("plot_score_bestEuclidean_text.pdf"))
{
  save_plot("plot_score_bestEuclidean_text.pdf",plot_score_bestEuclidean_text,base_height=6,base_aspect_ratio=1.6)
}

plot_sim_bestEuclidean_text=ggplot(data=data,aes(y=filtNABcovB_propU*100,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_euclidean],aes(ymin=filtNABcovB_propU*100,ymax=filtNABcovB_propU*100,linetype="solid"),color="#c51b7d",size=0.35)+geom_point(data=datanofilt,aes(x=sample,y=sim*100,shape=""),size=2.5,stroke=1)+scale_y_continuous(name="Similarity (%)")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_shape_manual(values=c(4),labels="Default\npipeline",drop=FALSE,name="")+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))

if (!file.exists("plot_sim_bestEuclidean_text.pdf"))
{
  save_plot("plot_sim_bestEuclidean_text.pdf",plot_sim_bestEuclidean_text,base_height=6,base_aspect_ratio=1.6)
}


plot_vars_bestEuclidean_text=ggplot(data=data,aes(y=filtNABcovB_TU,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_euclidean],aes(ymin=filtNABcovB_TU,ymax=filtNABcovB_TU,linetype="solid"),color="#c51b7d",size=0.35)+geom_point(data=datanofilt,aes(x=sample,y=X.AN+X.BN+X.comun,shape=""),size=2.5,stroke=1)+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_shape_manual(values=c(4),labels="Default\npipeline",drop=FALSE,name="")+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))


if (!file.exists("plot_vars_bestEuclidean_text.pdf"))
{
  save_plot("plot_vars_bestEuclidean_text.pdf",plot_vars_bestEuclidean_text,base_height=6,base_aspect_ratio=1.6)
}

plot_vars_log_bestEuclidean_text=ggplot(data=data,aes(y=filtNABcovB_TU,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_euclidean],aes(ymin=filtNABcovB_TU,ymax=filtNABcovB_TU,linetype="solid"),color="#c51b7d",size=0.35)+geom_point(data=datanofilt,aes(x=sample,y=X.AN+X.BN+X.comun,shape=""),size=2.5,stroke=1)+scale_y_log10(name="Log10 number of variants")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_shape_manual(values=c(4),labels="Default\npipeline",drop=FALSE,name="")+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))


if (!file.exists("plot_vars_log_bestEuclidean_text.pdf"))
{
  save_plot("plot_vars_log_bestEuclidean_text.pdf",plot_vars_log_bestEuclidean_text,base_height=6,base_aspect_ratio=1.6)
}

#plot_fpaf_bestEuclidean_text=ggplot(data=data,aes(y=filtNABcovB_U_PAF_0.05*100,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_euclidean],aes(ymin=filtNABcovB_U_PAF_0.05*100,ymax=filtNABcovB_U_PAF_0.05*100,linetype="solid"),color="#c51b7d",size=0.35)+geom_point(data=datanofilt,aes(x=sample,y=fpaf*100,shape=""),size=2.5,stroke=1)+scale_y_continuous(name="Common samples with PAF<=0.05 (%)")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_shape_manual(values=c(4),labels="Default\npipeline",drop=FALSE,name="")+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))

##I do not have the default yet
plot_fpaf_bestEuclidean_text=ggplot(data=data,aes(y=filtNABcovB_U_PAF_0.05*100,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_euclidean],aes(ymin=filtNABcovB_U_PAF_0.05*100,ymax=filtNABcovB_U_PAF_0.05*100,linetype="solid"),color="#c51b7d",size=0.35)+scale_y_continuous(name="Common samples with PAF<=0.05 (%)")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))

if (!file.exists("plot_fpaf_bestEuclidean_text.pdf"))
{
  save_plot("plot_fpaf_bestEuclidean_text.pdf",plot_fpaf_bestEuclidean_text,base_height=6,base_aspect_ratio=1.6)
}

#MeanPAF
plot_mpaf_bestEuclidean_text=ggplot(data=data,aes(y=filtNABcovB_U_PAF_Mean,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_euclidean],aes(ymin=filtNABcovB_U_PAF_Mean,ymax=filtNABcovB_U_PAF_Mean,linetype="solid"),color="#c51b7d",size=0.35)+scale_y_continuous(name="Mean PAF")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))

if (!file.exists("plot_mpaf_bestEuclidean_text.pdf"))
{
  save_plot("plot_mpaf_bestEuclidean_text.pdf",plot_mpaf_bestEuclidean_text,base_height=6,base_aspect_ratio=1.6)
}

##Best co-optimizationPAF
######################

ordered_samples=data[,.(median_scorefiltNABcovB_PAFFPAF0.05=median(scorefiltNABcovB_PAFFPAF0.05)),by=Sample][order(median_scorefiltNABcovB_PAFFPAF0.05),Sample]
data[,`:=`(SampleF=factor(Sample,ordered=TRUE,levels=ordered_samples))]
names(data)[which(names(data)=="filtNABcovBPAF_#U")]="filtNABcovBPAF_TU"

##I do not have the default for the score yet
plot_score_bestEuclideanPAF_text=ggplot(data=data,aes(y=scorefiltNABcovB_PAFFPAF0.05,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_euclideanPAF],aes(ymin=scorefiltNABcovB_PAFFPAF0.05,ymax=scorefiltNABcovB_PAFFPAF0.05,linetype="solid"),color="#c51b7d",size=0.35)+scale_y_continuous(name="Score (complement of prop. max. theoretical distance)")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))


##plot_score_bestEuclidean_text=ggplot(data=data,aes(y=scorefiltNABcovB_FPAF0.05,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_euclidean],aes(ymin=scorefiltNABcovB_FPAF0.05,ymax=scorefiltNABcovB_FPAF0.05,linetype="solid"),color="#c51b7d",size=0.35)+geom_point(data=datanofilt,aes(x=sample,y=score*100,shape=""),size=2.5,stroke=1)+scale_y_continuous(name="Similarity (%)")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_shape_manual(values=c(4),labels="Default\npipeline",drop=FALSE,name="")+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))

if (!file.exists("plot_score_bestEuclideanPAF_text.pdf"))
{
  save_plot("plot_score_bestEuclideanPAF_text.pdf",plot_score_bestEuclideanPAF_text,base_height=6,base_aspect_ratio=1.6)
}

plot_sim_bestEuclideanPAF_text=ggplot(data=data,aes(y=filtNABcovBPAF_propU*100,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_euclideanPAF],aes(ymin=filtNABcovBPAF_propU*100,ymax=filtNABcovBPAF_propU*100,linetype="solid"),color="#c51b7d",size=0.35)+geom_point(data=datanofilt,aes(x=sample,y=sim*100,shape=""),size=2.5,stroke=1)+scale_y_continuous(name="Similarity (%)")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_shape_manual(values=c(4),labels="Default\npipeline",drop=FALSE,name="")+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))

if (!file.exists("plot_sim_bestEuclideanPAF_text.pdf"))
{
  save_plot("plot_sim_bestEuclideanPAF_text.pdf",plot_sim_bestEuclideanPAF_text,base_height=6,base_aspect_ratio=1.6)
}


plot_vars_bestEuclideanPAF_text=ggplot(data=data,aes(y=filtNABcovBPAF_TU,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_euclideanPAF],aes(ymin=filtNABcovBPAF_TU,ymax=filtNABcovBPAF_TU,linetype="solid"),color="#c51b7d",size=0.35)+geom_point(data=datanofilt,aes(x=sample,y=X.AN+X.BN+X.comun,shape=""),size=2.5,stroke=1)+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_shape_manual(values=c(4),labels="Default\npipeline",drop=FALSE,name="")+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))


if (!file.exists("plot_vars_bestEuclideanPAF_text.pdf"))
{
  save_plot("plot_vars_bestEuclideanPAF_text.pdf",plot_vars_bestEuclideanPAF_text,base_height=6,base_aspect_ratio=1.6)
}

plot_vars_log_bestEuclideanPAF_text=ggplot(data=data,aes(y=filtNABcovBPAF_TU,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_euclideanPAF],aes(ymin=filtNABcovBPAF_TU,ymax=filtNABcovBPAF_TU,linetype="solid"),color="#c51b7d",size=0.35)+geom_point(data=datanofilt,aes(x=sample,y=X.AN+X.BN+X.comun,shape=""),size=2.5,stroke=1)+scale_y_log10(name="Log10 number of variants")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_shape_manual(values=c(4),labels="Default\npipeline",drop=FALSE,name="")+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))


if (!file.exists("plot_vars_log_bestEuclideanPAF_text.pdf"))
{
  save_plot("plot_vars_log_bestEuclideanPAF_text.pdf",plot_vars_log_bestEuclideanPAF_text,base_height=6,base_aspect_ratio=1.6)
}

#plot_fpaf_bestEuclidean_text=ggplot(data=data,aes(y=filtNABcovB_U_PAF_0.05*100,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_euclidean],aes(ymin=filtNABcovB_U_PAF_0.05*100,ymax=filtNABcovB_U_PAF_0.05*100,linetype="solid"),color="#c51b7d",size=0.35)+geom_point(data=datanofilt,aes(x=sample,y=fpaf*100,shape=""),size=2.5,stroke=1)+scale_y_continuous(name="Common samples with PAF<=0.05 (%)")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_shape_manual(values=c(4),labels="Default\npipeline",drop=FALSE,name="")+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))

##I do not have the default yet
plot_fpaf_bestEuclideanPAF_text=ggplot(data=data,aes(y=filtNABcovB_U_PAF_0.05*100,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_euclideanPAF],aes(ymin=filtNABcovB_U_PAF_0.05*100,ymax=filtNABcovB_U_PAF_0.05*100,linetype="solid"),color="#c51b7d",size=0.35)+scale_y_continuous(name="Common samples with PAF<=0.05 (%)")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))

if (!file.exists("plot_fpaf_bestEuclideanPAF_text.pdf"))
{
  save_plot("plot_fpaf_bestEuclideanPAF_text.pdf",plot_fpaf_bestEuclideanPAF_text,base_height=6,base_aspect_ratio=1.6)
}

#MeanPAF
plot_mpaf_bestEuclideanPAF_text=ggplot(data=data,aes(y=filtNABcovBPAF_U_PAF_Mean,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_euclideanPAF],aes(ymin=filtNABcovBPAF_U_PAF_Mean,ymax=filtNABcovBPAF_U_PAF_Mean,linetype="solid"),color="#c51b7d",size=0.35)+scale_y_continuous(name="Mean PAF")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))

if (!file.exists("plot_mpaf_bestEuclideanPAF_text.pdf"))
{
  save_plot("plot_mpaf_bestEuclideanPAF_text.pdf",plot_mpaf_bestEuclideanPAF_text,base_height=6,base_aspect_ratio=1.6)
}

##Best co-optimization meanPAF
##############################

ordered_samples=data[,.(median_scorefiltNABcovBPAF_MPAF=median(scorefiltNABcovBPAF_MPAF)),by=Sample][order(median_scorefiltNABcovBPAF_MPAF),Sample]
data[,`:=`(SampleF=factor(Sample,ordered=TRUE,levels=ordered_samples))]
names(data)[which(names(data)=="filtNABcovBPAF_#U")]="filtNABcovBPAF_TU"

##I do not have the default for the score yet
plot_score_bestEuclideanMPAF_text=ggplot(data=data,aes(y=scorefiltNABcovBPAF_MPAF,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_euclideanMPAF],aes(ymin=scorefiltNABcovBPAF_MPAF,ymax=scorefiltNABcovBPAF_MPAF,linetype="solid"),color="#c51b7d",size=0.35)+scale_y_continuous(name="Score (complement of prop. max. theoretical distance)")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))


##plot_score_bestEuclidean_text=ggplot(data=data,aes(y=scorefiltNABcovB_FPAF0.05,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_euclidean],aes(ymin=scorefiltNABcovB_FPAF0.05,ymax=scorefiltNABcovB_FPAF0.05,linetype="solid"),color="#c51b7d",size=0.35)+geom_point(data=datanofilt,aes(x=sample,y=score*100,shape=""),size=2.5,stroke=1)+scale_y_continuous(name="Similarity (%)")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_shape_manual(values=c(4),labels="Default\npipeline",drop=FALSE,name="")+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))

if (!file.exists("plot_score_bestEuclideanMPAF_text.pdf"))
{
  save_plot("plot_score_bestEuclideanMPAF_text.pdf",plot_score_bestEuclideanMPAF_text,base_height=6,base_aspect_ratio=1.6)
}

plot_sim_bestEuclideanMPAF_text=ggplot(data=data,aes(y=filtNABcovBPAF_propU*100,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_euclideanMPAF],aes(ymin=filtNABcovBPAF_propU*100,ymax=filtNABcovBPAF_propU*100,linetype="solid"),color="#c51b7d",size=0.35)+geom_point(data=datanofilt,aes(x=sample,y=sim*100,shape=""),size=2.5,stroke=1)+scale_y_continuous(name="Similarity (%)")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_shape_manual(values=c(4),labels="Default\npipeline",drop=FALSE,name="")+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))

if (!file.exists("plot_sim_bestEuclideanMPAF_text.pdf"))
{
  save_plot("plot_sim_bestEuclideanMPAF_text.pdf",plot_sim_bestEuclideanMPAF_text,base_height=6,base_aspect_ratio=1.6)
}


plot_vars_bestEuclideanMPAF_text=ggplot(data=data,aes(y=filtNABcovBPAF_TU,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_euclideanMPAF],aes(ymin=filtNABcovBPAF_TU,ymax=filtNABcovBPAF_TU,linetype="solid"),color="#c51b7d",size=0.35)+geom_point(data=datanofilt,aes(x=sample,y=X.AN+X.BN+X.comun,shape=""),size=2.5,stroke=1)+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_shape_manual(values=c(4),labels="Default\npipeline",drop=FALSE,name="")+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))


if (!file.exists("plot_vars_bestEuclideanMPAF_text.pdf"))
{
  save_plot("plot_vars_bestEuclideanMPAF_text.pdf",plot_vars_bestEuclideanMPAF_text,base_height=6,base_aspect_ratio=1.6)
}

plot_vars_log_bestEuclideanMPAF_text=ggplot(data=data,aes(y=filtNABcovBPAF_TU,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_euclideanMPAF],aes(ymin=filtNABcovBPAF_TU,ymax=filtNABcovBPAF_TU,linetype="solid"),color="#c51b7d",size=0.35)+geom_point(data=datanofilt,aes(x=sample,y=X.AN+X.BN+X.comun,shape=""),size=2.5,stroke=1)+scale_y_log10(name="Log10 number of variants")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_shape_manual(values=c(4),labels="Default\npipeline",drop=FALSE,name="")+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))


if (!file.exists("plot_vars_log_bestEuclideanMPAF_text.pdf"))
{
  save_plot("plot_vars_log_bestEuclideanMPAF_text.pdf",plot_vars_log_bestEuclideanMPAF_text,base_height=6,base_aspect_ratio=1.6)
}

#plot_fpaf_bestEuclidean_text=ggplot(data=data,aes(y=filtNABcovB_U_PAF_0.05*100,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_euclidean],aes(ymin=filtNABcovB_U_PAF_0.05*100,ymax=filtNABcovB_U_PAF_0.05*100,linetype="solid"),color="#c51b7d",size=0.35)+geom_point(data=datanofilt,aes(x=sample,y=fpaf*100,shape=""),size=2.5,stroke=1)+scale_y_continuous(name="Common samples with PAF<=0.05 (%)")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_shape_manual(values=c(4),labels="Default\npipeline",drop=FALSE,name="")+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))

##I do not have the default yet
plot_fpaf_bestEuclideanMPAF_text=ggplot(data=data,aes(y=filtNABcovB_U_PAF_0.05*100,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_euclideanMPAF],aes(ymin=filtNABcovB_U_PAF_0.05*100,ymax=filtNABcovB_U_PAF_0.05*100,linetype="solid"),color="#c51b7d",size=0.35)+scale_y_continuous(name="Common samples with PAF<=0.05 (%)")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))

if (!file.exists("plot_fpaf_bestEuclideanMPAF_text.pdf"))
{
  save_plot("plot_fpaf_bestEuclideanMPAF_text.pdf",plot_fpaf_bestEuclideanMPAF_text,base_height=6,base_aspect_ratio=1.6)
}

#MeanPAF
plot_mpaf_bestEuclideanMPAF_text=ggplot(data=data,aes(y=filtNABcovBPAF_U_PAF_Mean,x=SampleF))+geom_violin(aes(fill=dnagroupF,color=dnagroupF))+geom_crossbar(data=data[best_cond_euclideanMPAF],aes(ymin=filtNABcovBPAF_U_PAF_Mean,ymax=filtNABcovBPAF_U_PAF_Mean,linetype="solid"),color="#c51b7d",size=0.35)+scale_y_continuous(name="Mean PAF")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+scale_color_brewer(guide=FALSE, type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))

if (!file.exists("plot_mpaf_bestEuclideanMPAF_text.pdf"))
{
  save_plot("plot_mpaf_bestEuclideanMPAF_text.pdf",plot_mpaf_bestEuclideanMPAF_text,base_height=6,base_aspect_ratio=1.6)
}
#
##Let's get the 10 best unique conditions
#########################################
#conditions_ciMean=head(lowciMean_cond_uniq,10)$repcond
#names(conditions_ciMean)=head(lowciMean_cond_uniq,10)$condition
#conditions_stquant=head(stquantile_cond_uniq,10)$repcond
#names(conditions_stquant)=head(stquantile_cond_uniq,10)$condition
#
#conditions_ciMean=factor(conditions_ciMean,levels=conditions_ciMean)
#conditions_stquant=factor(conditions_stquant,levels=conditions_stquant)
#conditions_ciMean %in% conditions_stquant ##Are they the same?
#
##Combined: Plots without scaling
##################################
#
#data$Sample=with(data,reorder(Sample,filtNABcovB_propU,FUN = function(x){return(median(x))}))
#
#plot_sim_ciMean=ggplot(data=data,aes(y=filtNABcovB_propU,x=Sample))+geom_violin(aes(fill=dnagroup))+geom_point(data=data[data$condition %in% conditions_ciMean,],aes(color=factor(condition,levels=conditions_ciMean)),size=2,stat="summary",fun.y="mean")+geom_point(data=datanofilt,aes(x=sample,y=sim,shape=""),size=2.5,stroke=1.5)+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_brewer(name="Filter",type = "qual",palette = 3,labels=names(conditions_ciMean))+labs(title="Somatic similarity. CI90 Mean selection")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_shape_manual(values=c(4),labels="No optimization",drop=FALSE,name="")+guides(colour=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))
#
#if (!file.exists("sim_ciMean.png"))
#{
#  save_plot("sim_ciMean.png",plot_sim_ciMean,base_height=15,base_aspect_ratio=1.5)
#}
#
###For a text
#
#plot_sim_ciMean_text=ggplot(data=data,aes(y=filtNABcovB_propU*100,x=Sample))+geom_violin(aes(fill=dnagroup))+geom_crossbar(data=data[data$condition == conditions_ciMean[1],],aes(ymin=filtNABcovB_propU*100,ymax=filtNABcovB_propU*100,linetype="solid"),color="#c51b7d",size=0.35)+geom_point(data=datanofilt,aes(x=sample,y=sim*100,shape=""),size=2.5,stroke=1)+scale_y_continuous(name="Similarity (%)")+scale_x_discrete(name="Technical replicate")+labs(title="Pipeline optimization across 28 technical replicates")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_shape_manual(values=c(4),labels="Default\npipeline",drop=FALSE,name="")+scale_linetype_manual(values="solid",label="Optimized\npipeline",name="")+guides(linetype=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))
#
#if (!file.exists("sim_ciMean_text.pdf"))
#{
#  save_plot("sim_ciMean_text.pdf",plot_sim_ciMean_text,base_height=6,base_aspect_ratio=1.6)
#}
#
#plot_sim_stquant=ggplot(data=data,aes(y=filtNABcovB_propU,x=Sample))+geom_violin(aes(fill=dnagroup))+geom_point(data=data[data$condition %in% conditions_stquant,],aes(color=factor(condition,levels=conditions_stquant)),size=2,stat="summary",fun.y="mean")+geom_point(data=datanofilt,aes(x=sample,y=sim,shape=""),size=2.5,stroke=1.5)+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_brewer(name="Filter",type="qual",palette=3,labels=names(conditions_stquant))+labs(title="Somatic similarity. Q25 selection")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_shape_manual(values=c(4),labels="No optimization",drop=FALSE,name="")+guides(colour=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))
#
#if (!file.exists("sim_stquant.png"))
#{                 
#  save_plot("sim_stquant.png",plot_sim_stquant,base_height=15,base_aspect_ratio=1.5)
#}
#
###Should I reorder here or not? I do not think so
##data$Sample=with(data,reorder(Sample,filtNABcovB_.U,FUN = function(x){return(median(x))}))
#plot_var_ciMean=ggplot(data=data,aes(y=filtNABcovB_.U,x=Sample))+geom_violin(aes(fill=dnagroup))+geom_point(data=data[data$condition %in% conditions_ciMean,],aes(color=factor(condition,levels=conditions_ciMean)),size=2,stat="summary",fun.y="mean")+geom_point(data=datanofilt,aes(x=sample,y=X.AN+X.BN+X.comun,shape=""),size=2.5,stroke=1.5)+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Sample")+scale_color_brewer(name="Filter",type="qual",palette=3,labels=names(conditions_ciMean))+labs(title="Number of variants. CI90 Mean selection")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_shape_manual(values=c(4),labels="No optimization",drop=FALSE,name="")+guides(colour=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))
#plot_var_ciMean_log=ggplot(data=data,aes(y=filtNABcovB_.U,x=Sample))+geom_violin(aes(fill=dnagroup))+geom_point(data=data[data$condition %in% conditions_ciMean,],aes(color=factor(condition,levels=conditions_ciMean)),size=2,stat="summary",fun.y="mean")+geom_point(data=datanofilt,aes(x=sample,y=X.AN+X.BN+X.comun,shape=""),size=2.5,stroke=1.5)+scale_y_log10(name="Number of variants")+scale_x_discrete(name="Sample")+scale_color_brewer(name="Filter",type="qual",palette=3,labels=names(conditions_ciMean))+labs(title="Number of variants. CI90 Mean selection")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_shape_manual(values=c(4),labels="No optimization",drop=FALSE,name="")+guides(colour=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))
#
#if(!file.exists("var_ciMean.png"))
#{
#  save_plot("var_ciMean.png",plot_var_ciMean,base_height=15,base_aspect_ratio=1.5)
#}
#
#if(!file.exists("var_ciMean2.png"))
#{
#  save_plot("var_ciMean_log.png",plot_var_ciMean_log,base_height=15,base_aspect_ratio=1.5)
#}
#
#plot_var_stquant=ggplot(data=data,aes(y=filtNABcovB_.U,x=Sample))+geom_violin(aes(fill=dnagroup))+geom_point(data=data[data$condition %in% conditions_stquant,],aes(color=factor(condition,levels=conditions_stquant)),size=2,stat="summary",fun.y="mean")+geom_point(data=datanofilt,aes(x=sample,y=X.AN+X.BN+X.comun,shape=""),size=2.5,stroke=1.5)+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Sample")+scale_color_brewer(name="Filter",type="qual",palette=3,labels=names(conditions_stquant))+labs(title="Number of variants. Q25 selection")+scale_fill_brewer(type = "seq",palette=8, name="DNA (ng)")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_shape_manual(values=c(4),labels="No optimization",drop=FALSE,name="")+guides(colour=guide_legend(order=1),shape=guide_legend(order=2),fill=guide_legend(order=3))
#
#if(!file.exists("var_stquant.png"))
#{
#  save_plot("var_stquant.png",plot_var_stquant,base_height=15,base_aspect_ratio=1.5)
#}
#
###Plots crossvalidation
#rescrossvalbest=resultscrossvalidation[resultscrossvalidation$condition==crossvalidation_lowciMean$condition[1],]
#rescrossvalbest=gather(rescrossvalbest,"lowciMean_test","lowciMean_training",key="set",value="sim")
#set=c("Test"="lowciMean_test","Training"="lowciMean_trainin")
#crossval_best=ggplot(data=rescrossvalbest,aes(x=as.numeric(k),y=sim,color=set))+geom_point()+scale_y_continuous(limits=c(0,1),name="Similarity (lowciMean)")+geom_smooth()+scale_x_continuous(name="Fold")+scale_color_brewer(type = "qual",palette = 1,labels=names(set),name="Set")
#
#if(!file.exists("crossval_best.png"))
#{
#  save_plot("crossval_best.png",crossval_best,base_height=15,base_aspect_ratio=1.5)
#}
#
###NAB problem
#
####Best filter for every condition given the cross-validation data
#bestfiltcombNABCrossTest=crossvalidation_lowciMean %>% mutate(conditionNAB=paste(NABmin_freq_alt,NABmin_reads_alternate,NABmax_coverage,sep="_")) %>% filter(cval=="test") %>% group_by(conditionNAB) %>% top_n(n=1,wt=mean) %>% summarise(condition=condition[1],mean=mean[1]) %>% arrange(desc(mean))
#if(!file.exists("bestfilters.csv"))
#{
#  write.table(bestfiltcombNABCrossTest,file="bestfilters.csv",quote = FALSE,row.names = FALSE,dec = ".",sep = ";")
#}
#
####Mean and sd for all the data at the same time
#sumdatabestfiltcombNABCrossTest=data %>% filter(condition %in% bestfiltcombNABCrossTest$condition)  %>% group_by(condition) %>% summarise(NABmin_freq_alt=dplyr::first(NABmin_freq_alt),NABmin_reads_alternate=dplyr::first(NABmin_reads_alternate), NABmax_coverage=dplyr::first(NABmax_coverage),meanfiltNABcovB_propU=mean(filtNABcovB_propU), meanfiltNABcovB_NU=mean(filtNABcovB_NU),sdfiltNABcovB_propU=sd(filtNABcovB_propU), sdfiltNABcovB_NU=sd(filtNABcovB_NU)) %>% arrange(desc(meanfiltNABcovB_propU))
#if(!file.exists("summbestfiltsNABallsamples.csv"))
#{
#  write.table(sumdatabestfiltcombNABCrossTest,file="summbestfiltsNABallsamples.csv",quote = FALSE,row.names = FALSE,dec = ".",sep = ";")
#}
#
#if(!file.exists("best_NABs_plot.png"))
#{
#  bestNABs=ggplot(data=sumdatabestfiltcombNABCrossTest,aes(y=meanfiltNABcovB_propU,x=meanfiltNABcovB_NU,shape=as.factor(NABmin_freq_alt)))+geom_point(stroke=1,size=2)+facet_grid(NABmin_reads_alternate~NABmax_coverage)+scale_x_continuous(name="Mean number of variants")+scale_y_continuous(name="Mean similarity")+ theme_bw() + labs(title="Best NAB conditions. Faceted by min coverage in N (columns) and max reads alternative (rows)") + scale_shape_manual(values = c(0,3,4),name="Max freq alternative")
#  save_plot("best_NABs_plot.png",bestNABs,base_height=12,base_aspect_ratio=1.78)
#}
#
####All data for the best filter for every NAB condition
#databestfiltcombNABCrossTest=data %>% filter(condition %in% bestfiltcombNABCrossTest$condition)
#if(!file.exists("databestfiltsNABCrossValidation.csv"))
#{
#  write.table(databestfiltcombNABCrossTest,file="databestfiltsNABCrossValidation.csv",quote = FALSE,row.names = FALSE,dec = ".",sep = ";")
#}

####OLD
################################
################################
# 
# 
# #Selection of filtering options attending to FiltN
# ##################################################
# 
# lowciMean_cond_FiltN=aggregate(filtN_prop_mean~conditionFilt,data,function(x) ciMean(x,conf = .9)[1])
# stquantile_cond_FiltN=aggregate(filtN_prop_mean~conditionFilt,data,function(x) quantile(x,p=.25))
# stquantile_cond_FiltN=stquantile_cond_FiltN[order(-stquantile_cond_FiltN$filtN_prop_mean),]
# lowciMean_cond_FiltN=lowciMean_cond_FiltN[order(-lowciMean_cond_FiltN$filtN_prop_mean),]
# conditions_FiltN_ciMean=head(lowciMean_cond_FiltN,10)$conditionFilt
# conditions_FiltN_stquant=head(stquantile_cond_FiltN,10)$conditionFilt
# conditions_FiltN_ciMean=factor(conditions_FiltN_ciMean,levels=conditions_FiltN_ciMean)
# conditions_FiltN_stquant=factor(conditions_FiltN_stquant,levels=conditions_FiltN_stquant)
# conditions_FiltN_ciMean %in% conditions_FiltN_stquant ##Are they the same?
# 
# #Selection of filtering options + NAB attending to FiltNAB
# ####################################################
# 
# lowciMean_cond_FiltNAB=aggregate(filtNAB_prop_mean~condition,data,function(x) ciMean(x,conf = .9)[1])
# stquantile_cond_FiltNAB=aggregate(filtNAB_prop_mean~condition,data,function(x) quantile(x,p=.25))
# stquantile_cond_FiltNAB=stquantile_cond_FiltNAB[order(-stquantile_cond_FiltNAB$filtNAB_prop_mean),]
# lowciMean_cond_FiltNAB=lowciMean_cond_FiltNAB[order(-lowciMean_cond_FiltNAB$filtNAB_prop_mean),]
# conditions_FiltNAB_ciMean=head(lowciMean_cond_FiltNAB,10)$condition
# conditions_FiltNAB_stquant=head(stquantile_cond_FiltNAB,10)$condition
# conditions_FiltNAB_ciMean=factor(conditions_FiltNAB_ciMean,levels=conditions_FiltNAB_ciMean)
# conditions_FiltNAB_stquant=factor(conditions_FiltNAB_stquant,levels=conditions_FiltNAB_stquant)
# conditions_FiltNAB_ciMean %in% conditions_FiltNAB_stquant ##Are they the same?
# 
# #Selection of filtering options + NAB attending to FiltNAB (no mean)
# ####################################################################
# 
# lowciMean_cond_FiltNAB=aggregate(filtNAB_prop~condition,data,function(x) ciMean(x,conf = .9)[1])
# stquantile_cond_FiltNAB=aggregate(filtNAB_prop~condition,data,function(x) quantile(x,p=.25))
# stquantile_cond_FiltNAB=stquantile_cond_FiltNAB[order(-stquantile_cond_FiltNAB$filtNAB_prop),]
# lowciMean_cond_FiltNAB=lowciMean_cond_FiltNAB[order(-lowciMean_cond_FiltNAB$filtNAB_prop),]
# conditions_FiltNAB_ciMean=head(lowciMean_cond_FiltNAB,10)$condition
# conditions_FiltNAB_stquant=head(stquantile_cond_FiltNAB,10)$condition
# conditions_FiltNAB_ciMean=factor(conditions_FiltNAB_ciMean,levels=conditions_FiltNAB_ciMean)
# conditions_FiltNAB_stquant=factor(conditions_FiltNAB_stquant,levels=conditions_FiltNAB_stquant)
# conditions_FiltNAB_ciMean %in% conditions_FiltNAB_stquant ##Are they the same?
# 
# #Percentile per Sample of the best filtering condition
# ######################################################
# Percentiles=c()
# condition=lowciMean_cond_FiltNAB[1,1]
# Samples=unique(data$Sample)
# for (i in 1:length(unique(data$Sample))){
# temp_data=(data[data$Sample==Samples[i],])[,c("filtNAB_prop","condition")]
# temp_data=temp_data[order(-temp_data$filtNAB_prop),]
# Percentiles[i]=100-(which(temp_data$condition==condition)/nrow(temp_data)*100)
# }
# results_percentiles=data.frame(Samples,Percentiles)
# 
# sink("percentiles.tex")
# cat("\\documentclass{article}
# \\usepackage[landscape]{geometry}
# \\usepackage{adjustbox}
# \\begin{document}
# \\begin{adjustbox}{width={\\textwidth},totalheight={\\textheight},keepaspectratio}%
# \\begin{tabular}{",paste0(rep("c",ncol(results_percentiles)),collapse=""),"}")
# print(xtable(results_percentiles),include.rownames=FALSE,only.contents=TRUE)
# cat("\\end{tabular}
# \\end{adjustbox}
# \\end{document}")
# sink()
# 
# 
# #EXPERIMENTAL: Selection of filtering options for NAB
# #####################################################
# 
# lowciMean_condNAB_FiltNAB=aggregate(filtNAB_prop_mean~conditionNAB,data,function(x) ciMean(x,conf = .9)[1])
# stquantile_condNAB_FiltNAB=aggregate(filtNAB_prop_mean~conditionNAB,data,function(x) quantile(x,p=.25))
# stquantile_condNAB_FiltNAB=stquantile_condNAB_FiltNAB[order(-stquantile_condNAB_FiltNAB$filtNAB_prop_mean),]
# lowciMean_condNAB_FiltNAB=lowciMean_condNAB_FiltNAB[order(-lowciMean_condNAB_FiltNAB$filtNAB_prop_mean),]
# conditionsNAB_FiltNAB_ciMean=head(lowciMean_condNAB_FiltNAB,10)$conditionNAB
# conditionsNAB_FiltNAB_stquant=head(stquantile_condNAB_FiltNAB,10)$conditionNAB
# conditionsNAB_FiltNAB_ciMean=factor(conditionsNAB_FiltNAB_ciMean,levels=conditionsNAB_FiltNAB_ciMean)
# conditionsNAB_FiltNAB_stquant=factor(conditionsNAB_FiltNAB_stquant,levels=conditionsNAB_FiltNAB_stquant)
# conditionsNAB_FiltNAB_ciMean %in% conditionsNAB_FiltNAB_stquant ##Are they the same?
# 
# 
# #Combined: Plots without scaling
# #################################
# 
# #FiltN
# plot_sim_FiltN_ciMean=ggplot(data=data,aes(y=filtN_prop_mean,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$conditionFilt %in% conditions_FiltN_ciMean,],aes(color=factor(conditionFilt,levels=conditions_FiltN_ciMean)),size=2,stat="summary",fun.y="mean")+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Somatic similarity. CI90 Mean selection")
# save_plot("sim_FiltN_ciMean.png",plot_sim_FiltN_ciMean,base_height=10,base_aspect_ratio=1.5)
# plot_sim_FiltN_stquant=ggplot(data=data,aes(y=filtN_prop_mean,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$conditionFilt %in% conditions_FiltN_stquant,],aes(color=factor(conditionFilt,levels=conditions_FiltN_stquant)),size=2,stat="summary",fun.y="mean")+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Somatic similarity. Q25 selection")
# save_plot("sim_FiltN_stquant.png",plot_sim_FiltN_stquant,base_height=10,base_aspect_ratio=1.5)
# 
# plot_var_FiltN_ciMean=ggplot(data=data,aes(y=filtN_N,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$conditionFilt %in% conditions_FiltN_ciMean,],aes(color=factor(conditionFilt,levels=conditions_FiltN_ciMean)),size=2,stat="summary",fun.y="mean")+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Number of variants. CI90 Mean selection")
# save_plot("var_FiltN_ciMean.png",plot_var_FiltN_ciMean,base_height=10,base_aspect_ratio=1.5)
# plot_var_FiltN_stquant=ggplot(data=data,aes(y=filtN_N,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$conditionFilt %in% conditions_FiltN_stquant,],aes(color=factor(conditionFilt,levels=conditions_FiltN_stquant)),size=2,stat="summary",fun.y="mean")+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Number of variants. Q25 selection")
# save_plot("var_FiltN_stquant.png",plot_var_FiltN_stquant,base_height=10,base_aspect_ratio=1.5)
# 
# #FiltNAB
# 
# plot_sim_FiltNAB_ciMean=ggplot(data=data,aes(y=filtNAB_prop_mean,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$condition %in% conditions_FiltNAB_ciMean,],aes(color=factor(condition,levels=conditions_FiltNAB_ciMean)),size=2)+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Somatic similarity -NAB . CI90 Mean selection")
# save_plot("sim_FiltNAB_ciMean.png",plot_sim_FiltNAB_ciMean,base_height=10,base_aspect_ratio=1.5)
# plot_sim_FiltNAB_stquant=ggplot(data=data,aes(y=filtNAB_prop_mean,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$condition %in% conditions_FiltNAB_stquant,],aes(color=factor(condition,levels=conditions_FiltNAB_stquant)),size=2)+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Somatic similarity -NAB. Q25 selection")
# save_plot("sim_FiltNAB_stquant.png",plot_sim_FiltNAB_stquant,base_height=10,base_aspect_ratio=1.5)
# 
# plot_var_FiltNAB_ciMean=ggplot(data=data,aes(y=filtNAB_N,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$condition %in% conditions_FiltNAB_ciMean,],aes(color=factor(condition,levels=conditions_FiltNAB_ciMean)),size=2)+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Number of variants -NAB. CI90 Mean selection")
# save_plot("var_FiltNAB_ciMean.png",plot_var_FiltNAB_ciMean,base_height=10,base_aspect_ratio=1.5)
# plot_var_FiltNAB_stquant=ggplot(data=data,aes(y=filtNAB_N,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$condition %in% conditions_FiltNAB_stquant,],aes(color=factor(condition,levels=conditions_FiltNAB_stquant)),size=2)+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Number of variants -NAB. Q25 selection")
# save_plot("var_FiltNAB_stquant.png",plot_var_FiltNAB_stquant,base_height=10,base_aspect_ratio=1.5)
# 
# #NAB
# 
# plot_sim_NAB_ciMean=ggplot(data=data,aes(y=filtNAB_prop_mean,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$conditionNAB %in% conditionsNAB_FiltNAB_ciMean,],aes(color=factor(conditionNAB,conditionsNAB_FiltNAB_ciMean)),size=2,stat="summary",fun.y="mean")+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Mean somatic similarity, NAB filter. CI90 Mean selection")
# save_plot("sim_NAB_ciMean.png",plot_sim_NAB_ciMean,base_height=10,base_aspect_ratio=1.5)
# plot_sim_NAB_stquant=ggplot(data=data,aes(y=filtNAB_prop_mean,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$conditionNAB %in% conditionsNAB_FiltNAB_stquant,],aes(color=factor(conditionNAB,conditionsNAB_FiltNAB_stquant)),size=2,stat="summary",fun.y="mean")+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Mean somatic similarity, NAB filter. Q25 selection")
# save_plot("sim_NAB_stquant.png",plot_sim_NAB_stquant,base_height=10,base_aspect_ratio=1.5)
# 
# plot_var_NAB_ciMean=ggplot(data=data,aes(y=filtNAB_N,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$conditionNAB %in% conditionsNAB_FiltNAB_ciMean,],aes(color=factor(conditionNAB,conditionsNAB_FiltNAB_ciMean)),size=2,stat="summary",fun.y="mean")+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Mean number of variants, NAB filter. CI90 Mean selection")
# save_plot("var_NAB_ciMean.png",plot_var_NAB_ciMean,base_height=10,base_aspect_ratio=1.5)
# plot_var_NAB_stquant=ggplot(data=data,aes(y=filtNAB_N,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$conditionNAB %in% conditionsNAB_FiltNAB_stquant,],aes(color=factor(conditionNAB,conditionsNAB_FiltNAB_stquant)),size=2,stat="summary",fun.y="mean")+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Mean number of variants, NAB filter. Q25 selection")
# save_plot("var_NAB_stquant.png",plot_var_NAB_stquant,base_height=10,base_aspect_ratio=1.5)
# 
# ##################################################################################
# #WORKING NOW
# ##################################################################################
# 
# #Plots TsTv
# ############
# 
# lm_plot=ggplot(data=data[round(runif(1000000,1,nrow(data))),],aes(x=filtNAB_prop_mean,y=TsTv_filtNAB,color=as.factor(Sample)))+geom_point()+geom_smooth(aes(fill=as.factor(Sample)),color="black",alpha=0.8,method="lm")+geom_smooth(color="red",method="lm")+scale_x_continuous(name="Similarity")+scale_y_continuous(name="TsTv")+scale_fill_discrete(name="Sample")+scale_color_discrete(guide=FALSE)
# loess_plot=ggplot(data=data[round(runif(1000000,1,nrow(data))),],aes(x=filtNAB_prop_mean,y=TsTv_filtNAB,color=as.factor(Sample)))+geom_point()+geom_smooth(aes(fill=as.factor(Sample)),color="black",alpha=0.8)+geom_smooth(color="red")+scale_x_continuous(name="Similarity")+scale_y_continuous(name="TsTv")+scale_fill_discrete(name="Sample")+scale_color_discrete(guide=FALSE)
# save_plot("lm.png",lm_plot,base_height=10,base_aspect_ratio=1.5)
# save_plot("loess.png",loess_plot,base_height=10,base_aspect_ratio=1.5)
# 
# #Selection of filtering options + NAB attending to FiltNAB, using both similarity and tstv
# ###########################################################################################
# 
# filtNAB_prop_mean_maxima=aggregate(filtNAB_prop_mean~Sample,data,function(x) max(x))
# filtNAB_prop_mean_minima=aggregate(filtNAB_prop_mean~Sample,data,function(x) min(x))
# data$filtNAB_prop_mean_scaled=apply(data,1,function(x){(as.numeric(x[33])-filtNAB_prop_mean_minima[filtNAB_prop_mean_minima$Sample==x[1],2])/(filtNAB_prop_mean_maxima[filtNAB_prop_mean_maxima$Sample==x[1],2]-filtNAB_prop_mean_minima[filtNAB_prop_mean_minima$Sample==x[1],2])})
# 
# TsTv_filtNAB_maxima=aggregate(TsTv_filtNAB~Sample,data,function(x) max(x))
# TsTv_filtNAB_minima=aggregate(TsTv_filtNAB~Sample,data,function(x) min(x))
# data$TsTv_filtNAB_scaled=apply(data,1,function(x){(as.numeric(x[39])-TsTv_filtNAB_minima[TsTv_filtNAB_minima$Sample==x[1],2])/(TsTv_filtNAB_maxima[TsTv_filtNAB_maxima$Sample==x[1],2]-TsTv_filtNAB_minima[TsTv_filtNAB_minima$Sample==x[1],2])})
# 
# select_variable=(data$TsTv_filtNAB_scaled+data$filtNAB_prop_mean_scaled)/2.0
# data$select_variable=select_variable
# 
# lowciMean_cond_FiltNABTsTv=aggregate(select_variable~condition,data,function(x) ciMean(x,conf = .9)[1])
# stquantile_cond_FiltNABTsTv=aggregate(select_variable~condition,data,function(x) quantile(x,p=.25))
# stquantile_cond_FiltNABTsTv=stquantile_cond_FiltNABTsTv[order(-stquantile_cond_FiltNABTsTv$select_variable),]
# lowciMean_cond_FiltNABTsTv=lowciMean_cond_FiltNABTsTv[order(-lowciMean_cond_FiltNABTsTv$select_variable),]
# conditions_FiltNABTsTv_ciMean=head(lowciMean_cond_FiltNABTsTv,10)$condition
# conditions_FiltNABTsTv_stquant=head(stquantile_cond_FiltNABTsTv,10)$condition
# conditions_FiltNABTsTv_ciMean=factor(conditions_FiltNABTsTv_ciMean,levels=conditions_FiltNABTsTv_ciMean)
# conditions_FiltNABTsTv_stquant=factor(conditions_FiltNABTsTv_stquant,levels=conditions_FiltNABTsTv_stquant)
# conditions_FiltNABTsTv_ciMean %in% conditions_FiltNABTsTv_stquant ##Are they the same?
# 
# plot_sim_FiltNABTsTv_ciMean=ggplot(data=data,aes(y=filtNAB_prop_mean,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$condition %in% conditions_FiltNABTsTv_ciMean,],aes(color=factor(condition,levels=conditions_FiltNABTsTv_ciMean)),size=2)+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Somatic similarity -NAB . CI90 Mean selection")
# save_plot("sim_FiltNABTsTv_ciMean.png",plot_sim_FiltNABTsTv_ciMean,base_height=10,base_aspect_ratio=1.5)
# plot_sim_FiltNABTsTv_stquant=ggplot(data=data,aes(y=filtNAB_prop_mean,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$condition %in% conditions_FiltNABTsTv_stquant,],aes(color=factor(condition,levels=conditions_FiltNABTsTv_stquant)),size=2)+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Somatic similarity -NAB, Q25 selection")
# save_plot("sim_FiltNABTsTv_stquant.png",plot_sim_FiltNABTsTv_stquant,base_height=10,base_aspect_ratio=1.5)
# 
# plot_var_FiltNABTsTv_ciMean=ggplot(data=data,aes(y=filtNAB_N,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$condition %in% conditions_FiltNABTsTv_ciMean,],aes(color=factor(condition,levels=conditions_FiltNABTsTv_ciMean)),size=2)+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Number of variants -NAB. CI90 Mean selection")
# save_plot("var_FiltNABTsTv_ciMean.png",plot_var_FiltNABTsTv_ciMean,base_height=10,base_aspect_ratio=1.5)
# plot_var_FiltNABTsTv_stquant=ggplot(data=data,aes(y=filtNAB_N,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$condition %in% conditions_FiltNABTsTv_stquant,],aes(color=factor(condition,levels=conditions_FiltNABTsTv_stquant)),size=2)+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Number of variants -NAB. Q25 selection")
# save_plot("var_FiltNABTsTv_stquant.png",plot_var_FiltNABTsTv_stquant,base_height=10,base_aspect_ratio=1.5)
# 
# 
# ##Selection of filtering options attending to TsTvNAB per replicate
# ###################################################################
# 
# filts=aggregate(TsTv_filtNAB~Sample,data,function(x) sort(x,decreasing=TRUE)[1])
# 
# output_dataframe=data[0,]
# for (i in 1:length(filts$Sample)){
# 	temp_dataframe=(data[data$Sample==filts[i,1],])[which(data[data$Sample==filts[i,1],]$TsTv_filtNAB==filts[filts[i,1],2]),]
# 	temp_dataframe=temp_dataframe[order(temp_dataframe$filtNAB_N,decreasing=TRUE),][1,]
# 	output_dataframe=rbind(output_dataframe,temp_dataframe)
# }
# ciMean(output_dataframe$filtNAB_prop_mean,conf=.9)[1]
# 
# ##Independent selection of filtering options
# ################################################
# 
# conds=names(data)[2:(n_sim_parameters+1)]
# cond_cimean=""
# cond_stquant=""
# for (i in 1:length(conds)){
# temp_list1=aggregate(as.formula(paste0("filtNAB_prop_mean~",conds[i])),data,function(x) ciMean(x,conf = .9)[1])
# temp_list2=aggregate(as.formula(paste0("filtNAB_prop_mean~",conds[i])),data,function(x) ciMean(x,conf = .9)[1])
# temp_list1[which.max(temp_list1[,2]),1]
# temp_list2[which.max(temp_list1[,2]),1]
# cond_cimean=paste0(cond_cimean,temp_list1[which.max(temp_list1[,2]),1],",")
# cond_stquant=paste0(cond_stquant,temp_list2[which.max(temp_list2[,2]),1],",")
# }
# cond_cimean=substr(cond_cimean, 1, nchar(cond_cimean)-1)
# cond_stquant=substr(cond_stquant, 1, nchar(cond_stquant)-1)
# ciMean(data[gsub(" ","",data$condition)==cond_cimean,"filtNAB_prop_mean"],conf=.9)[1]
# quantile(data[gsub(" ","",data$condition)==cond_stquant,"filtNAB_prop_mean"],p=.25)
# 
# 
# ##Comparative via bootstrapping ##ATTENTION: very computationally intesive
# ########################################################################
# 
# total_samples=unique(data$Sample)
# red_data=data[,c("Sample","filtNAB_prop_mean","condition","select_variable")]
# 
# #mybs=function(data,indices,realdata,tstvdata)
# mybs=function(data,indices,realdata)
# {
# 	my_samples=data[indices]
# 	#tstv_result=ciMean((tstvdata[indices,])$filtNAB_prop_mean,conf = .9)[1]
# 	#return(tstv_result)
# 	output_dataframe=realdata[0,]
# 	for (i in 1:length(my_samples)){
# 		temp_dataframe=realdata[realdata$Sample==my_samples[i],]
# 		temp_dataframe$Sample=i
# 		output_dataframe=rbind(output_dataframe,temp_dataframe)
# 	}
# 	lowciMean_cond_FiltNAB=aggregate(filtNAB_prop_mean~condition,output_dataframe,function(x) ciMean(x,conf = .9)[1])
# 	stquantile_cond_FiltNAB=aggregate(filtNAB_prop_mean~condition,output_dataframe,function(x) quantile(x,p=.25))
# 	cond_mean=lowciMean_cond_FiltNAB[which.max(lowciMean_cond_FiltNAB$filtNAB_prop_mean),1]
# 	cond_stquant=stquantile_cond_FiltNAB[which.max(stquantile_cond_FiltNAB$filtNAB_prop_mean),1]
# 	rm(lowciMean_cond_FiltNAB)
# 	rm(stquantile_cond_FiltNAB)
# 	cimean=ciMean(realdata[realdata$condition==cond_mean,"filtNAB_prop_mean"],conf=.9)[1]
# 	stquant=quantile(realdata[realdata$condition==cond_stquant,"filtNAB_prop_mean"],p=.25)
# 
# 	lowciMean_cond_FiltNABTsTv=aggregate(select_variable~condition,output_dataframe,function(x) ciMean(x,conf = .9)[1])
# 	stquantile_cond_FiltNABTsTv=aggregate(select_variable~condition,output_dataframe,function(x) quantile(x,p=.25))
# 	cond_mean=lowciMean_cond_FiltNABTsTv[which.max(lowciMean_cond_FiltNABTsTv$select_variable),1]
# 	cond_stquant=stquantile_cond_FiltNABTsTv[which.max(stquantile_cond_FiltNABTsTv$select_variable),1]
# 	rm(lowciMean_cond_FiltNABTsTv)
# 	rm(stquantile_cond_FiltNABTsTv)
# 	citstvmean=ciMean(realdata[realdata$condition==cond_mean,"filtNAB_prop_mean"],conf=.9)[1]
# 	sttstvquant=quantile(realdata[realdata$condition==cond_stquant,"filtNAB_prop_mean"],p=.25)
# 
# 
# 	conds=names(realdata)[2:(n_sim_parameters+1)]
# 	cond_cimean=""
# 	cond_stquant=""
# 	for (i in 1:length(conds)){
# 	temp_list1=aggregate(as.formula(paste0("filtNAB_prop_mean~",conds[i])),output_dataframe,function(x) ciMean(x,conf = .9)[1])
# 	temp_list2=aggregate(as.formula(paste0("filtNAB_prop_mean~",conds[i])),output_dataframe,function(x) ciMean(x,conf = .9)[1])
# 	temp_list1[which.max(temp_list1[,2]),1]
# 	temp_list2[which.max(temp_list2[,2]),1]
# 	cond_cimean=paste0(cond_cimean,temp_list1[which.max(temp_list1[,2]),1],",")
# 	cond_stquant=paste0(cond_stquant,temp_list2[which.max(temp_list2[,2]),1],",")
# 	}
# 	cond_cimean=substr(cond_cimean, 1, nchar(cond_cimean)-1)
# 	cond_stquant=substr(cond_stquant, 1, nchar(cond_stquant)-1)
# 	ciindmean=ciMean(realdata[gsub(" ","",realdata$condition)==cond_cimean,"filtNAB_prop_mean"],conf=.9)[1]
# 	stindquant=quantile(realdata[gsub(" ","",realdata$condition)==cond_stquant,"filtNAB_prop_mean"],p=.25)
# 
# 	return(c(cimean,stquant,citstvmean,sttstvquant,ciindmean,stindquant))
# }
# 
# results_bootstrapping=boot(data=total_samples,statistic=mybs,R=104,realdata=data,parallel="multicore",ncpus=8)
# 
# 
# #Summarized results
# #output_dataframe_angelo=data[0,]
# #for (i in 1:length(filts$Sample)){
# #	temp_dataframe=data[data$Sample==filts[i,1],]
# #	output_dataframe_angelo=rbind(output_dataframe_angelo,(temp_dataframe[order(temp_dataframe$filtNAB_prop,decreasing=TRUE),])[1:10000,])
# #}
# 
# 
# ######################################################
# #Best parameters per Sample
# ######################################################
# Samples=unique(data$Sample)
# sorted_filtN=data[order(data$filtN_prop_mean,decreasing=TRUE),]
# sorted_filtNAB=data[order(data$filtNAB_prop_mean,decreasing=TRUE),]
# 
# ####By ctbrown, obtained in stackoverflow
# cbind.all <- function (...)
# {
#     nm <- list(...)
#     nm <- lapply(nm, as.matrix)
#     n <- max(sapply(nm, nrow))
#     do.call(cbind, lapply(nm, function(x) rbind(x, matrix(, n -
#         nrow(x), ncol(x)))))
# }
# 
# table_filtN=data.frame()
# table_filtNAB=data.frame()
# 
# for (Sample in Samples)
# {
# 	table_filtN=cbind.all(table_filtN,head(unique(sorted_filtN[sorted_filtN$Sample==Sample,]$conditionFilt),n=10))
# 	table_filtNAB=cbind.all(table_filtNAB,head(unique(sorted_filtNAB[sorted_filtNAB$Sample==Sample,]$condition),n=10))
# }
# colnames(table_filtN)=Samples
# colnames(table_filtNAB)=Samples
# 
# sink("10best_filtN.tex")
# cat("\\documentclass{article}
# \\usepackage[landscape]{geometry}
# \\usepackage{adjustbox}
# \\begin{document}
# \\begin{adjustbox}{width={\\textwidth},totalheight={\\textheight},keepaspectratio}%
# \\begin{tabular}{",paste0(rep("c",ncol(table_filtN)),collapse=""),"}")
# print(xtable(table_filtN),include.rownames=FALSE,only.contents=TRUE)
# cat("\\end{tabular}
# \\end{adjustbox}
# \\end{document}")
# sink()
# 
# sink("10best_filtNAB.tex")
# cat("\\documentclass{article}
# \\usepackage[landscape]{geometry}
# \\usepackage{adjustbox}
# \\begin{document}
# \\begin{adjustbox}{width={\\textwidth},totalheight={\\textheight},keepaspectratio}%
# \\begin{tabular}{",paste0(rep("c",ncol(table_filtNAB)),collapse=""),"}")
# print(xtable(table_filtNAB),include.rownames=FALSE,only.contents=TRUE)
# cat("\\end{tabular}
# \\end{adjustbox}
# \\end{document}")
# sink()
# 
# #####################################################
# #Logit regression
# #####################################################
# 
# data=data[order(data$condition),]
# 
# xrepNAB=data[,1:n_sim_parameters] ##The parameters are repeated, for the six Samples.
# xrepN=data[,1:(n_sim_parameters-n_sim_parameters_NAB)]
# 
# yrepfiltNAB=data$filtNAB_prop_mean
# yrepfiltN=data$filtN_prop_mean
# 
# xrepfiltNAB=apply(xrepfiltNAB,2,function(x){as.numeric(gsub("-1","0",gsub("0_","",x)))})
# xrepfiltNAB=apply(xrepfiltNAB,2,function(x){(x-min(x))/(max(x)-min(x))})
# xrepfiltNAB=data.frame(xrepfiltNAB)
# xrepfiltNAB=cbind(xrepfiltNAB,data$Sample)
# 
# xrepfiltN=apply(xrepfiltN,2,function(x){as.numeric(gsub("-1","0",gsub("0_","",x)))})
# xrepfiltN=apply(xrepfiltN,2,function(x){(x-min(x))/(max(x)-min(x))})
# xrepfiltN=data.frame(xrepfiltN)
# xrepfiltN=cbind(xrepfiltN,data$Sample)
# 
# ##Stepwise glm
# ###############
# initial_aic_glm_filtN_model=glm(yrepfiltN~1,family=binomial(link='logit'))
# aic_glm_filtN_model=step(model,as.formula(paste0("~ .+ (",paste0("xrepfiltN[,\"",colnames(xrep),"\"]",collapse="+"),")^2")))
# initial_aic_glm_filtNAB_model=glm(yrepfiltNAB~1,family=binomial(link='logit'))
# aic_glm_filtNAB_model=step(model,as.formula(paste0("~ .+ (",paste0("xrepfiltNAB[,\"",colnames(xrep),"\"]",collapse="+"),")^2")))
# 
# ##glm full model
# ############################
# glm_filtN_all=glm(yrepfiltN~xrepfiltN,family=binomial(link='logit'))
# glm_filtNAB_all=glm(yrepfiltNAB~xrepfiltNAB,family=binomial(link='logit'))
# #glm_filtN_all_interact=glm(as.formula(paste0("yrepfiltN ~  (",paste0("xrepfiltN[,\"",colnames(xrepfiltN),"\"]",collapse="+"),")^2")),family=binomial(link='logit'))
# #glm_filtNAB_all_interact=glm(as.formula(paste0("yrepfiltNAB ~  (",paste0("xrepfiltNAB[,\"",colnames(xrepfiltNAB),"\"]",collapse="+"),")^2")),family=binomial(link='logit'))
# 
# 
# #########################################################
# #Regression by Sample
# #########################################################
# 
# 
# for (Sample in Samples) {
# 	xfiltNAB=data[data$Sample==Sample,1:n_sim_parameters] ##The parameters are repeated, for the six Samples, just keeping one
# 	xfiltN=data[data$Sample==Sample,1:(n_sim_parameters-n_sim_parameters_NAB)]
# 	xfiltNAB=apply(xfiltNAB,2,function(x){as.numeric(gsub("-1","0",gsub("0_","",x)))})
# 	xfiltNAB=apply(xfiltNAB,2,function(x){(x-min(x))/(max(x)-min(x))})
# 	xfiltN=apply(xfiltN,2,function(x){as.numeric(gsub("-1","0",gsub("0_","",x)))})
# 	xfiltN=apply(xfiltN,2,function(x){(x-min(x))/(max(x)-min(x))})
# 	#assign(paste0("glm_filtN_interact_",Sample),glm(as.formula(paste0("data[data$Sample==Sample,]$filtN_prop_mean ~  (",paste0("xfiltN[,\"",colnames(xfiltN),"\"]",collapse="+"),")^2")),family=binomial(link='logit')))
# 	#assign(paste0("glm_filtNAB_interact_",Sample),glm(as.formula(paste0("data[data$Sample==Sample,]$filtNAB_prop_mean ~  (",paste0("xfiltNAB[,\"",colnames(xfiltNAB),"\"]",collapse="+"),")^2")),family=binomial(link='logit')))
# 	assign(paste0("glm_filtN_",Sample),glm(data[data$Sample==Sample,]$filtN_prop_mean ~ xfiltN,family=binomial(link='logit')))
# 	assign(paste0("glm_filtNAB_",Sample),glm(data[data$Sample==Sample,]$filtNAB_prop_mean ~ xfiltNAB,family=binomial(link='logit')))
# }
# 
# ###Plots
# ########
# 
# parameters_filtNAB=names(data[,1:n_sim_parameters])
# parameters_filtN=names(data[,1:(n_sim_parameters-n_sim_parameters_NAB)])
# 
# for (Sample in Samples) {
# 	for (parameter in parameters_filtN) {
# 			temp=as.numeric((data[data$Sample==Sample,])[parameter][[1]])
# 			x=seq(min(temp),max(temp),length.out=1000)
# 			myglm=eval(parse(text=paste0("glm_filtN_",Sample)))
# 			y=1/(1+exp(-(myglm$coefficients[1]+myglm$coefficients[paste0("xfiltN",parameter)]*x)))
# 			logit_data=data.frame(x,y)
# 			plot=ggplot(data=data[data$Sample==Sample,],aes_string(x=parameter,y="filtN_prop_mean"))+geom_violin(aes_string(group=parameter))+geom_boxplot(aes_string(group=parameter),width=.05)+geom_line(data=logit_data,aes(x=x,y=y),color="red",size=1)+labs(title=paste0(Sample,"_",parameter,"_filtN"))+scale_y_continuous(name="Similarity")
# 			save_plot(paste0(Sample,"_",parameter,"_filtN",".png"),plot,base_height=10,base_aspect_ratio=1.5)
# 		}
# 	for (parameter in parameters_filtNAB) {
# 			temp=as.numeric((data[data$Sample==Sample,])[parameter][[1]])
# 			x=seq(min(temp),max(temp),length.out=1000)
# 			myglm=eval(parse(text=paste0("glm_filtNAB_",Sample)))
# 			y=1/(1+exp(-(myglm$coefficients[1]+myglm$coefficients[paste0("xfiltNAB",parameter)]*x)))
# 			logit_data=data.frame(x,y)
# 			plot=ggplot(data=data[data$Sample==Sample,],aes_string(x=parameter,y="filtNAB_prop_mean"))+geom_violin(aes_string(group=parameter))+geom_boxplot(aes_string(group=parameter),width=.05)+geom_line(data=logit_data,aes(x=x,y=y),color="red",size=1)+labs(title=paste0(Sample,"_",parameter,"_filtNAB"))+scale_y_continuous(name="Similarity")
# 			save_plot(paste0(Sample,"_",parameter,"_filtNAB",".png"),plot,base_height=10,base_aspect_ratio=1.5)
# 		}
# }
# 
# 
# #########################################################
# # Grouping Samples
# #########################################################
# res_dunn.test=dunn.test(data$filtNAB_prop_mean, g=as.factor(data$Sample)) ##Kluskal-wallis + post-hocs. It takes forever!!!!
# 
# 
# 
# 
# ####OLD, may be outdated
# ##########################
# 
# #Selection of filtering options scaling the distributions per Sample
# ####################################################################
# 
# maxima=data.frame(t(tapply(data$filtN_prop_mean,data$Sample,max)))
# minima=data.frame(t(tapply(data$filtN_prop_mean,data$Sample,min)))
# data$filtN_prop_mean_norm=apply(data,1,function(x) {as.numeric((as.numeric(x[8])-minima[as.character(x[9])])/(maxima[as.character(x[9])]-minima[as.character(x[9])]))})
# lowciMean_cond_norm=aggregate(filtN_prop_mean_norm~condition,data,function(x) ciMean(x,conf = .9)[1])
# stquantile_cond_norm=aggregate(filtN_prop_mean_norm~condition,data,function(x) quantile(x,p=.25))
# stquantile_cond_norm=stquantile_cond_norm[order(-stquantile_cond_norm$filtN_prop_mean_norm),]
# lowciMean_cond_norm=lowciMean_cond_norm[order(-lowciMean_cond_norm$filtN_prop_mean_norm),]
# conditions_norm=head(lowciMean_cond_norm,10)$condition
# conditions_norm2=head(stquantile_cond_norm,10)$condition
# ucondnorm=union(conditions_norm,conditions_norm2)
# icondnorm=intersect(conditions_norm,conditions_norm2)
# 
# #Combined: Plots with scaling
# ###############################
# 
# plot_norm=ggplot(data=data,aes(y=filtN_prop_mean,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$condition %in% icondnorm,],aes(color=condition),size=2)+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")
# plot_norm2=ggplot(data=data,aes(y=filtN_N,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$condition %in% icondnorm,],aes(color=condition),size=2)+scale_y_continuous(name="Variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")
# save_plot("similarity_norm.png",plot_norm,base_height=10,base_aspect_ratio=1.5)
# save_plot("variants_norm.png",plot_norm2,base_height=10,base_aspect_ratio=1.5)
# 
# 
# #A and B separated
# #############################################################################
# dataA=data[,c(seq(1,n_sim_parameters+4),n_sim_parameters+7,n_sim_parameters+8)]
# dataA$Sample=paste(dataA$Sample,"A",sep="_")
# dataB=data[,c(seq(1,n_sim_parameters+2),seq(n_sim_parameters+5,n_sim_parameters+8))]
# dataB$Sample=paste(dataB$Sample,"B",sep="_")
# 
# mynames=names(dataA)
# mynames[n_sim_parameters+1]="filtN_N_comb"
# mynames[n_sim_parameters+2]="filtN_prop_mean_comb"
# mynames[n_sim_parameters+4]="filtN_N"
# mynames[n_sim_parameters+3]="filtN_prop_mean"
# names(dataA)=mynames
# names(dataB)=mynames
# datasep=rbind(dataA,dataB)
# 
# #Selection of filtering options separated
# ##########################################
# 
# lowciMean_cond_sep=aggregate(filtN_prop_mean~condition,datasep,function(x) ciMean(x,conf = .9)[1])
# stquantile_cond_sep=aggregate(filtN_prop_mean~condition,datasep,function(x) quantile(x,p=.25))
# stquantile_cond_sep=stquantile_cond_sep[order(-stquantile_cond_sep$filtN_prop_mean),]
# lowciMean_cond_sep=lowciMean_cond_sep[order(-lowciMean_cond_sep$filtN_prop_mean),]
# conditions_sep=head(lowciMean_cond_sep,10)$condition
# conditions2_sep=head(stquantile_cond_sep,10)$condition
# sort(conditions_sep)==sort(conditions2_sep) ##Different!
# 
# #By replicate: Plots without scaling
# ####################################
# 
# plot_byrep_combconds=ggplot(data=datasep,aes(y=filtN_prop_mean,x=as.factor(Sample)))+geom_violin()+geom_point(data=datasep[datasep$condition %in% conditions,],aes(color=condition),size=2)+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Combined filter. Similarity")
# save_plot("similarity_byrep_combonds.png",plot_byrep_combconds,base_height=10,base_aspect_ratio=1.5)
# plot_byrep_ciconds=ggplot(data=datasep,aes(y=filtN_prop_mean,x=as.factor(Sample)))+geom_violin()+geom_point(data=datasep[datasep$condition %in% conditions_sep,],aes(color=condition),size=2)+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="ci 90% filter. Similarity")
# save_plot("similarity_byrep_ciconds.png",plot_byrep_ciconds,base_height=10,base_aspect_ratio=1.5)
# plot_byrep_qconds=ggplot(data=datasep,aes(y=filtN_prop_mean,x=as.factor(Sample)))+geom_violin()+geom_point(data=datasep[datasep$condition %in% conditions2_sep,],aes(color=condition),size=2)+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="1st quantile filter. Similarity")
# save_plot("similarity_byrep_qconds.png",plot_byrep_qconds,base_height=10,base_aspect_ratio=1.5)
# plot_byrep_combconds_var=ggplot(data=datasep,aes(y=filtN_N,x=as.factor(Sample)))+geom_violin()+geom_point(data=datasep[datasep$condition %in% conditions,],aes(color=condition),size=2)+scale_y_continuous(name="N common variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Combined filter. Variants")
# save_plot("variants_byrep_combonds.png",plot_byrep_combconds_var,base_height=10,base_aspect_ratio=1.5)
# plot_byrep_ciconds_var=ggplot(data=datasep,aes(y=filtN_N,x=as.factor(Sample)))+geom_violin()+geom_point(data=datasep[datasep$condition %in% conditions_sep,],aes(color=condition),size=2)+scale_y_continuous(name="N common variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="ci 90% filter. Variants")
# save_plot("variants_byrep_ciconds.png",plot_byrep_ciconds_var,base_height=10,base_aspect_ratio=1.5)
# plot_byrep_qconds_var=ggplot(data=datasep,aes(y=filtN_N,x=as.factor(Sample)))+geom_violin()+geom_point(data=datasep[datasep$condition %in% conditions2_sep,],aes(color=condition),size=2)+scale_y_continuous(name="N common variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="1st quantile filter. Variants")
# save_plot("variants_byrep_qconds.png",plot_byrep_qconds_var,base_height=10,base_aspect_ratio=1.5)
# 
# 
# 
# #################################################################################################################################################################################
# #################################################################################################################################################################################
# #PROBLEM DATA
# #################################################################################################################################################################################
# #################################################################################################################################################################################
# data_problem=read.csv("results_problem.csv")
# 
# ##DANGER: HARDCODED DATA!!!!
# data_problem$type=c("Pure DCIS","Pure DCIS","Pure DCIS","Pure DCIS","Pure DCIS","Pure DCIS","Pure DCIS","Pure DCIS","Pure DCIS","Pure DCIS","Pure DCIS","DCIS with adjacent invasive","DCIS with adjacent invasive","DCIS with adjacent invasive","DCIS with adjacent invasive","DCIS with adjacent invasive","DCIS with adjacent invasive","DCIS with adjacent invasive","DCIS with adjacent invasive","DCIS with adjacent invasive","DCIS with adjacent invasive","DCIS with adjacent invasive")
# data_problem$type=as.factor(data_problem$type)
# types=levels(data_problem$type)
# types=rev(types)
# #data_problem$filtNAB_NU_somatic=c(29,2,13,4,7,7,11,44,6,3,1,98,13,73,16,0,6,4,5,1,51,89)
# data_problem=data_problem[-2,]
# data_problem=data_problem[-15,]
# data_problem$filtNAB_.U_somatic=c(35,22,370,56,12,24,61,366,18,20,122,170,97,45,55,26,41,16,71,133)
# 
# ###
# 
# #Violin plot mean mutational burden
# ###################################
# 
# violin_mut_burden=ggplot(data=data_problem, aes(y=filtNAB_._mean,x=type,fill=type))+geom_violin()+geom_boxplot(width=.2,outlier.size=0.5)+scale_y_continuous(name=NULL)+labs(title="Number of mutations")+scale_x_discrete(name=NULL,limits=types,labels=c("Pure DCIS", "DCIS with\n adjacent invasive"))+scale_fill_manual(values=c("#009E73","#0072B2"))+guides(fill=FALSE)+theme(text=element_text(size=6,face="bold"),axis.text.y=element_text(size=4),axis.text.x=element_text(size=4,face="bold"),plot.margin=unit(c(0.2,0.1,0.1,0.1),"cm"),plot.title=element_text(size=5))
# save_plot("mut_burden.png",violin_mut_burden,base_height=2,base_aspect_ratio=.5)
# 
# #Violin plot protein coding mutational burden (total)
# #############################################
# 
# p_value=t.test(filtNAB_.U_somatic~type,data=data_problem)$p.value
# p_valuew=wilcox.test(filtNAB_.U_somatic~type,data=data_problem)$p.value
# 
# #violin_coding_burden=ggplot(data=data_problem, aes(y=filtNAB_.U_somatic,x=type,fill=type))+geom_violin()+geom_boxplot(width=.08,outlier.size=0.5)+scale_y_continuous(name=NULL)+labs(title="Number of somatic\n non-synonymous\n mutations")+scale_x_discrete(name=NULL,limits=types,labels=c("Pure DCIS", "DCIS with\n adjacent invasive"))+scale_fill_manual(values=c("#009E73","#0072B2"))+guides(fill=FALSE)+theme(text=element_text(size=6,face="bold"),axis.text.y=element_text(size=4),axis.text.x=element_text(size=4,face="bold"),plot.margin=unit(c(0.2,0.1,0.1,0.1),"cm"),plot.title=element_text(size=5))+annotate("text",x=1.2,y=100,label=paste0("T-test p-value = ",round(p_value,digits=3)),size=1.1)+annotate("text",x=1.2,y=95,label=paste0("Wilcoxon p-value =",round(p_valuew,digits=3)),size=1.1)
# violin_coding_burden=ggplot(data=data_problem, aes(y=filtNAB_.U_somatic,x=type,fill=type))+geom_violin()+geom_boxplot(width=.08,outlier.size=0.5)+scale_y_continuous(name=NULL)+labs(title="Number of somatic\n non-synonymous\n mutations")+scale_x_discrete(name=NULL,limits=types,labels=c("Pure DCIS", "DCIS with\n adjacent invasive"))+scale_fill_manual(values=c("#009E73","#0072B2"))+guides(fill=FALSE)+theme(text=element_text(size=6,face="bold"),axis.text.y=element_text(size=4),axis.text.x=element_text(size=4,face="bold"),plot.margin=unit(c(0.2,0.1,0.1,0.1),"cm"),plot.title=element_text(size=5))
# 
# save_plot("coding_burden.pdf",violin_coding_burden,base_height=2,base_aspect_ratio=.5)
# 
# 
# #Similarity plot
# ###################################
# 
# violin_similarity=ggplot(data=data_problem, aes(y=1-filtNAB_propU,x=type,fill=type))+geom_violin()+geom_boxplot(width=.2,outlier.size=0.5)+scale_y_continuous(name=NULL)+labs(title="Divergence\n between regions")+scale_x_discrete(name=NULL,limits=types,labels=c("Pure DCIS", "DCIS with\n adjacent invasive"))+scale_fill_manual(values=c("#009E73","#0072B2"))+guides(fill=FALSE)+theme(text=element_text(size=6,face="bold"),axis.text.y=element_text(size=4),axis.text.x=element_text(size=4,face="bold"),plot.margin=unit(c(0.2,0.1,0.1,0.1),"cm"),plot.title=element_text(size=5))
# save_plot("similarity.png",violin_similarity,base_height=2,base_aspect_ratio=.5)
# 
# 
# #Folds mutational burden
# ####################################
# 
# folds=function(x, y) {
# 	results=c()
# 	if(length(x) != length(y))
# 	{
# 		stop("The two vectors must have the same length")
# 	}
# 	for (i in 1:length(x))
# 	{
# 		if(x[i]>=y[i])
# 		{
# 			results[i]=x[i]/y[i]
# 		}
# 		else
# 		{
# 			results[i]=y[i]/x[i]
# 		}
# 	}
# 	return(results)
# }
# 
# data_problem$foldsfiltNAB_.=folds(data_problem$AfiltNAB_.,data_problem$BfiltNAB_.)
# p_value=t.test(foldsfiltNAB_.~type,data=data_problem)$p.value
# p_valuew=wilcox.test(foldsfiltNAB_.~type,data=data_problem)$p.value
# 
# violin_folds_mutburden=ggplot(data=data_problem, aes(y=foldsfiltNAB_.,x=type,fill=type))+geom_violin()+geom_boxplot(width=.05,outlier.size=0.5)+scale_y_continuous(name=NULL)+labs(title="Fold-differences in \n number of mutations \n between regions")+scale_x_discrete(name=NULL,limits=types,labels=c("Pure DCIS", "DCIS with\n adjacent invasive"))+scale_fill_manual(values=c("#009E73","#0072B2"))+guides(fill=FALSE)+theme(text=element_text(size=6,face="bold"),axis.text.y=element_text(size=4),axis.text.x=element_text(size=4,face="bold"),plot.margin=unit(c(0.2,0.1,0.1,0.1),"cm"),plot.title=element_text(size=5))+annotate("text",x=1.2,y=6,label=paste0("T-test p-value = ",round(p_value,digits=3)),size=1.1)+annotate("text",x=1.2,y=5.7,label=paste0("Wilcoxon p-value = ",round(p_valuew,digits=3)),size=1.1)
# save_plot("folds.png",violin_folds_mutburden,base_height=2,base_aspect_ratio=.5)
# 
# 
# 
# ########################################################
# ########################################################
# #ROIS
# ########################################################
# ########################################################
# setwd("~/Desktop/dcis")
# data=read.csv("rois.csv")
# out_data=as.data.frame(list(id = vector("character", nrow(data)), test = vector("character",nrow(data)), p_value=vector("numeric",nrow(data)), t_test.p_value = vector("numeric", nrow(data)),wilcox.p_value=vector("numeric",nrow(data)), f_test.p_value = vector("numeric", nrow(data)),shapiro_test_A.p_value=vector("numeric",nrow(data)),shapiro_test_B.p_value=vector("numeric",nrow(data))),stringsAsFactors=FALSE)
# 
# for (i in 3:length(names(data)))
# {
# 	name=as.character(names(data)[i])
# 	shapiro_test_A.p_value=as.numeric(tryCatch({shapiro.test(data[data$Type=="Pure",names(data)[i]])$p.value},error=function(err){return(0)}))
# 	shapiro_test_B.p_value=as.numeric(tryCatch({shapiro.test(data[data$Type=="Inv",names(data)[i]])$p.value},error=function(err){return(0)}))
# 	f_test.p_value=as.numeric(var.test(as.formula(paste0((names(data))[i],"~Type")),data=data)$p.value)
# 	wilcox.p_value=as.numeric(wilcox.test(as.formula(paste0((names(data))[i],"~Type")),data=data)$p.value)
# 	if(is.nan(f_test.p_value) || f_test.p_value>0.05)
# 	{
# 	t_test.p_value=as.numeric(t.test(as.formula(paste0((names(data))[i],"~Type")),data=data,var.equal=TRUE)$p.value)
# 	}
# 	else
# 	{
# 	t_test.p_value=as.numeric(t.test(as.formula(paste0((names(data))[i],"~Type")),data=data)$p.value)
# 	}
# 	if(shapiro_test_A.p_value>0.05 && shapiro_test_B.p_value>0.05)
# 	{
# 		p_value=t_test.p_value
# 		test="t_test"
# 	}
# 	else
# 	{
# 		p_value=wilcox.p_value
# 		test="wilcoxon"
# 	}
# 	  row=c(name,test,p_value,t_test.p_value,wilcox.p_value,f_test.p_value,shapiro_test_A.p_value,shapiro_test_B.p_value)
# 	length(shapiro_test_A.p_value)
#   out_data[i-2,]=row
# }
# write.csv(out_data,file="rois_out.csv",row.names=FALSE)
# 
# ##
# #Violin plots
# ##
# 
# data$Type=as.factor(data$Type)
# types=levels(data$Type)
# types=rev(types)
# 
# p_value=as.numeric(out_data[out_data$id=="stroma.ALDH.0",]$p_value)
# plot_aldh=ggplot(data=data, aes(y=stroma.ALDH.0,x=Type,fill=Type))+geom_violin()+geom_boxplot(width=.05,outlier.size=0.5)+scale_y_continuous(name=NULL)+labs(title="Difference in % of\n slide without ALDH stain\n in the stroma")+scale_x_discrete(name=NULL,limits=types,labels=c("Pure DCIS", "DCIS with\n adjacent invasive"))+scale_fill_manual(values=c("#009E73","#0072B2"))+guides(fill=FALSE)+theme(text=element_text(size=6,face="bold"),axis.text.y=element_text(size=4),axis.text.x=element_text(size=4,face="bold"),plot.margin=unit(c(0.2,0.1,0.1,0.1),"cm"),plot.title=element_text(size=5))+annotate("text",x=1.8,y=85,label=paste0("Wilcoxon \np-value = ",round(p_value,digits=3)),size=1.1)
# save_plot("aldh0.png",plot_aldh,base_height=2,base_aspect_ratio=.5)
# 
# p_value=as.numeric(out_data[out_data$id=="stroma.ALDH.2",]$p_value)
# plot_aldh2=ggplot(data=data, aes(y=stroma.ALDH.2,x=Type,fill=Type))+geom_violin()+geom_boxplot(width=.05,outlier.size=0.5)+scale_y_continuous(name=NULL)+labs(title="Difference in % of\n slide with intensity 2\nof ALDH stain\n in the stroma")+scale_x_discrete(name=NULL,limits=types,labels=c("Pure DCIS", "DCIS with\n adjacent invasive"))+scale_fill_manual(values=c("#009E73","#0072B2"))+guides(fill=FALSE)+theme(text=element_text(size=6,face="bold"),axis.text.y=element_text(size=4),axis.text.x=element_text(size=4,face="bold"),plot.margin=unit(c(0.2,0.1,0.1,0.1),"cm"),plot.title=element_text(size=5))+annotate("text",x=1.8,y=85,label=paste0("Wilcoxon \np-value = ",round(p_value,digits=3)),size=1.1)
# save_plot("aldh2.png",plot_aldh2,base_height=2,base_aspect_ratio=.5)
# 
# p_value=as.numeric(out_data[out_data$id=="Mean.Number.Cancer.HS",]$p_value)
# plot_mean_cancer_hs=ggplot(data=data, aes(y=Mean.Number.Cancer.HS,x=Type,fill=Type))+geom_violin()+geom_boxplot(width=.05,outlier.size=0.5)+scale_y_continuous(name=NULL)+labs(title="Mean of the number\nof cancer hotspots")+scale_x_discrete(name=NULL,limits=types,labels=c("Pure DCIS", "DCIS with\n adjacent invasive"))+scale_fill_manual(values=c("#009E73","#0072B2"))+guides(fill=FALSE)+theme(text=element_text(size=6,face="bold"),axis.text.y=element_text(size=4),axis.text.x=element_text(size=4,face="bold"),plot.margin=unit(c(0.2,0.1,0.1,0.1),"cm"),plot.title=element_text(size=5))+annotate("text",x=1.8,y=85,label=paste0("Wilcoxon \np-value = ",round(p_value,digits=3)),size=1.1)
# save_plot("meanhs.png",plot_mean_cancer_hs,base_height=2,base_aspect_ratio=.5)
# 
# p_valuef=as.numeric(out_data[out_data$id=="StD.number.Cancer.HS",]$f_test.p_value)
# plot_sd_cancer_hs=ggplot(data=data, aes(y=StD.number.Cancer.HS,x=Type,fill=Type))+geom_violin()+geom_boxplot(width=.05,outlier.size=0.5)+scale_y_continuous(name=NULL)+labs(title="SD of the numbern\nof cancer hotspots")+scale_x_discrete(name=NULL,limits=types,labels=c("Pure DCIS", "DCIS with\n adjacent invasive"))+scale_fill_manual(values=c("#009E73","#0072B2"))+guides(fill=FALSE)+theme(text=element_text(size=6,face="bold"),axis.text.y=element_text(size=4),axis.text.x=element_text(size=4,face="bold"),plot.margin=unit(c(0.2,0.1,0.1,0.1),"cm"),plot.title=element_text(size=5))+annotate("text",x=1.8,y=85,label=paste0("F-test \np-value = ",round(p_valuef,digits=3)),size=1.1)
# save_plot("sdhs.png",plot_sd_cancer_hs,base_height=2,base_aspect_ratio=.5)
# 
# 
# #p_value=as.numeric(out_data[out_data$id=="Mean.Number.Cancer.HS",]$p_value)
# #plot_mean_cancer_hs=ggplot(data=data, aes(y=Mean.Number.Cancer.HS,x=Type,fill=Type))+geom_violin()+geom_boxplot(width=.05,outlier.size=0.5)+scale_y_continuous(name=NULL)+labs(title="Mean of the number\nof cancer hotspots")+scale_x_discrete(name=NULL,limits=types,labels=c("Pure DCIS", "DCIS with\n adjacent invasive"))+scale_fill_manual(values=c("#fbb4b9","#7a0177"))+guides(fill=FALSE)+theme(text=element_text(size=6,face="bold"),axis.text.y=element_text(size=4),axis.text.x=element_text(size=4,face="bold"),plot.margin=unit(c(0.2,0.1,0.1,0.1),"cm"),plot.title=element_text(size=5))+annotate("text",x=1.8,y=85,label=paste0("Wilcoxon \np-value = ",round(p_value,digits=3)),size=1.1)
# #save_plot("meanhs.pdf",plot_mean_cancer_hs,base_height=2,base_aspect_ratio=.5)
# #
# #p_valuef=as.numeric(out_data[out_data$id=="StD.number.Cancer.HS",]$f_test.p_value)
# #plot_sd_cancer_hs=ggplot(data=data, aes(y=StD.number.Cancer.HS,x=Type,fill=Type))+geom_violin()+geom_boxplot(width=.05,outlier.size=0.5)+scale_y_continuous(name=NULL)+labs(title="SD of the numbern\nof cancer hotspots")+scale_x_discrete(name=NULL,limits=types,labels=c("Pure DCIS", "DCIS with\n adjacent invasive"))+scale_fill_manual(values=c("#fbb4b9","#7a0177"))+guides(fill=FALSE)+theme(text=element_text(size=6,face="bold"),axis.text.y=element_text(size=4),axis.text.x=element_text(size=4,face="bold"),plot.margin=unit(c(0.2,0.1,0.1,0.1),"cm"),plot.title=element_text(size=5))+annotate("text",x=1.8,y=85,label=paste0("F-test \np-value = ",round(p_valuef,digits=4)),size=1.1)
# #save_plot("sdhs.pdf",plot_sd_cancer_hs,base_height=2,base_aspect_ratio=.5)
# 
# p_value=as.numeric(out_data[out_data$id=="DCIS.P.FAK.3",]$p_value)
# plot_pfak=ggplot(data=data, aes(y=DCIS.P.FAK.3,x=Type,fill=Type))+geom_violin()+geom_boxplot(width=.05,outlier.size=0.5)+scale_y_continuous(name=NULL)+labs(title="Difference in % of\nslide with the maximum\n intensity of pFAK stain")+scale_x_discrete(name=NULL,limits=types,labels=c("Pure DCIS", "DCIS with\n adjacent invasive"))+scale_fill_manual(values=c("#009E73","#0072B2"))+guides(fill=FALSE)+theme(text=element_text(size=6,face="bold"),axis.text.y=element_text(size=4),axis.text.x=element_text(size=4,face="bold"),plot.margin=unit(c(0.2,0.1,0.1,0.1),"cm"),plot.title=element_text(size=5))+annotate("text",x=1.8,y=85,label=paste0("Wilcoxon \np-value = ",round(p_value,digits=3)),size=1.1)
# save_plot("pfak.png",plot_pfak,base_height=2,base_aspect_ratio=.5)
# 
# p_value=as.numeric(out_data[out_data$id=="DCIS.RANK.1",]$p_value)
# plot_rank=ggplot(data=data, aes(y=DCIS.RANK.1,x=Type,fill=Type))+geom_violin()+geom_boxplot(width=.05,outlier.size=0.5)+scale_y_continuous(name=NULL)+labs(title="Difference in % of\nslide with the DCIS rank 1")+scale_x_discrete(name=NULL,limits=types,labels=c("Pure DCIS", "DCIS with\n adjacent invasive"))+scale_fill_manual(values=c("#009E73","#0072B2"))+guides(fill=FALSE)+theme(text=element_text(size=6,face="bold"),axis.text.y=element_text(size=4),axis.text.x=element_text(size=4,face="bold"),plot.margin=unit(c(0.2,0.1,0.1,0.1),"cm"),plot.title=element_text(size=5))+annotate("text",x=1.8,y=95,label=paste0("Wilcoxon \np-value = ",round(p_value,digits=3)),size=1.1)
# save_plot("rank.png",plot_rank,base_height=2,base_aspect_ratio=.5)
# 
# p_value=as.numeric(out_data[out_data$id=="DCIS.RANK.0",]$p_value)
# plot_rank0=ggplot(data=data, aes(y=DCIS.RANK.0,x=Type,fill=Type))+geom_violin()+geom_boxplot(width=.05,outlier.size=0.5)+scale_y_continuous(name=NULL)+labs(title="Difference in % of\nslide with the DCIS rank 0")+scale_x_discrete(name=NULL,limits=types,labels=c("Pure DCIS", "DCIS with\n adjacent invasive"))+scale_fill_manual(values=c("#009E73","#0072B2"))+guides(fill=FALSE)+theme(text=element_text(size=6,face="bold"),axis.text.y=element_text(size=4),axis.text.x=element_text(size=4,face="bold"),plot.margin=unit(c(0.2,0.1,0.1,0.1),"cm"),plot.title=element_text(size=5))+annotate("text",x=1.8,y=95,label=paste0("Wilcoxon \np-value = ",round(p_value,digits=3)),size=1.1)
# save_plot("rank0.png",plot_rank0,base_height=2,base_aspect_ratio=.5)
# 
# p_value=as.numeric(out_data[out_data$id=="DCIS.PR.1",]$p_value)
# plot_pr=ggplot(data=data, aes(y=DCIS.PR.1,x=Type,fill=Type))+geom_violin()+geom_boxplot(width=.05,outlier.size=0.5)+scale_y_continuous(name=NULL)+labs(title="Difference in % of\nslide with intensity 1\nof PR stain")+scale_x_discrete(name=NULL,limits=types,labels=c("Pure DCIS", "DCIS with\n adjacent invasive"))+scale_fill_manual(values=c("#009E73","#0072B2"))+guides(fill=FALSE)+theme(text=element_text(size=6,face="bold"),axis.text.y=element_text(size=4),axis.text.x=element_text(size=4,face="bold"),plot.margin=unit(c(0.2,0.1,0.1,0.1),"cm"),plot.title=element_text(size=5))+annotate("text",x=1.9,y=7,label=paste0("Wilcoxon \np-value = ",round(p_value,digits=3)),size=1.1)
# save_plot("pr.png",plot_pr,base_height=2,base_aspect_ratio=.5)
# 
# p_value=as.numeric(out_data[out_data$id=="DCIS.PR.2",]$p_value)
# plot_pr2=ggplot(data=data, aes(y=DCIS.PR.2,x=Type,fill=Type))+geom_violin()+geom_boxplot(width=.05,outlier.size=0.5)+scale_y_continuous(name=NULL)+labs(title="Difference in % of\nslide with intensity 2\nof PR stain")+scale_x_discrete(name=NULL,limits=types,labels=c("Pure DCIS", "DCIS with\n adjacent invasive"))+scale_fill_manual(values=c("#009E73","#0072B2"))+guides(fill=FALSE)+theme(text=element_text(size=6,face="bold"),axis.text.y=element_text(size=4),axis.text.x=element_text(size=4,face="bold"),plot.margin=unit(c(0.2,0.1,0.1,0.1),"cm"),plot.title=element_text(size=5))+annotate("text",x=1.9,y=7,label=paste0("Wilcoxon \np-value = ",round(p_value,digits=3)),size=1.1)
# save_plot("pr2.png",plot_pr2,base_height=2,base_aspect_ratio=.5)
# 
# p_value=as.numeric(out_data[out_data$id=="DCIS.CA9.1",]$p_value)
# plot_ca9=ggplot(data=data, aes(y=DCIS.CA9.1,x=Type,fill=Type))+geom_violin()+geom_boxplot(width=.05,outlier.size=0.5)+scale_y_continuous(name=NULL)+labs(title="Difference in % of\nslide with intensity 1\nof CA9 stain")+scale_x_discrete(name=NULL,limits=types,labels=c("Pure DCIS", "DCIS with\n adjacent invasive"))+scale_fill_manual(values=c("#009E73","#0072B2"))+guides(fill=FALSE)+theme(text=element_text(size=6,face="bold"),axis.text.y=element_text(size=4),axis.text.x=element_text(size=4,face="bold"),plot.margin=unit(c(0.2,0.1,0.1,0.1),"cm"),plot.title=element_text(size=5))+annotate("text",x=1.9,y=40,label=paste0("Wilcoxon \np-value = ",round(p_value,digits=3)),size=1.1)
# save_plot("ca9.png",plot_ca9,base_height=2,base_aspect_ratio=.5)
