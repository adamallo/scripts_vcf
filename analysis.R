library(lsr)
library(cowplot)
library(xtable)
library(dunn.test)
library(boot)

# Data parsing #HARDCODED
############################################

###Conf variables
#dir="~/res3_tab"
dir="~/Desktop/dcis"
n_sim_parameters=9 ###WARNING: CHANGE THIS ACCORDINGLY WITH YOUR DATA. TOTAL NUMBER OF PARAMETERS
n_sim_parameters_NAB=3 ###WARNING: CHANGE THIS ACCORDINGLY WITH YOUR DATA. PARAMETERS RELATIVE TO NAB
setwd(dir)
#################

data=read.csv("results.csv")
#B3=read.csv("B3_withtstv_tab.csv")
#B6=read.csv("B6_withtstv_tab.csv")
#D5=read.csv("D5_withtstv_tab.csv")
#D8=read.csv("D8_withtstv_tab.csv")
#DCIS64=read.csv("DCIS64_withtstv_tab.csv")
#K12=read.csv("K12_withtstv_tab.csv")
#B3$Sample=rep("B3",nrow(B3))
#B6$Sample=rep("B6",nrow(B6))
#D5$Sample=rep("D5",nrow(D5))
#D8$Sample=rep("D8",nrow(D8))
#DCIS64$Sample=rep("DCIS64",nrow(DCIS64))
#K12$Sample=rep("K12",nrow(K12))

#data=rbind(B3,B6,D5,D8,DCIS64,K12)
#data$condition=apply(data,1,function(x) paste(x[1:n_sim_parameters],collapse=","))
#data$conditionFilt=apply(data,1,function(x) paste(x[1:(n_sim_parameters-n_sim_parameters_NAB)],collapse=","))
#data$conditionNAB=apply(data,1,function(x) paste(x[(n_sim_parameters-n_sim_parameters_NAB+1):n_sim_parameters],collapse=","))
data$condition=apply(data,1,function(x) paste(x[2:(n_sim_parameters+1)],collapse=","))
data$conditionFilt=apply(data,1,function(x) paste(x[2:(n_sim_parameters-n_sim_parameters_NAB+1)],collapse=","))
data$conditionNAB=apply(data,1,function(x) paste(x[(n_sim_parameters-n_sim_parameters_NAB+2):(n_sim_parameters+1)],collapse=","))

#Selection of filtering options attending to FiltN
##################################################

lowciMean_cond_FiltN=aggregate(filtN_prop_mean~conditionFilt,data,function(x) ciMean(x,conf = .9)[1])
stquantile_cond_FiltN=aggregate(filtN_prop_mean~conditionFilt,data,function(x) quantile(x,p=.25))
stquantile_cond_FiltN=stquantile_cond_FiltN[order(-stquantile_cond_FiltN$filtN_prop_mean),]
lowciMean_cond_FiltN=lowciMean_cond_FiltN[order(-lowciMean_cond_FiltN$filtN_prop_mean),]
conditions_FiltN_ciMean=head(lowciMean_cond_FiltN,10)$conditionFilt
conditions_FiltN_stquant=head(stquantile_cond_FiltN,10)$conditionFilt
conditions_FiltN_ciMean=factor(conditions_FiltN_ciMean,levels=conditions_FiltN_ciMean)
conditions_FiltN_stquant=factor(conditions_FiltN_stquant,levels=conditions_FiltN_stquant)
conditions_FiltN_ciMean %in% conditions_FiltN_stquant ##Are they the same?

#Selection of filtering options + NAB attending to FiltNAB
####################################################

lowciMean_cond_FiltNAB=aggregate(filtNAB_prop_mean~condition,data,function(x) ciMean(x,conf = .9)[1])
stquantile_cond_FiltNAB=aggregate(filtNAB_prop_mean~condition,data,function(x) quantile(x,p=.25))
stquantile_cond_FiltNAB=stquantile_cond_FiltNAB[order(-stquantile_cond_FiltNAB$filtNAB_prop_mean),]
lowciMean_cond_FiltNAB=lowciMean_cond_FiltNAB[order(-lowciMean_cond_FiltNAB$filtNAB_prop_mean),]
conditions_FiltNAB_ciMean=head(lowciMean_cond_FiltNAB,10)$condition
conditions_FiltNAB_stquant=head(stquantile_cond_FiltNAB,10)$condition
conditions_FiltNAB_ciMean=factor(conditions_FiltNAB_ciMean,levels=conditions_FiltNAB_ciMean)
conditions_FiltNAB_stquant=factor(conditions_FiltNAB_stquant,levels=conditions_FiltNAB_stquant)
conditions_FiltNAB_ciMean %in% conditions_FiltNAB_stquant ##Are they the same?

#Selection of filtering options + NAB attending to FiltNAB (no mean)
####################################################################

lowciMean_cond_FiltNAB=aggregate(filtNAB_prop~condition,data,function(x) ciMean(x,conf = .9)[1])
stquantile_cond_FiltNAB=aggregate(filtNAB_prop~condition,data,function(x) quantile(x,p=.25))
stquantile_cond_FiltNAB=stquantile_cond_FiltNAB[order(-stquantile_cond_FiltNAB$filtNAB_prop),]
lowciMean_cond_FiltNAB=lowciMean_cond_FiltNAB[order(-lowciMean_cond_FiltNAB$filtNAB_prop),]
conditions_FiltNAB_ciMean=head(lowciMean_cond_FiltNAB,10)$condition
conditions_FiltNAB_stquant=head(stquantile_cond_FiltNAB,10)$condition
conditions_FiltNAB_ciMean=factor(conditions_FiltNAB_ciMean,levels=conditions_FiltNAB_ciMean)
conditions_FiltNAB_stquant=factor(conditions_FiltNAB_stquant,levels=conditions_FiltNAB_stquant)
conditions_FiltNAB_ciMean %in% conditions_FiltNAB_stquant ##Are they the same?

#Percentile per Sample of the best filtering condition
######################################################
Percentiles=c()
condition=lowciMean_cond_FiltNAB[1,1]
Samples=unique(data$Sample)
for (i in 1:length(unique(data$Sample))){
temp_data=(data[data$Sample==Samples[i],])[,c("filtNAB_prop","condition")]
temp_data=temp_data[order(-temp_data$filtNAB_prop),]
Percentiles[i]=100-(which(temp_data$condition==condition)/nrow(temp_data)*100)
}
results_percentiles=data.frame(Samples,Percentiles)

sink("percentiles.tex")
cat("\\documentclass{article}
\\usepackage[landscape]{geometry}
\\usepackage{adjustbox}
\\begin{document}
\\begin{adjustbox}{width={\\textwidth},totalheight={\\textheight},keepaspectratio}%
\\begin{tabular}{",paste0(rep("c",ncol(results_percentiles)),collapse=""),"}")
print(xtable(results_percentiles),include.rownames=FALSE,only.contents=TRUE)
cat("\\end{tabular}
\\end{adjustbox}
\\end{document}")
sink()


#EXPERIMENTAL: Selection of filtering options for NAB
#####################################################

lowciMean_condNAB_FiltNAB=aggregate(filtNAB_prop_mean~conditionNAB,data,function(x) ciMean(x,conf = .9)[1])
stquantile_condNAB_FiltNAB=aggregate(filtNAB_prop_mean~conditionNAB,data,function(x) quantile(x,p=.25))
stquantile_condNAB_FiltNAB=stquantile_condNAB_FiltNAB[order(-stquantile_condNAB_FiltNAB$filtNAB_prop_mean),]
lowciMean_condNAB_FiltNAB=lowciMean_condNAB_FiltNAB[order(-lowciMean_condNAB_FiltNAB$filtNAB_prop_mean),]
conditionsNAB_FiltNAB_ciMean=head(lowciMean_condNAB_FiltNAB,10)$conditionNAB
conditionsNAB_FiltNAB_stquant=head(stquantile_condNAB_FiltNAB,10)$conditionNAB
conditionsNAB_FiltNAB_ciMean=factor(conditionsNAB_FiltNAB_ciMean,levels=conditionsNAB_FiltNAB_ciMean)
conditionsNAB_FiltNAB_stquant=factor(conditionsNAB_FiltNAB_stquant,levels=conditionsNAB_FiltNAB_stquant)
conditionsNAB_FiltNAB_ciMean %in% conditionsNAB_FiltNAB_stquant ##Are they the same?


#Combined: Plots without scaling
#################################

#FiltN
plot_sim_FiltN_ciMean=ggplot(data=data,aes(y=filtN_prop_mean,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$conditionFilt %in% conditions_FiltN_ciMean,],aes(color=factor(conditionFilt,levels=conditions_FiltN_ciMean)),size=2,stat="summary",fun.y="mean")+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Somatic similarity. CI90 Mean selection")
save_plot("sim_FiltN_ciMean.png",plot_sim_FiltN_ciMean,base_height=10,base_aspect_ratio=1.5)
plot_sim_FiltN_stquant=ggplot(data=data,aes(y=filtN_prop_mean,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$conditionFilt %in% conditions_FiltN_stquant,],aes(color=factor(conditionFilt,levels=conditions_FiltN_stquant)),size=2,stat="summary",fun.y="mean")+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Somatic similarity. Q25 selection")
save_plot("sim_FiltN_stquant.png",plot_sim_FiltN_stquant,base_height=10,base_aspect_ratio=1.5)

plot_var_FiltN_ciMean=ggplot(data=data,aes(y=filtN_N,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$conditionFilt %in% conditions_FiltN_ciMean,],aes(color=factor(conditionFilt,levels=conditions_FiltN_ciMean)),size=2,stat="summary",fun.y="mean")+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Number of variants. CI90 Mean selection")
save_plot("var_FiltN_ciMean.png",plot_var_FiltN_ciMean,base_height=10,base_aspect_ratio=1.5)
plot_var_FiltN_stquant=ggplot(data=data,aes(y=filtN_N,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$conditionFilt %in% conditions_FiltN_stquant,],aes(color=factor(conditionFilt,levels=conditions_FiltN_stquant)),size=2,stat="summary",fun.y="mean")+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Number of variants. Q25 selection")
save_plot("var_FiltN_stquant.png",plot_var_FiltN_stquant,base_height=10,base_aspect_ratio=1.5)

#FiltNAB

plot_sim_FiltNAB_ciMean=ggplot(data=data,aes(y=filtNAB_prop_mean,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$condition %in% conditions_FiltNAB_ciMean,],aes(color=factor(condition,levels=conditions_FiltNAB_ciMean)),size=2)+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Somatic similarity -NAB . CI90 Mean selection")
save_plot("sim_FiltNAB_ciMean.png",plot_sim_FiltNAB_ciMean,base_height=10,base_aspect_ratio=1.5)
plot_sim_FiltNAB_stquant=ggplot(data=data,aes(y=filtNAB_prop_mean,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$condition %in% conditions_FiltNAB_stquant,],aes(color=factor(condition,levels=conditions_FiltNAB_stquant)),size=2)+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Somatic similarity -NAB. Q25 selection")
save_plot("sim_FiltNAB_stquant.png",plot_sim_FiltNAB_stquant,base_height=10,base_aspect_ratio=1.5)

plot_var_FiltNAB_ciMean=ggplot(data=data,aes(y=filtNAB_N,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$condition %in% conditions_FiltNAB_ciMean,],aes(color=factor(condition,levels=conditions_FiltNAB_ciMean)),size=2)+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Number of variants -NAB. CI90 Mean selection")
save_plot("var_FiltNAB_ciMean.png",plot_var_FiltNAB_ciMean,base_height=10,base_aspect_ratio=1.5)
plot_var_FiltNAB_stquant=ggplot(data=data,aes(y=filtNAB_N,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$condition %in% conditions_FiltNAB_stquant,],aes(color=factor(condition,levels=conditions_FiltNAB_stquant)),size=2)+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Number of variants -NAB. Q25 selection")
save_plot("var_FiltNAB_stquant.png",plot_var_FiltNAB_stquant,base_height=10,base_aspect_ratio=1.5)

#NAB

plot_sim_NAB_ciMean=ggplot(data=data,aes(y=filtNAB_prop_mean,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$conditionNAB %in% conditionsNAB_FiltNAB_ciMean,],aes(color=factor(conditionNAB,conditionsNAB_FiltNAB_ciMean)),size=2,stat="summary",fun.y="mean")+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Mean somatic similarity, NAB filter. CI90 Mean selection")
save_plot("sim_NAB_ciMean.png",plot_sim_NAB_ciMean,base_height=10,base_aspect_ratio=1.5)
plot_sim_NAB_stquant=ggplot(data=data,aes(y=filtNAB_prop_mean,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$conditionNAB %in% conditionsNAB_FiltNAB_stquant,],aes(color=factor(conditionNAB,conditionsNAB_FiltNAB_stquant)),size=2,stat="summary",fun.y="mean")+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Mean somatic similarity, NAB filter. Q25 selection")
save_plot("sim_NAB_stquant.png",plot_sim_NAB_stquant,base_height=10,base_aspect_ratio=1.5)

plot_var_NAB_ciMean=ggplot(data=data,aes(y=filtNAB_N,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$conditionNAB %in% conditionsNAB_FiltNAB_ciMean,],aes(color=factor(conditionNAB,conditionsNAB_FiltNAB_ciMean)),size=2,stat="summary",fun.y="mean")+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Mean number of variants, NAB filter. CI90 Mean selection")
save_plot("var_NAB_ciMean.png",plot_var_NAB_ciMean,base_height=10,base_aspect_ratio=1.5)
plot_var_NAB_stquant=ggplot(data=data,aes(y=filtNAB_N,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$conditionNAB %in% conditionsNAB_FiltNAB_stquant,],aes(color=factor(conditionNAB,conditionsNAB_FiltNAB_stquant)),size=2,stat="summary",fun.y="mean")+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Mean number of variants, NAB filter. Q25 selection")
save_plot("var_NAB_stquant.png",plot_var_NAB_stquant,base_height=10,base_aspect_ratio=1.5)

##################################################################################
#WORKING NOW
##################################################################################

#Plots TsTv
############

lm_plot=ggplot(data=data[round(runif(1000000,1,nrow(data))),],aes(x=filtNAB_prop_mean,y=TsTv_filtNAB,color=as.factor(Sample)))+geom_point()+geom_smooth(aes(fill=as.factor(Sample)),color="black",alpha=0.8,method="lm")+geom_smooth(color="red",method="lm")+scale_x_continuous(name="Similarity")+scale_y_continuous(name="TsTv")+scale_fill_discrete(name="Sample")+scale_color_discrete(guide=FALSE)
loess_plot=ggplot(data=data[round(runif(1000000,1,nrow(data))),],aes(x=filtNAB_prop_mean,y=TsTv_filtNAB,color=as.factor(Sample)))+geom_point()+geom_smooth(aes(fill=as.factor(Sample)),color="black",alpha=0.8)+geom_smooth(color="red")+scale_x_continuous(name="Similarity")+scale_y_continuous(name="TsTv")+scale_fill_discrete(name="Sample")+scale_color_discrete(guide=FALSE)
save_plot("lm.png",lm_plot,base_height=10,base_aspect_ratio=1.5)
save_plot("loess.png",loess_plot,base_height=10,base_aspect_ratio=1.5)

#Selection of filtering options + NAB attending to FiltNAB, using both similarity and tstv
###########################################################################################

filtNAB_prop_mean_maxima=aggregate(filtNAB_prop_mean~Sample,data,function(x) max(x))
filtNAB_prop_mean_minima=aggregate(filtNAB_prop_mean~Sample,data,function(x) min(x))
data$filtNAB_prop_mean_scaled=apply(data,1,function(x){(as.numeric(x[33])-filtNAB_prop_mean_minima[filtNAB_prop_mean_minima$Sample==x[1],2])/(filtNAB_prop_mean_maxima[filtNAB_prop_mean_maxima$Sample==x[1],2]-filtNAB_prop_mean_minima[filtNAB_prop_mean_minima$Sample==x[1],2])})

TsTv_filtNAB_maxima=aggregate(TsTv_filtNAB~Sample,data,function(x) max(x))
TsTv_filtNAB_minima=aggregate(TsTv_filtNAB~Sample,data,function(x) min(x))
data$TsTv_filtNAB_scaled=apply(data,1,function(x){(as.numeric(x[39])-TsTv_filtNAB_minima[TsTv_filtNAB_minima$Sample==x[1],2])/(TsTv_filtNAB_maxima[TsTv_filtNAB_maxima$Sample==x[1],2]-TsTv_filtNAB_minima[TsTv_filtNAB_minima$Sample==x[1],2])})

select_variable=(data$TsTv_filtNAB_scaled+data$filtNAB_prop_mean_scaled)/2.0
data$select_variable=select_variable

lowciMean_cond_FiltNABTsTv=aggregate(select_variable~condition,data,function(x) ciMean(x,conf = .9)[1])
stquantile_cond_FiltNABTsTv=aggregate(select_variable~condition,data,function(x) quantile(x,p=.25))
stquantile_cond_FiltNABTsTv=stquantile_cond_FiltNABTsTv[order(-stquantile_cond_FiltNABTsTv$select_variable),]
lowciMean_cond_FiltNABTsTv=lowciMean_cond_FiltNABTsTv[order(-lowciMean_cond_FiltNABTsTv$select_variable),]
conditions_FiltNABTsTv_ciMean=head(lowciMean_cond_FiltNABTsTv,10)$condition
conditions_FiltNABTsTv_stquant=head(stquantile_cond_FiltNABTsTv,10)$condition
conditions_FiltNABTsTv_ciMean=factor(conditions_FiltNABTsTv_ciMean,levels=conditions_FiltNABTsTv_ciMean)
conditions_FiltNABTsTv_stquant=factor(conditions_FiltNABTsTv_stquant,levels=conditions_FiltNABTsTv_stquant)
conditions_FiltNABTsTv_ciMean %in% conditions_FiltNABTsTv_stquant ##Are they the same?

plot_sim_FiltNABTsTv_ciMean=ggplot(data=data,aes(y=filtNAB_prop_mean,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$condition %in% conditions_FiltNABTsTv_ciMean,],aes(color=factor(condition,levels=conditions_FiltNABTsTv_ciMean)),size=2)+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Somatic similarity -NAB . CI90 Mean selection")
save_plot("sim_FiltNABTsTv_ciMean.png",plot_sim_FiltNABTsTv_ciMean,base_height=10,base_aspect_ratio=1.5)
plot_sim_FiltNABTsTv_stquant=ggplot(data=data,aes(y=filtNAB_prop_mean,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$condition %in% conditions_FiltNABTsTv_stquant,],aes(color=factor(condition,levels=conditions_FiltNABTsTv_stquant)),size=2)+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Somatic similarity -NAB, Q25 selection")
save_plot("sim_FiltNABTsTv_stquant.png",plot_sim_FiltNABTsTv_stquant,base_height=10,base_aspect_ratio=1.5)

plot_var_FiltNABTsTv_ciMean=ggplot(data=data,aes(y=filtNAB_N,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$condition %in% conditions_FiltNABTsTv_ciMean,],aes(color=factor(condition,levels=conditions_FiltNABTsTv_ciMean)),size=2)+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Number of variants -NAB. CI90 Mean selection")
save_plot("var_FiltNABTsTv_ciMean.png",plot_var_FiltNABTsTv_ciMean,base_height=10,base_aspect_ratio=1.5)
plot_var_FiltNABTsTv_stquant=ggplot(data=data,aes(y=filtNAB_N,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$condition %in% conditions_FiltNABTsTv_stquant,],aes(color=factor(condition,levels=conditions_FiltNABTsTv_stquant)),size=2)+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Number of variants -NAB. Q25 selection")
save_plot("var_FiltNABTsTv_stquant.png",plot_var_FiltNABTsTv_stquant,base_height=10,base_aspect_ratio=1.5)


##Selection of filtering options attending to TsTvNAB per replicate
###################################################################

filts=aggregate(TsTv_filtNAB~Sample,data,function(x) sort(x,decreasing=TRUE)[1])

output_dataframe=data[0,]
for (i in 1:length(filts$Sample)){
	temp_dataframe=(data[data$Sample==filts[i,1],])[which(data[data$Sample==filts[i,1],]$TsTv_filtNAB==filts[filts[i,1],2]),]
	temp_dataframe=temp_dataframe[order(temp_dataframe$filtNAB_N,decreasing=TRUE),][1,]
	output_dataframe=rbind(output_dataframe,temp_dataframe)
}
ciMean(output_dataframe$filtNAB_prop_mean,conf=.9)[1]

##Independent selection of filtering options
################################################

conds=names(data)[2:(n_sim_parameters+1)]
cond_cimean=""
cond_stquant=""
for (i in 1:length(conds)){
temp_list1=aggregate(as.formula(paste0("filtNAB_prop_mean~",conds[i])),data,function(x) ciMean(x,conf = .9)[1])
temp_list2=aggregate(as.formula(paste0("filtNAB_prop_mean~",conds[i])),data,function(x) ciMean(x,conf = .9)[1])
temp_list1[which.max(temp_list1[,2]),1]
temp_list2[which.max(temp_list1[,2]),1]
cond_cimean=paste0(cond_cimean,temp_list1[which.max(temp_list1[,2]),1],",")
cond_stquant=paste0(cond_stquant,temp_list2[which.max(temp_list2[,2]),1],",")
}
cond_cimean=substr(cond_cimean, 1, nchar(cond_cimean)-1)
cond_stquant=substr(cond_stquant, 1, nchar(cond_stquant)-1)
ciMean(data[gsub(" ","",data$condition)==cond_cimean,"filtNAB_prop_mean"],conf=.9)[1]
quantile(data[gsub(" ","",data$condition)==cond_stquant,"filtNAB_prop_mean"],p=.25)


##Comparative via bootstrapping ##ATTENTION: very computationally intesive
########################################################################

total_samples=unique(data$Sample)
red_data=data[,c("Sample","filtNAB_prop_mean","condition","select_variable")]

#mybs=function(data,indices,realdata,tstvdata)
mybs=function(data,indices,realdata)
{
	my_samples=data[indices]
	#tstv_result=ciMean((tstvdata[indices,])$filtNAB_prop_mean,conf = .9)[1]
	#return(tstv_result)
	output_dataframe=realdata[0,]
	for (i in 1:length(my_samples)){
		temp_dataframe=realdata[realdata$Sample==my_samples[i],]
		temp_dataframe$Sample=i
		output_dataframe=rbind(output_dataframe,temp_dataframe)
	}
	lowciMean_cond_FiltNAB=aggregate(filtNAB_prop_mean~condition,output_dataframe,function(x) ciMean(x,conf = .9)[1])
	stquantile_cond_FiltNAB=aggregate(filtNAB_prop_mean~condition,output_dataframe,function(x) quantile(x,p=.25))
	cond_mean=lowciMean_cond_FiltNAB[which.max(lowciMean_cond_FiltNAB$filtNAB_prop_mean),1]
	cond_stquant=stquantile_cond_FiltNAB[which.max(stquantile_cond_FiltNAB$filtNAB_prop_mean),1]
	rm(lowciMean_cond_FiltNAB)
	rm(stquantile_cond_FiltNAB)
	cimean=ciMean(realdata[realdata$condition==cond_mean,"filtNAB_prop_mean"],conf=.9)[1]
	stquant=quantile(realdata[realdata$condition==cond_stquant,"filtNAB_prop_mean"],p=.25)

	lowciMean_cond_FiltNABTsTv=aggregate(select_variable~condition,output_dataframe,function(x) ciMean(x,conf = .9)[1])
	stquantile_cond_FiltNABTsTv=aggregate(select_variable~condition,output_dataframe,function(x) quantile(x,p=.25))
	cond_mean=lowciMean_cond_FiltNABTsTv[which.max(lowciMean_cond_FiltNABTsTv$select_variable),1]
	cond_stquant=stquantile_cond_FiltNABTsTv[which.max(stquantile_cond_FiltNABTsTv$select_variable),1]
	rm(lowciMean_cond_FiltNABTsTv)
	rm(stquantile_cond_FiltNABTsTv)
	citstvmean=ciMean(realdata[realdata$condition==cond_mean,"filtNAB_prop_mean"],conf=.9)[1]
	sttstvquant=quantile(realdata[realdata$condition==cond_stquant,"filtNAB_prop_mean"],p=.25)


	conds=names(realdata)[2:(n_sim_parameters+1)]
	cond_cimean=""
	cond_stquant=""
	for (i in 1:length(conds)){
	temp_list1=aggregate(as.formula(paste0("filtNAB_prop_mean~",conds[i])),output_dataframe,function(x) ciMean(x,conf = .9)[1])
	temp_list2=aggregate(as.formula(paste0("filtNAB_prop_mean~",conds[i])),output_dataframe,function(x) ciMean(x,conf = .9)[1])
	temp_list1[which.max(temp_list1[,2]),1]
	temp_list2[which.max(temp_list2[,2]),1]
	cond_cimean=paste0(cond_cimean,temp_list1[which.max(temp_list1[,2]),1],",")
	cond_stquant=paste0(cond_stquant,temp_list2[which.max(temp_list2[,2]),1],",")
	}
	cond_cimean=substr(cond_cimean, 1, nchar(cond_cimean)-1)
	cond_stquant=substr(cond_stquant, 1, nchar(cond_stquant)-1)
	ciindmean=ciMean(realdata[gsub(" ","",realdata$condition)==cond_cimean,"filtNAB_prop_mean"],conf=.9)[1]
	stindquant=quantile(realdata[gsub(" ","",realdata$condition)==cond_stquant,"filtNAB_prop_mean"],p=.25)

	return(c(cimean,stquant,citstvmean,sttstvquant,ciindmean,stindquant))
}

results_bootstrapping=boot(data=total_samples,statistic=mybs,R=104,realdata=data,parallel="multicore",ncpus=8)


#Summarized results
#output_dataframe_angelo=data[0,]
#for (i in 1:length(filts$Sample)){
#	temp_dataframe=data[data$Sample==filts[i,1],]
#	output_dataframe_angelo=rbind(output_dataframe_angelo,(temp_dataframe[order(temp_dataframe$filtNAB_prop,decreasing=TRUE),])[1:10000,])
#}


######################################################
#Best parameters per Sample
######################################################
Samples=unique(data$Sample)
sorted_filtN=data[order(data$filtN_prop_mean,decreasing=TRUE),]
sorted_filtNAB=data[order(data$filtNAB_prop_mean,decreasing=TRUE),]

####By ctbrown, obtained in stackoverflow
cbind.all <- function (...)
{
    nm <- list(...)
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow))
    do.call(cbind, lapply(nm, function(x) rbind(x, matrix(, n -
        nrow(x), ncol(x)))))
}

table_filtN=data.frame()
table_filtNAB=data.frame()

for (Sample in Samples)
{
	table_filtN=cbind.all(table_filtN,head(unique(sorted_filtN[sorted_filtN$Sample==Sample,]$conditionFilt),n=10))
	table_filtNAB=cbind.all(table_filtNAB,head(unique(sorted_filtNAB[sorted_filtNAB$Sample==Sample,]$condition),n=10))
}
colnames(table_filtN)=Samples
colnames(table_filtNAB)=Samples

sink("10best_filtN.tex")
cat("\\documentclass{article}
\\usepackage[landscape]{geometry}
\\usepackage{adjustbox}
\\begin{document}
\\begin{adjustbox}{width={\\textwidth},totalheight={\\textheight},keepaspectratio}%
\\begin{tabular}{",paste0(rep("c",ncol(table_filtN)),collapse=""),"}")
print(xtable(table_filtN),include.rownames=FALSE,only.contents=TRUE)
cat("\\end{tabular}
\\end{adjustbox}
\\end{document}")
sink()

sink("10best_filtNAB.tex")
cat("\\documentclass{article}
\\usepackage[landscape]{geometry}
\\usepackage{adjustbox}
\\begin{document}
\\begin{adjustbox}{width={\\textwidth},totalheight={\\textheight},keepaspectratio}%
\\begin{tabular}{",paste0(rep("c",ncol(table_filtNAB)),collapse=""),"}")
print(xtable(table_filtNAB),include.rownames=FALSE,only.contents=TRUE)
cat("\\end{tabular}
\\end{adjustbox}
\\end{document}")
sink()

#####################################################
#Logit regression
#####################################################

data=data[order(data$condition),]

xrepNAB=data[,1:n_sim_parameters] ##The parameters are repeated, for the six Samples.
xrepN=data[,1:(n_sim_parameters-n_sim_parameters_NAB)]

yrepfiltNAB=data$filtNAB_prop_mean
yrepfiltN=data$filtN_prop_mean

xrepfiltNAB=apply(xrepfiltNAB,2,function(x){as.numeric(gsub("-1","0",gsub("0_","",x)))})
xrepfiltNAB=apply(xrepfiltNAB,2,function(x){(x-min(x))/(max(x)-min(x))})
xrepfiltNAB=data.frame(xrepfiltNAB)
xrepfiltNAB=cbind(xrepfiltNAB,data$Sample)

xrepfiltN=apply(xrepfiltN,2,function(x){as.numeric(gsub("-1","0",gsub("0_","",x)))})
xrepfiltN=apply(xrepfiltN,2,function(x){(x-min(x))/(max(x)-min(x))})
xrepfiltN=data.frame(xrepfiltN)
xrepfiltN=cbind(xrepfiltN,data$Sample)

##Stepwise glm
###############
initial_aic_glm_filtN_model=glm(yrepfiltN~1,family=binomial(link='logit'))
aic_glm_filtN_model=step(model,as.formula(paste0("~ .+ (",paste0("xrepfiltN[,\"",colnames(xrep),"\"]",collapse="+"),")^2")))
initial_aic_glm_filtNAB_model=glm(yrepfiltNAB~1,family=binomial(link='logit'))
aic_glm_filtNAB_model=step(model,as.formula(paste0("~ .+ (",paste0("xrepfiltNAB[,\"",colnames(xrep),"\"]",collapse="+"),")^2")))

##glm full model
############################
glm_filtN_all=glm(yrepfiltN~xrepfiltN,family=binomial(link='logit'))
glm_filtNAB_all=glm(yrepfiltNAB~xrepfiltNAB,family=binomial(link='logit'))
#glm_filtN_all_interact=glm(as.formula(paste0("yrepfiltN ~  (",paste0("xrepfiltN[,\"",colnames(xrepfiltN),"\"]",collapse="+"),")^2")),family=binomial(link='logit'))
#glm_filtNAB_all_interact=glm(as.formula(paste0("yrepfiltNAB ~  (",paste0("xrepfiltNAB[,\"",colnames(xrepfiltNAB),"\"]",collapse="+"),")^2")),family=binomial(link='logit'))


#########################################################
#Regression by Sample
#########################################################


for (Sample in Samples) {
	xfiltNAB=data[data$Sample==Sample,1:n_sim_parameters] ##The parameters are repeated, for the six Samples, just keeping one
	xfiltN=data[data$Sample==Sample,1:(n_sim_parameters-n_sim_parameters_NAB)]
	xfiltNAB=apply(xfiltNAB,2,function(x){as.numeric(gsub("-1","0",gsub("0_","",x)))})
	xfiltNAB=apply(xfiltNAB,2,function(x){(x-min(x))/(max(x)-min(x))})
	xfiltN=apply(xfiltN,2,function(x){as.numeric(gsub("-1","0",gsub("0_","",x)))})
	xfiltN=apply(xfiltN,2,function(x){(x-min(x))/(max(x)-min(x))})
	#assign(paste0("glm_filtN_interact_",Sample),glm(as.formula(paste0("data[data$Sample==Sample,]$filtN_prop_mean ~  (",paste0("xfiltN[,\"",colnames(xfiltN),"\"]",collapse="+"),")^2")),family=binomial(link='logit')))
	#assign(paste0("glm_filtNAB_interact_",Sample),glm(as.formula(paste0("data[data$Sample==Sample,]$filtNAB_prop_mean ~  (",paste0("xfiltNAB[,\"",colnames(xfiltNAB),"\"]",collapse="+"),")^2")),family=binomial(link='logit')))
	assign(paste0("glm_filtN_",Sample),glm(data[data$Sample==Sample,]$filtN_prop_mean ~ xfiltN,family=binomial(link='logit')))
	assign(paste0("glm_filtNAB_",Sample),glm(data[data$Sample==Sample,]$filtNAB_prop_mean ~ xfiltNAB,family=binomial(link='logit')))
}

###Plots
########

parameters_filtNAB=names(data[,1:n_sim_parameters])
parameters_filtN=names(data[,1:(n_sim_parameters-n_sim_parameters_NAB)])

for (Sample in Samples) {
	for (parameter in parameters_filtN) {
			temp=as.numeric((data[data$Sample==Sample,])[parameter][[1]])
			x=seq(min(temp),max(temp),length.out=1000)
			myglm=eval(parse(text=paste0("glm_filtN_",Sample)))
			y=1/(1+exp(-(myglm$coefficients[1]+myglm$coefficients[paste0("xfiltN",parameter)]*x)))
			logit_data=data.frame(x,y)
			plot=ggplot(data=data[data$Sample==Sample,],aes_string(x=parameter,y="filtN_prop_mean"))+geom_violin(aes_string(group=parameter))+geom_boxplot(aes_string(group=parameter),width=.05)+geom_line(data=logit_data,aes(x=x,y=y),color="red",size=1)+labs(title=paste0(Sample,"_",parameter,"_filtN"))+scale_y_continuous(name="Similarity")
			save_plot(paste0(Sample,"_",parameter,"_filtN",".png"),plot,base_height=10,base_aspect_ratio=1.5)
		}
	for (parameter in parameters_filtNAB) {
			temp=as.numeric((data[data$Sample==Sample,])[parameter][[1]])
			x=seq(min(temp),max(temp),length.out=1000)
			myglm=eval(parse(text=paste0("glm_filtNAB_",Sample)))
			y=1/(1+exp(-(myglm$coefficients[1]+myglm$coefficients[paste0("xfiltNAB",parameter)]*x)))
			logit_data=data.frame(x,y)
			plot=ggplot(data=data[data$Sample==Sample,],aes_string(x=parameter,y="filtNAB_prop_mean"))+geom_violin(aes_string(group=parameter))+geom_boxplot(aes_string(group=parameter),width=.05)+geom_line(data=logit_data,aes(x=x,y=y),color="red",size=1)+labs(title=paste0(Sample,"_",parameter,"_filtNAB"))+scale_y_continuous(name="Similarity")
			save_plot(paste0(Sample,"_",parameter,"_filtNAB",".png"),plot,base_height=10,base_aspect_ratio=1.5)
		}
}


#########################################################
# Grouping Samples
#########################################################
res_dunn.test=dunn.test(data$filtNAB_prop_mean, g=as.factor(data$Sample)) ##Kluskal-wallis + post-hocs. It takes forever!!!!




####OLD, may be outdated
##########################

#Selection of filtering options scaling the distributions per Sample
####################################################################

maxima=data.frame(t(tapply(data$filtN_prop_mean,data$Sample,max)))
minima=data.frame(t(tapply(data$filtN_prop_mean,data$Sample,min)))
data$filtN_prop_mean_norm=apply(data,1,function(x) {as.numeric((as.numeric(x[8])-minima[as.character(x[9])])/(maxima[as.character(x[9])]-minima[as.character(x[9])]))})
lowciMean_cond_norm=aggregate(filtN_prop_mean_norm~condition,data,function(x) ciMean(x,conf = .9)[1])
stquantile_cond_norm=aggregate(filtN_prop_mean_norm~condition,data,function(x) quantile(x,p=.25))
stquantile_cond_norm=stquantile_cond_norm[order(-stquantile_cond_norm$filtN_prop_mean_norm),]
lowciMean_cond_norm=lowciMean_cond_norm[order(-lowciMean_cond_norm$filtN_prop_mean_norm),]
conditions_norm=head(lowciMean_cond_norm,10)$condition
conditions_norm2=head(stquantile_cond_norm,10)$condition
ucondnorm=union(conditions_norm,conditions_norm2)
icondnorm=intersect(conditions_norm,conditions_norm2)

#Combined: Plots with scaling
###############################

plot_norm=ggplot(data=data,aes(y=filtN_prop_mean,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$condition %in% icondnorm,],aes(color=condition),size=2)+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")
plot_norm2=ggplot(data=data,aes(y=filtN_N,x=as.factor(Sample)))+geom_violin()+geom_point(data=data[data$condition %in% icondnorm,],aes(color=condition),size=2)+scale_y_continuous(name="Variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")
save_plot("similarity_norm.png",plot_norm,base_height=10,base_aspect_ratio=1.5)
save_plot("variants_norm.png",plot_norm2,base_height=10,base_aspect_ratio=1.5)


#A and B separated
#############################################################################
dataA=data[,c(seq(1,n_sim_parameters+4),n_sim_parameters+7,n_sim_parameters+8)]
dataA$Sample=paste(dataA$Sample,"A",sep="_")
dataB=data[,c(seq(1,n_sim_parameters+2),seq(n_sim_parameters+5,n_sim_parameters+8))]
dataB$Sample=paste(dataB$Sample,"B",sep="_")

mynames=names(dataA)
mynames[n_sim_parameters+1]="filtN_N_comb"
mynames[n_sim_parameters+2]="filtN_prop_mean_comb"
mynames[n_sim_parameters+4]="filtN_N"
mynames[n_sim_parameters+3]="filtN_prop_mean"
names(dataA)=mynames
names(dataB)=mynames
datasep=rbind(dataA,dataB)

#Selection of filtering options separated
##########################################

lowciMean_cond_sep=aggregate(filtN_prop_mean~condition,datasep,function(x) ciMean(x,conf = .9)[1])
stquantile_cond_sep=aggregate(filtN_prop_mean~condition,datasep,function(x) quantile(x,p=.25))
stquantile_cond_sep=stquantile_cond_sep[order(-stquantile_cond_sep$filtN_prop_mean),]
lowciMean_cond_sep=lowciMean_cond_sep[order(-lowciMean_cond_sep$filtN_prop_mean),]
conditions_sep=head(lowciMean_cond_sep,10)$condition
conditions2_sep=head(stquantile_cond_sep,10)$condition
sort(conditions_sep)==sort(conditions2_sep) ##Different!

#By replicate: Plots without scaling
####################################

plot_byrep_combconds=ggplot(data=datasep,aes(y=filtN_prop_mean,x=as.factor(Sample)))+geom_violin()+geom_point(data=datasep[datasep$condition %in% conditions,],aes(color=condition),size=2)+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Combined filter. Similarity")
save_plot("similarity_byrep_combonds.png",plot_byrep_combconds,base_height=10,base_aspect_ratio=1.5)
plot_byrep_ciconds=ggplot(data=datasep,aes(y=filtN_prop_mean,x=as.factor(Sample)))+geom_violin()+geom_point(data=datasep[datasep$condition %in% conditions_sep,],aes(color=condition),size=2)+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="ci 90% filter. Similarity")
save_plot("similarity_byrep_ciconds.png",plot_byrep_ciconds,base_height=10,base_aspect_ratio=1.5)
plot_byrep_qconds=ggplot(data=datasep,aes(y=filtN_prop_mean,x=as.factor(Sample)))+geom_violin()+geom_point(data=datasep[datasep$condition %in% conditions2_sep,],aes(color=condition),size=2)+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="1st quantile filter. Similarity")
save_plot("similarity_byrep_qconds.png",plot_byrep_qconds,base_height=10,base_aspect_ratio=1.5)
plot_byrep_combconds_var=ggplot(data=datasep,aes(y=filtN_N,x=as.factor(Sample)))+geom_violin()+geom_point(data=datasep[datasep$condition %in% conditions,],aes(color=condition),size=2)+scale_y_continuous(name="N common variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Combined filter. Variants")
save_plot("variants_byrep_combonds.png",plot_byrep_combconds_var,base_height=10,base_aspect_ratio=1.5)
plot_byrep_ciconds_var=ggplot(data=datasep,aes(y=filtN_N,x=as.factor(Sample)))+geom_violin()+geom_point(data=datasep[datasep$condition %in% conditions_sep,],aes(color=condition),size=2)+scale_y_continuous(name="N common variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="ci 90% filter. Variants")
save_plot("variants_byrep_ciconds.png",plot_byrep_ciconds_var,base_height=10,base_aspect_ratio=1.5)
plot_byrep_qconds_var=ggplot(data=datasep,aes(y=filtN_N,x=as.factor(Sample)))+geom_violin()+geom_point(data=datasep[datasep$condition %in% conditions2_sep,],aes(color=condition),size=2)+scale_y_continuous(name="N common variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="1st quantile filter. Variants")
save_plot("variants_byrep_qconds.png",plot_byrep_qconds_var,base_height=10,base_aspect_ratio=1.5)



#################################################################################################################################################################################
#################################################################################################################################################################################
#PROBLEM DATA
#################################################################################################################################################################################
#################################################################################################################################################################################
data_problem=read.csv("results_problem.csv")

##DANGER: HARDCODED DATA!!!!
data_problem$type=c("Pure DCIS","Pure DCIS","Pure DCIS","Pure DCIS","Pure DCIS","Pure DCIS","Pure DCIS","Pure DCIS","Pure DCIS","Pure DCIS","Pure DCIS","DCIS with adjacent invasive","DCIS with adjacent invasive","DCIS with adjacent invasive","DCIS with adjacent invasive","DCIS with adjacent invasive","DCIS with adjacent invasive","DCIS with adjacent invasive","DCIS with adjacent invasive","DCIS with adjacent invasive","DCIS with adjacent invasive","DCIS with adjacent invasive")
data_problem$type=as.factor(data_problem$type)
types=levels(data_problem$type)
types=rev(types)
data_problem$filtNAB_NU_somatic=c(29,2,13,4,7,7,11,44,6,3,1,98,13,73,16,0,6,4,5,1,51,89)
data_problem$filtNAB_.U_somatic=c(31,8,18,14,20,11,20,55,17,8,6,114,43,85,36,3,20,6,15,6,62,124)

###

#Violin plot mean mutational burden
###################################

violin_mut_burden=ggplot(data=data_problem, aes(y=filtNAB_._mean,x=type,fill=type))+geom_violin()+geom_boxplot(width=.2,outlier.size=0.5)+scale_y_continuous(name=NULL)+labs(title="Number of mutations")+scale_x_discrete(name=NULL,limits=types,labels=c("Pure DCIS", "DCIS with\n adjacent invasive"))+scale_fill_manual(values=c("#009E73","#0072B2"))+guides(fill=FALSE)+theme(text=element_text(size=6,face="bold"),axis.text.y=element_text(size=4),axis.text.x=element_text(size=4,face="bold"),plot.margin=unit(c(0.2,0.1,0.1,0.1),"cm"),plot.title=element_text(size=5))
save_plot("mut_burden.png",violin_mut_burden,base_height=2,base_aspect_ratio=.5)

#Violin plot protein coding mutational burden (total)
#############################################

p_value=t.test(filtNAB_.U_somatic~type,data=data_problem)$p.value
p_valuew=wilcox.test(filtNAB_.U_somatic~type,data=data_problem)$p.value

violin_coding_burden=ggplot(data=data_problem, aes(y=filtNAB_.U_somatic,x=type,fill=type))+geom_violin()+geom_boxplot(width=.08,outlier.size=0.5)+scale_y_continuous(name=NULL)+labs(title="Number of somatic\n non-synonymous\n mutations")+scale_x_discrete(name=NULL,limits=types,labels=c("Pure DCIS", "DCIS with\n adjacent invasive"))+scale_fill_manual(values=c("#009E73","#0072B2"))+guides(fill=FALSE)+theme(text=element_text(size=6,face="bold"),axis.text.y=element_text(size=4),axis.text.x=element_text(size=4,face="bold"),plot.margin=unit(c(0.2,0.1,0.1,0.1),"cm"),plot.title=element_text(size=5))+annotate("text",x=1.2,y=100,label=paste0("T-test p-value = ",round(p_value,digits=3)),size=1.1)+annotate("text",x=1.2,y=95,label=paste0("Wilcoxon p-value =",round(p_valuew,digits=3)),size=1.1)
save_plot("coding_burden.png",violin_coding_burden,base_height=2,base_aspect_ratio=.5)


#Similarity plot
###################################

violin_similarity=ggplot(data=data_problem, aes(y=1-filtNAB_propU,x=type,fill=type))+geom_violin()+geom_boxplot(width=.2,outlier.size=0.5)+scale_y_continuous(name=NULL)+labs(title="Divergence\n between regions")+scale_x_discrete(name=NULL,limits=types,labels=c("Pure DCIS", "DCIS with\n adjacent invasive"))+scale_fill_manual(values=c("#009E73","#0072B2"))+guides(fill=FALSE)+theme(text=element_text(size=6,face="bold"),axis.text.y=element_text(size=4),axis.text.x=element_text(size=4,face="bold"),plot.margin=unit(c(0.2,0.1,0.1,0.1),"cm"),plot.title=element_text(size=5))
save_plot("similarity.png",violin_similarity,base_height=2,base_aspect_ratio=.5)


#Folds mutational burden
####################################

folds=function(x, y) {
	results=c()
	if(length(x) != length(y))
	{
		stop("The two vectors must have the same length")
	}
	for (i in 1:length(x))
	{
		if(x[i]>=y[i])
		{
			results[i]=x[i]/y[i]
		}
		else
		{
			results[i]=y[i]/x[i]
		}
	}
	return(results)
}

data_problem$foldsfiltNAB_.=folds(data_problem$AfiltNAB_.,data_problem$BfiltNAB_.)
p_value=t.test(foldsfiltNAB_.~type,data=data_problem)$p.value
p_valuew=wilcox.test(foldsfiltNAB_.~type,data=data_problem)$p.value

violin_folds_mutburden=ggplot(data=data_problem, aes(y=foldsfiltNAB_.,x=type,fill=type))+geom_violin()+geom_boxplot(width=.05,outlier.size=0.5)+scale_y_continuous(name=NULL)+labs(title="Fold-differences in \n number of mutations \n between regions")+scale_x_discrete(name=NULL,limits=types,labels=c("Pure DCIS", "DCIS with\n adjacent invasive"))+scale_fill_manual(values=c("#009E73","#0072B2"))+guides(fill=FALSE)+theme(text=element_text(size=6,face="bold"),axis.text.y=element_text(size=4),axis.text.x=element_text(size=4,face="bold"),plot.margin=unit(c(0.2,0.1,0.1,0.1),"cm"),plot.title=element_text(size=5))+annotate("text",x=1.2,y=6,label=paste0("T-test p-value = ",round(p_value,digits=3)),size=1.1)+annotate("text",x=1.2,y=5.7,label=paste0("Wilcoxon p-value = ",round(p_valuew,digits=3)),size=1.1)
save_plot("folds.png",violin_folds_mutburden,base_height=2,base_aspect_ratio=.5)
