library(lsr)
library(cowplot)
library(xtable)
library(dunn.test)

# Data parsing #HARDCODED
############################################

###Conf variables
dir="~/res3_tab"
n_sim_parameters=9 ###WARNING: CHANGE THIS ACCORDINGLY WITH YOUR DATA. TOTAL NUMBER OF PARAMETERS
n_sim_parameters_NAB=3 ###WARNING: CHANGE THIS ACCORDINGLY WITH YOUR DATA. PARAMETERS RELATIVE TO NAB
setwd(dir)
#################

B3=read.csv("B3.csv")
B6=read.csv("B6.csv")
D5=read.csv("D5.csv")
D8=read.csv("D8.csv")
DCIS64=read.csv("DCIS64.csv")
K12=read.csv("K12.csv")
B3$sample=rep("B3",nrow(B3))
B6$sample=rep("B6",nrow(B6))
D5$sample=rep("D5",nrow(D5))
D8$sample=rep("D8",nrow(D8))
DCIS64$sample=rep("DCIS64",nrow(DCIS64))
K12$sample=rep("K12",nrow(K12))


data=rbind(B3,B6,D5,D8,DCIS64,K12)
data$condition=apply(data,1,function(x) paste(x[1:n_sim_parameters],collapse=","))
data$conditionFilt=apply(data,1,function(x) paste(x[1:(n_sim_parameters-n_sim_parameters_NAB)],collapse=","))
data$conditionNAB=apply(data,1,function(x) paste(x[(n_sim_parameters-n_sim_parameters_NAB+1):n_sim_parameters],collapse=","))


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
plot_sim_FiltN_ciMean=ggplot(data=data,aes(y=filtN_prop_mean,x=as.factor(sample)))+geom_violin()+geom_point(data=data[data$conditionFilt %in% conditions_FiltN_ciMean,],aes(color=factor(conditionFilt,levels=conditions_FiltN_ciMean)),size=2,stat="summary",fun.y="mean")+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Somatic similarity. CI90 Mean selection")
save_plot("sim_FiltN_ciMean.png",plot_sim_FiltN_ciMean,base_height=10,base_aspect_ratio=1.5)
plot_sim_FiltN_stquant=ggplot(data=data,aes(y=filtN_prop_mean,x=as.factor(sample)))+geom_violin()+geom_point(data=data[data$conditionFilt %in% conditions_FiltN_stquant,],aes(color=factor(conditionFilt,levels=conditions_FiltN_stquant)),size=2,stat="summary",fun.y="mean")+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Somatic similarity. Q25 selection")
save_plot("sim_FiltN_stquant.png",plot_sim_FiltN_stquant,base_height=10,base_aspect_ratio=1.5)

plot_var_FiltN_ciMean=ggplot(data=data,aes(y=filtN_N,x=as.factor(sample)))+geom_violin()+geom_point(data=data[data$conditionFilt %in% conditions_FiltN_ciMean,],aes(color=factor(conditionFilt,levels=conditions_FiltN_ciMean)),size=2,stat="summary",fun.y="mean")+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Number of variants. CI90 Mean selection")
save_plot("var_FiltN_ciMean.png",plot_sim_FiltN_ciMean,base_height=10,base_aspect_ratio=1.5)
plot_var_FiltN_stquant=ggplot(data=data,aes(y=filtN_N,x=as.factor(sample)))+geom_violin()+geom_point(data=data[data$conditionFilt %in% conditions_FiltN_stquant,],aes(color=factor(conditionFilt,levels=conditions_FiltN_stquant)),size=2,stat="summary",fun.y="mean")+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Number of variants. Q25 selection")
save_plot("var_FiltN_stquant.png",plot_sim_FiltN_stquant,base_height=10,base_aspect_ratio=1.5)

#FiltNAB

plot_sim_FiltNAB_ciMean=ggplot(data=data,aes(y=filtNAB_prop_mean,x=as.factor(sample)))+geom_violin()+geom_point(data=data[data$condition %in% conditions_FiltNAB_ciMean,],aes(color=factor(condition,levels=conditions_FiltNAB_ciMean)),size=2)+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Somatic similarity -NAB . CI90 Mean selection")
save_plot("sim_FiltNAB_ciMean.png",plot_sim_FiltNAB_ciMean,base_height=10,base_aspect_ratio=1.5)
plot_sim_FiltNAB_stquant=ggplot(data=data,aes(y=filtNAB_prop_mean,x=as.factor(sample)))+geom_violin()+geom_point(data=data[data$condition %in% conditions_FiltNAB_stquant,],aes(color=factor(condition,levels=conditions_FiltNAB_stquant)),size=2)+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Somatic similarity -NAB. Q25 selection")
save_plot("sim_FiltNAB_stquant.png",plot_sim_FiltNAB_stquant,base_height=10,base_aspect_ratio=1.5)

plot_var_FiltNAB_ciMean=ggplot(data=data,aes(y=filtNAB_N,x=as.factor(sample)))+geom_violin()+geom_point(data=data[data$condition %in% conditions_FiltNAB_ciMean,],aes(color=factor(condition,levels=conditions_FiltNAB_ciMean)),size=2)+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Number of variants -NAB. CI90 Mean selection")
save_plot("var_FiltNAB_ciMean.png",plot_sim_FiltNAB_ciMean,base_height=10,base_aspect_ratio=1.5)
plot_var_FiltNAB_stquant=ggplot(data=data,aes(y=filtNAB_N,x=as.factor(sample)))+geom_violin()+geom_point(data=data[data$condition %in% conditions_FiltNAB_stquant,],aes(color=factor(condition,levels=conditions_FiltNAB_stquant)),size=2)+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Number of variants -NAB. Q25 selection")
save_plot("var_FiltNAB_stquant.png",plot_sim_FiltNAB_stquant,base_height=10,base_aspect_ratio=1.5)

#NAB

plot_sim_NAB_ciMean=ggplot(data=data,aes(y=filtNAB_prop_mean,x=as.factor(sample)))+geom_violin()+geom_point(data=data[data$conditionNAB %in% conditionsNAB_FiltNAB_ciMean,],aes(color=factor(conditionNAB,conditionsNAB_FiltNAB_ciMean)),size=2,stat="summary",fun.y="mean")+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Mean somatic similarity, NAB filter. CI90 Mean selection")
save_plot("sim_NAB_ciMean.png",plot_sim_NAB_ciMean,base_height=10,base_aspect_ratio=1.5)
plot_sim_NAB_stquant=ggplot(data=data,aes(y=filtNAB_prop_mean,x=as.factor(sample)))+geom_violin()+geom_point(data=data[data$conditionNAB %in% conditionsNAB_FiltNAB_stquant,],aes(color=factor(conditionNAB,conditionsNAB_FiltNAB_stquant)),size=2,stat="summary",fun.y="mean")+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Mean somatic similarity, NAB filter. Q25 selection")
save_plot("sim_NAB_stquant.png",plot_sim_NAB_stquant,base_height=10,base_aspect_ratio=1.5)

plot_var_NAB_ciMean=ggplot(data=data,aes(y=filtNAB_N,x=as.factor(sample)))+geom_violin()+geom_point(data=data[data$conditionNAB %in% conditionsNAB_FiltNAB_ciMean,],aes(color=factor(conditionNAB,conditionsNAB_FiltNAB_ciMean)),size=2,stat="summary",fun.y="mean")+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Mean number of variants, NAB filter. CI90 Mean selection")
save_plot("var_NAB_ciMean.png",plot_sim_NAB_ciMean,base_height=10,base_aspect_ratio=1.5)
plot_var_NAB_stquant=ggplot(data=data,aes(y=filtNAB_N,x=as.factor(sample)))+geom_violin()+geom_point(data=data[data$conditionNAB %in% conditionsNAB_FiltNAB_stquant,],aes(color=factor(conditionNAB,conditionsNAB_FiltNAB_stquant)),size=2,stat="summary",fun.y="mean")+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Mean number of variants, NAB filter. Q25 selection")
save_plot("var_NAB_stquant.png",plot_sim_NAB_stquant,base_height=10,base_aspect_ratio=1.5)

######################################################
#Best parameters per Sample
######################################################
samples=unique(data$sample)
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

for (sample in samples)
{
	table_filtN=cbind.all(table_filtN,head(unique(sorted_filtN[sorted_filtN$sample==sample,]$conditionFilt),n=10))
	table_filtNAB=cbind.all(table_filtNAB,head(unique(sorted_filtNAB[sorted_filtNAB$sample==sample,]$condition),n=10))
}
colnames(table_filtN)=samples
colnames(table_filtNAB)=samples

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

xrepNAB=data[,1:n_sim_parameters] ##The parameters are repeated, for the six samples.
xrepN=data[,1:(n_sim_parameters-n_sim_parameters_NAB)]

yrepfiltNAB=data$filtNAB_prop_mean
yrepfiltN=data$filtN_prop_mean

xrepfiltNAB=apply(xrepfiltNAB,2,function(x){as.numeric(gsub("-1","0",gsub("0_","",x)))})
xrepfiltNAB=apply(xrepfiltNAB,2,function(x){(x-min(x))/(max(x)-min(x))})
xrepfiltNAB=data.frame(xrepfiltNAB)
xrepfiltNAB=cbind(xrepfiltNAB,data$sample)

xrepfiltN=apply(xrepfiltN,2,function(x){as.numeric(gsub("-1","0",gsub("0_","",x)))})
xrepfiltN=apply(xrepfiltN,2,function(x){(x-min(x))/(max(x)-min(x))})
xrepfiltN=data.frame(xrepfiltN)
xrepfiltN=cbind(xrepfiltN,data$sample)

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
#Regression by sample
#########################################################


for (sample in samples) {
	xfiltNAB=data[data$sample==sample,1:n_sim_parameters] ##The parameters are repeated, for the six samples, just keeping one
	xfiltN=data[data$sample==sample,1:(n_sim_parameters-n_sim_parameters_NAB)]
	xfiltNAB=apply(xfiltNAB,2,function(x){as.numeric(gsub("-1","0",gsub("0_","",x)))})
	xfiltNAB=apply(xfiltNAB,2,function(x){(x-min(x))/(max(x)-min(x))})
	xfiltN=apply(xfiltN,2,function(x){as.numeric(gsub("-1","0",gsub("0_","",x)))})
	xfiltN=apply(xfiltN,2,function(x){(x-min(x))/(max(x)-min(x))})
	#assign(paste0("glm_filtN_interact_",sample),glm(as.formula(paste0("data[data$sample==sample,]$filtN_prop_mean ~  (",paste0("xfiltN[,\"",colnames(xfiltN),"\"]",collapse="+"),")^2")),family=binomial(link='logit')))
	#assign(paste0("glm_filtNAB_interact_",sample),glm(as.formula(paste0("data[data$sample==sample,]$filtNAB_prop_mean ~  (",paste0("xfiltNAB[,\"",colnames(xfiltNAB),"\"]",collapse="+"),")^2")),family=binomial(link='logit')))
	assign(paste0("glm_filtN_",sample),glm(data[data$sample==sample,]$filtN_prop_mean ~ xfiltN,family=binomial(link='logit')))
	assign(paste0("glm_filtNAB_",sample),glm(data[data$sample==sample,]$filtNAB_prop_mean ~ xfiltNAB,family=binomial(link='logit')))
}

###Plots
########

parameters_filtNAB=names(data[,1:n_sim_parameters])
parameters_filtN=names(data[,1:(n_sim_parameters-n_sim_parameters_NAB)])

for (sample in samples) {
	for (parameter in parameters_filtN) {
			temp=as.numeric((data[data$sample==sample,])[parameter][[1]])
			x=seq(min(temp),max(temp),length.out=1000)
			myglm=eval(parse(text=paste0("glm_filtN_",sample)))
			y=1/(1+exp(-(myglm$coefficients[1]+myglm$coefficients[paste0("xfiltN",parameter)]*x)))
			logit_data=data.frame(x,y)
			plot=ggplot(data=data[data$sample==sample,],aes_string(x=parameter,y="filtN_prop_mean"))+geom_violin(aes_string(group=parameter))+geom_boxplot(aes_string(group=parameter),width=.05)+geom_line(data=logit_data,aes(x=x,y=y),color="red",size=1)+labs(title=paste0(sample,"_",parameter,"_filtN"))+scale_y_continuous(name="Similarity")
			save_plot(paste0(sample,"_",parameter,"_filtN",".png"),plot,base_height=10,base_aspect_ratio=1.5)
		}
	for (parameter in parameters_filtNAB) {
			temp=as.numeric((data[data$sample==sample,])[parameter][[1]])
			x=seq(min(temp),max(temp),length.out=1000)
			myglm=eval(parse(text=paste0("glm_filtNAB_",sample)))
			y=1/(1+exp(-(myglm$coefficients[1]+myglm$coefficients[paste0("xfiltNAB",parameter)]*x)))
			logit_data=data.frame(x,y)
			plot=ggplot(data=data[data$sample==sample,],aes_string(x=parameter,y="filtNAB_prop_mean"))+geom_violin(aes_string(group=parameter))+geom_boxplot(aes_string(group=parameter),width=.05)+geom_line(data=logit_data,aes(x=x,y=y),color="red",size=1)+labs(title=paste0(sample,"_",parameter,"_filtNAB"))+scale_y_continuous(name="Similarity")
			save_plot(paste0(sample,"_",parameter,"_filtNAB",".png"),plot,base_height=10,base_aspect_ratio=1.5)
		}
}


#########################################################
# Grouping samples
#########################################################
res_dunn.test=dunn.test(data$filtNAB_prop_mean, g=as.factor(data$sample)) ##Kluskal-wallis + post-hocs. It takes forever!!!!




####OLD, may be outdated
##########################

#Selection of filtering options scaling the distributions per sample
####################################################################

maxima=data.frame(t(tapply(data$filtN_prop_mean,data$sample,max)))
minima=data.frame(t(tapply(data$filtN_prop_mean,data$sample,min)))
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

plot_norm=ggplot(data=data,aes(y=filtN_prop_mean,x=as.factor(sample)))+geom_violin()+geom_point(data=data[data$condition %in% icondnorm,],aes(color=condition),size=2)+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")
plot_norm2=ggplot(data=data,aes(y=filtN_N,x=as.factor(sample)))+geom_violin()+geom_point(data=data[data$condition %in% icondnorm,],aes(color=condition),size=2)+scale_y_continuous(name="Variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")
save_plot("similarity_norm.png",plot_norm,base_height=10,base_aspect_ratio=1.5)
save_plot("variants_norm.png",plot_norm2,base_height=10,base_aspect_ratio=1.5)


#A and B separated
#############################################################################
dataA=data[,c(seq(1,n_sim_parameters+4),n_sim_parameters+7,n_sim_parameters+8)]
dataA$sample=paste(dataA$sample,"A",sep="_")
dataB=data[,c(seq(1,n_sim_parameters+2),seq(n_sim_parameters+5,n_sim_parameters+8))]
dataB$sample=paste(dataB$sample,"B",sep="_")

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

plot_byrep_combconds=ggplot(data=datasep,aes(y=filtN_prop_mean,x=as.factor(sample)))+geom_violin()+geom_point(data=datasep[datasep$condition %in% conditions,],aes(color=condition),size=2)+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Combined filter. Similarity")
save_plot("similarity_byrep_combonds.png",plot_byrep_combconds,base_height=10,base_aspect_ratio=1.5)
plot_byrep_ciconds=ggplot(data=datasep,aes(y=filtN_prop_mean,x=as.factor(sample)))+geom_violin()+geom_point(data=datasep[datasep$condition %in% conditions_sep,],aes(color=condition),size=2)+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="ci 90% filter. Similarity")
save_plot("similarity_byrep_ciconds.png",plot_byrep_ciconds,base_height=10,base_aspect_ratio=1.5)
plot_byrep_qconds=ggplot(data=datasep,aes(y=filtN_prop_mean,x=as.factor(sample)))+geom_violin()+geom_point(data=datasep[datasep$condition %in% conditions2_sep,],aes(color=condition),size=2)+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="1st quantile filter. Similarity")
save_plot("similarity_byrep_qconds.png",plot_byrep_qconds,base_height=10,base_aspect_ratio=1.5)
plot_byrep_combconds_var=ggplot(data=datasep,aes(y=filtN_N,x=as.factor(sample)))+geom_violin()+geom_point(data=datasep[datasep$condition %in% conditions,],aes(color=condition),size=2)+scale_y_continuous(name="N common variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="Combined filter. Variants")
save_plot("variants_byrep_combonds.png",plot_byrep_combconds_var,base_height=10,base_aspect_ratio=1.5)
plot_byrep_ciconds_var=ggplot(data=datasep,aes(y=filtN_N,x=as.factor(sample)))+geom_violin()+geom_point(data=datasep[datasep$condition %in% conditions_sep,],aes(color=condition),size=2)+scale_y_continuous(name="N common variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="ci 90% filter. Variants")
save_plot("variants_byrep_ciconds.png",plot_byrep_ciconds_var,base_height=10,base_aspect_ratio=1.5)
plot_byrep_qconds_var=ggplot(data=datasep,aes(y=filtN_N,x=as.factor(sample)))+geom_violin()+geom_point(data=datasep[datasep$condition %in% conditions2_sep,],aes(color=condition),size=2)+scale_y_continuous(name="N common variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")+labs(title="1st quantile filter. Variants")
save_plot("variants_byrep_qconds.png",plot_byrep_qconds_var,base_height=10,base_aspect_ratio=1.5)