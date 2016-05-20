library(ltr)
library(cowplot)

# Data parsing #HARDCODED
############################################

###Conf variables
dir="~/res_angelo_tabulated2"
n_sim_parameters=6
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

#Selection of filtering options
###############################

lowciMean_cond=aggregate(filtN_prop_mean~condition,data,function(x) ciMean(x,conf = .9)[1])
stquantile_cond=aggregate(filtN_prop_mean~condition,data,function(x) quantile(x,p=.25))
stquantile_cond=stquantile_cond[order(-stquantile_cond$filtN_prop_mean),]
lowciMean_cond=lowciMean_cond[order(-lowciMean_cond$filtN_prop_mean),]
conditions=head(lowciMean_cond,10)$condition
conditions2=head(stquantile_cond,10)$condition
sort(conditions)==sort(conditions2) ##They all are the same

lowciMean_cond_sep=aggregate(filtN_prop_mean~condition,datasep,function(x) ciMean(x,conf = .9)[1])
stquantile_cond_sep=aggregate(filtN_prop_mean~condition,datasep,function(x) quantile(x,p=.25))
stquantile_cond_sep=stquantile_cond_sep[order(-stquantile_cond_sep$filtN_prop_mean),]
lowciMean_cond_sep=lowciMean_cond_sep[order(-lowciMean_cond_sep$filtN_prop_mean),]
conditions_sep=head(lowciMean_cond_sep,10)$condition
conditions2_sep=head(stquantile_cond_sep,10)$condition
sort(conditions_sep)==sort(conditions2_sep) ##Different!

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
icondnorm=intersect(conditions_norm,conditions_norm2

#Combined: Plots without scaling
#################################

plot=ggplot(data=data,aes(y=filtN_prop_mean,x=as.factor(sample)))+geom_violin()+geom_point(data=data[data$condition %in% conditions,],aes(color=condition),size=2)+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")
save_plot("similarity.png",plot,base_height=10,base_aspect_ratio=1.5)
plot2=ggplot(data=data,aes(y=filtN_N,x=as.factor(sample)))+geom_violin()+geom_point(data=data[data$condition %in% conditions,],aes(color=condition),size=2)+scale_y_continuous(name="Number of variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")
save_plot("variants.png",plot2,base_height=10,base_aspect_ratio=1.5)

#Combined: Plots with scaling
###############################

plot_norm=ggplot(data=data,aes(y=filtN_prop_mean,x=as.factor(sample)))+geom_violin()+geom_point(data=data[data$condition %in% icondnorm,],aes(color=condition),size=2)+scale_y_continuous(name="Similarity")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")
plot_norm2=ggplot(data=data,aes(y=filtN_N,x=as.factor(sample)))+geom_violin()+geom_point(data=data[data$condition %in% icondnorm,],aes(color=condition),size=2)+scale_y_continuous(name="Variants")+scale_x_discrete(name="Sample")+scale_color_discrete(name="Filter")
save_plot("similarity_norm.png",plot_norm,base_height=10,base_aspect_ratio=1.5)
save_plot("variants_norm.png",plot_norm2,base_height=10,base_aspect_ratio=1.5)

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