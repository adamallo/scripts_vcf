library(cowplot)
library(data.table)

thisdata=fread("/Users/Diego/Desktop/dcis/replicates_withPAF/resultsSensitivity.csv")
theplot=ggplot(thisdata,aes(x=gsize,y=rAc))+geom_point()+scale_x_continuous(name="Technical replicates",breaks=scales::pretty_breaks(n = 10))+stat_summary(fun.y="mean",colour="red",geom="point",shape=3,size=2,stroke=1.5)+scale_y_continuous(name="Score, relative to best using all samples")
save_plot(theplot,file="/Users/Diego/Desktop/dcis/replicates_withPAF/sensitivity.pdf",base_height = 6)
