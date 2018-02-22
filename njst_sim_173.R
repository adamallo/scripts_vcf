library(phangorn)
sim_data=read.csv("Desktop/dcis/similarity_dcis173_N.csv")
dist_data=1-sim_data
dist_m=as.matrix(dist_data)
rownames(dist_m)=colnames(dist_m)
tree=nj(dist_m)
write.tree(phy = tree,file = "Desktop/dcis/nj_similarity_dcis173_N.tree")