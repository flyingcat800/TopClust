#install.packages("clustree")
library(clustree)

setwd("C:/Users/Administrator/Desktop/pbmc3k/")
data = read.table("pbmc3k_layers_annotation.csv",sep=",", row.names = 1, header=T)

png("PBMC_3K_Layers_subclusters.png", width = 1920*24, height = 1080*6, res=600)
clustree(data, prefix = "Layer")
dev.off()