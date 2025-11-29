#install.packages("devtools")
#devtools::install_github("igrabski/sc-SHC")
library(scSHC)
setwd("C:/Users/Administrator/Desktop/pbmc3k/")
# read.label
table = read.table("S4_Layer7_Annotation.csv",header = T,sep=",",row.names = 1)
clust = table[,2001:2015]
# read.counts
count = read.table("pbmc3k_raw.txt",header=T,sep="\t",row.names = 1)
data  = as.matrix(count)




new_clusters = testClusters(data, as.character(clust$L7AMap), num_features = 2000, num_PCs = 30, parallel = F, cores=1)
new_labels = new_clusters[[1]]
new_labels = data.frame(new_labels)
write.table(new_labels, "pbmc_new_label.txt",row.names = F,quote=F,sep="\t")



