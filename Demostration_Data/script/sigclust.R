library(sigclust)
setwd("C:/Users/Administrator/Desktop/pbmc3k/")

dat = read.table("S4_Layer2_Annotation.csv", header = T, sep = ",", row.names = 1)
#Layer1

label = c(1,2,3,4,5,6,7,8)
Len = length(label)
for(i in 1:Len)
{
  for(j in i:Len)
  {    
    if(i!=j)
    {
      temp = dat[dat$L1P == label[i] | dat$L1P == label[j] , ]
      temp = temp[,1:2000]
      pvalue = sigclust(temp, nsim = 1000, nrep=1, labflag=0, label=0, icovest=1)
      png(paste(c("sigclust_layer1_cluster_",i,"-",j,".png"),collapse = ""), width = 1920,height = 1080)
      plot(pvalue)
      dev.off()
    }
  }
}



label = unique(dat$L2PMap)
clean = c()
for(i in label)
{
  if(substr(i,nchar(i)-1,nchar(i)) != "_0")
  {
    clean = c(clean,i)
  }
}

Len = length(clean)
for(i in 1:Len)
{
  for(j in i:Len)
  {    
    if(substr(clean[i],1,nchar(clean[i])-2) == substr(clean[j],1,nchar(clean[j])-2))
    {
      a = dat[dat$L2PMap == clean[i] | dat$L2PMap == clean[j] , ]
      a = a[,1:2000]
      if(nrow(a)>0)
      {
        print(paste(c(clean[i],"_",clean[j])))
        pvalue = sigclust(a, nsim = 1000, nrep=1, labflag=0, label=0, icovest=1)
        png(paste(c("sigclust_layer2_subcluster",clean[i],"-",clean[j],".png"),collapse = ""), width = 1920,height = 1080)
        plot(pvalue)
        dev.off()
      }

    }
  }
}








label = unique(dat$L2PMap)
clean = c()
for(i in label)
{
  if(substr(i,nchar(i)-1,nchar(i)) != "_0")
  {
    clean = c(clean,i)
  }
}
clean = sort(clean)
ids = dat$L2PMap


for(i in c(1,2,3,4,5,6,7,8))
{
  pattern = paste(c(i,"_"),collapse = "")
  subset = grep("1_",ids)
  a = dat[subset , ]
  a = a[,1:2000]
  print(i)
  pvalue = sigclust(a, nsim = 1000, nrep=1, labflag=0, label=0, icovest=1)
  png(paste(c("sigclust_layer2_subclusters_of_Cluster",i,".png"),collapse = ""), width = 1920,height = 1080)
  plot(pvalue)
  dev.off()
}