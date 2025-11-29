#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 17:24:28 2022

@author: cshu
"""

### Step-0 Import packages
import os
os.chdir("C:/Users/Administrator/Desktop/pbmc3k/")
from script.functions import *


### Step-1 Specify dataset
starttime = datetime.datetime.now()
file = "pbmc3K_raw.txt"
name = file.split("_")[0]

if not os.path.exists(name):
    os.makedirs(name)
    
### Step-2 Data processing
S1 = pd.read_csv(file, header=0, index_col=0, sep="\t")
S2 = Quality_Control(S1, cellNum=0, UMINumLowLimit= 500,UMINumUpLimit=10000,mtPct=0.1)
S3 = Normalize_Rawdata(S2)
S4 = Select_HEGenes(S3,ngenes=2000)
S5 = Scaling_Data(S4)
S6 = S5.T   

truelabel = pd.read_csv(name+'_label.txt', header=0, index_col=0, sep="\t")
truelabel = truelabel["pbmc.active.ident"]


### Step-3 PCA dimension reduction for density peaks finding
Dimension = 100
pca = PCA(n_components=Dimension,  random_state = 1) 
PCA_s = pd.DataFrame(pca.fit_transform(S6), index=S6.index, columns=["PC%d" % k for k in range(1,Dimension+1)]).iloc[:,:Dimension]
ElbowPlot(S6,name)
Dim = PCA_Permutation(S6, name, nperm = 100, nPC = 40)


PCA_Data = np.round(PCA_s.iloc[:,0:Dim],4)
CellID = pd.DataFrame(PCA_Data.index, columns=["ID"])
nCells = PCA_Data.shape[0]


### Step-4 tSNE dimension reduction for visualization
# tsne = TSNE(n_components=2, random_state=1).fit_transform(PCA_Data)
# tSNE = pd.DataFrame(tsne,columns=["tSNE-1","tSNE-2"]) 
tSNE = pd.read_csv("tSNE_Ref.csv", header=0, index_col=0, sep=",")
Scatter_Plot(tSNE,name+"/"+name+"_tSNE","tSNE-1","tSNE-2")
endtime = datetime.datetime.now()
print("\nPreprocessing Time:")
print((endtime-starttime))  
tSNEPlot(tSNE,name +"/Label_of_Reference",truelabel,legend={"loc":'center right', "bbox_to_anchor":(1.3,0.5),"fontsize":5, "ncol":1})
### Step-5 Node-weighted SNN Graph
starttime = datetime.datetime.now()
K = K_Estimation(PCA_Data,5)
print("K = " + str(K))


endtime = datetime.datetime.now()
print("\nNeighborhood Size Esitmation Time:")
print((endtime-starttime))  

# the least K that makes SNN a 1-connected graph is estimated at 13
starttime = datetime.datetime.now()
Edges, Den = KNN_Density_Ball_SNN(PCA_Data, K)
#Edges, Den = SNN_Graph(PCA_Data, K)
endtime = datetime.datetime.now()
print("\nSNN network construction Time:")
print((endtime-starttime))  


### Step-6 Persistenct homology on superlevel set filtration
starttime = datetime.datetime.now()
PH, Cluster, Cluster_U = PersistentHomology(PCA_Data, tSNE, Edges, Den, name, pd_name = "raw", Stratified = True, iter_plot = False, cluster_plot = False)   
endtime = datetime.datetime.now()
print("\nPeaks Identification Time:")
print((endtime-starttime))  

### Step-7 Optional: Cluster stability evaluation via downsampling 
# starttime = datetime.datetime.now()
# Summary = DownSamplingStability(PCA_Data, PH, Cluster, CellID, tSNE, K, name, iteration = 20)
# candidatePeaks = Summary.loc[Summary.Stable_Fraction >= 0,"Barcode"].values
# endtime = datetime.datetime.now()
# print("\nPeaks Identification Time:")
# print((endtime-starttime))  


### Step-8 Report peaks_fraction of each cluster
starttime = datetime.datetime.now()
PH_Peaks, Cluster_Peaks = Peaks_Cells_Label(PH, Cluster, Cluster_U)
uniq_label = list(np.unique(Cluster_Peaks.Label.values))
endtime = datetime.datetime.now()
print("\nPeaks Identification Time:")
print((endtime-starttime))  

### Step-9 Rename peaks
rename_dict = dict(zip(uniq_label, np.arange(len(uniq_label))))
Cluster_Peaks.loc[:,"Label"] =  [rename_dict[i] for i in Cluster_Peaks.Label.values]


### Step-9 Iterative KNN to call all cells back
Cluster_All = Iterative_KNN_Calling(Cluster_Peaks, PCA_Data, tSNE, name, iter_plot=False)
nClust = len(list(np.unique(Cluster_All.Label.values)))

predictlabel = Cluster_All.Label.values
ARI = np.round(adjusted_rand_score(truelabel,predictlabel),3)
endtime = datetime.datetime.now()
print("\nAssign All Other Cells Time:")
print((endtime-starttime))

### Step-10 Output results
starttime = datetime.datetime.now()
# report parameters
g=open(name + "/summary.txt","w")
g.write("The estimated number of PCs is: " + str(Dim) + "\n")
g.write("The estimated neighborhood is: " + str(K)+ "\n")
g.write("The estimated cluster number is: " +  str(nClust)+ "\n")
g.write("The ARI index is: " +  str(ARI)+ "\n")
g.close()

# output PCA table with clustering labels
#PCA_Data["RealLabel"] = truelabel
PCA_Data["PeakLabel"] = Cluster_Peaks.Label.values
PCA_Data["AllLabel"] = Cluster_All.Label.values
PCA_Data["TrueLabel"] = truelabel.values
PCA_Data.to_csv(name + "/PCA.csv")

# output normalized gene expression table with clustering labels
S4 = S4.T
#S4["RealLabel"] = truelabel
S4["PeakLabel"] = Cluster_Peaks.Label.values
S4["AllLabel"]  = Cluster_All.Label.values
S4["TrueLabel"] = truelabel.values
S4.to_csv(name + "/S4.csv")

# output other important results
tSNE["PeakLabel"] = Cluster_Peaks.Label.values
tSNE["AllLabel"]  = Cluster_All.Label.values
tSNE.to_csv(name + "/tSNE.csv")

# output clustering results
PH.to_csv(name + "/PersistenceDiagram.csv")
PH_Peaks.to_csv(name + "/PersistenceDiagramPeaks.csv")
Cluster_All.to_csv(name + "/ClusteringLabel.csv")
Cluster_Peaks.to_csv(name + "/ClusteringLabelPeaks.csv")

### Step-11 Further analysis
# intra-cluster stability(pair-wise pearson correlation)
PCA_Data = pd.read_csv(name + "/PCA.csv", header=0, index_col=0, sep=",")
IntraClusterCorrelation(PCA_Data, "AllLabel", name)
IntraClusterCorrelation(PCA_Data, "PeakLabel", name)
IntraClusterCorrelation(PCA_Data, "TrueLabel", name)
PCA_Data["MergedLabel"] = PCA_Data["TrueLabel"]
PCA_Data.loc[PCA_Data.MergedLabel == "CD14+ Mono","MergedLabel"] = "Merged_7"
PCA_Data.loc[PCA_Data.MergedLabel == "DC","MergedLabel"] = "Merged_7"
IntraClusterCorrelation(PCA_Data, "MergedLabel", name)
Markers = ClusterSpecificGeneSet(S4, name, top=10)
#Markers_all = ClusterSpecificGeneSetAllCell(S4, name, top=10)

PeakLabel = Cluster_Peaks.Label.values
S4["PeakLabel"] = Cluster_Peaks.Label.values
S4["AllLabel"]  = Cluster_All.Label.values
S4.to_csv(name + "/S4.csv")

tSNE.to_csv(name + "/tSNE.csv")
PH.to_csv(name + "/PersistenceDiagram.csv")
PH_Peaks.to_csv(name + "/PersistenceDiagramPeaks.csv")
Cluster_All.to_csv(name + "/ClusteringLabel.csv")
PredictLabel = Cluster_All.Label.values
LabelID = np.unique(PredictLabel)
peak_tpm = RepresentativeTranscriptome2(S4,"PeakLabel")
seurat_tpm = RepresentativeTranscriptome2(S4,"TrueLabel")



# 对 cluster 1 进行细分
name="cluster1"
if not os.path.exists(name):
    os.makedirs(name)
subset = S4.loc[S4.AllLabel==1,:]
S1_sub = S1.loc[:,subset.index]
S2_sub = Quality_Control(S1_sub, cellNum=0, UMINumLowLimit= 500,UMINumUpLimit=10000,mtPct=0.1)
S3_sub = Normalize_Rawdata(S2_sub)
S4_sub = Select_HEGenes(S3_sub,ngenes=2000)
S5_sub = Scaling_Data(S4_sub)
S6_sub = S5_sub.T   


### Step-3 PCA dimension reduction for density peaks finding
Dimension = 100
pca = PCA(n_components=Dimension,  random_state = 1) 
PCA_s_sub = pd.DataFrame(pca.fit_transform(S6_sub), index=S6_sub.index, columns=["PC%d" % k for k in range(1,Dimension+1)]).iloc[:,:Dimension]
ElbowPlot(S6_sub,name)
Dim_sub = PCA_Permutation(S6_sub, name, nperm = 100, nPC = 40)
#Dim=7

PCA_Data_sub = np.round(PCA_s_sub.iloc[:,0:Dim_sub],4)
CellID_sub = pd.DataFrame(PCA_Data_sub.index, columns=["ID"])
nCells_sub = PCA_Data_sub.shape[0]


### Step-4 tSNE dimension reduction for visualization
tSNE_sub = tSNE.loc[tSNE.AllLabel==1,:]
tSNE_sub.index = PCA_Data_sub.index
Scatter_Plot(tSNE_sub,name+"/"+name+"_tSNE","tSNE-1","tSNE-2")
endtime = datetime.datetime.now()
print("\nPreprocessing Time:")
print((endtime-starttime))  


### Step-5 Node-weighted SNN Graph
starttime = datetime.datetime.now()
K = K_Estimation(PCA_Data_sub,5)
print("K = " + str(K))
#K=9

endtime = datetime.datetime.now()
print("\nNeighborhood Size Esitmation Time:")
print((endtime-starttime))  

# the least K that makes SNN a 1-connected graph is estimated at 20
starttime = datetime.datetime.now()
Edges_sub, Den_sub = KNN_Density_Ball_SNN(PCA_Data_sub, K)
#Edges, Den = SNN_Graph(PCA_Data, K)
endtime = datetime.datetime.now()
print("\nSNN network construction Time:")
print((endtime-starttime))  


### Step-6 Persistenct homology on superlevel set filtration
starttime = datetime.datetime.now()
PH_sub, Cluster_sub, Cluster_U_sub = PersistentHomology(PCA_Data_sub, tSNE_sub, Edges_sub, Den_sub, name, pd_name = "raw", Stratified = True, iter_plot = False, cluster_plot = False)   
endtime = datetime.datetime.now()
print("\nPeaks Identification Time:")
print((endtime-starttime))  


### Step-8 Report peaks_fraction of each cluster
starttime = datetime.datetime.now()
PH_Peaks_sub, Cluster_Peaks_sub = Peaks_Cells_Label(PH_sub, Cluster_sub, Cluster_U_sub)
uniq_label_sub = list(np.unique(Cluster_Peaks_sub.Label.values))
endtime = datetime.datetime.now()
print("\nPeaks Identification Time:")
print((endtime-starttime))  

### Step-9 Rename peaks
rename_dict = dict(zip(uniq_label_sub, np.arange(len(uniq_label_sub))))
Cluster_Peaks_sub.loc[:,"Label"] =  [rename_dict[i] for i in Cluster_Peaks_sub.Label.values]


### Step-9 Iterative KNN to call all cells back
Cluster_All_sub = Iterative_KNN_Calling(Cluster_Peaks_sub, PCA_Data_sub, tSNE_sub, name, iter_plot=False)
nClust_sub = len(list(np.unique(Cluster_All_sub.Label.values)))


### Step-10 Output results
starttime = datetime.datetime.now()
# report parameters
g=open(name + "/summary.txt","w")
g.write("The estimated number of PCs is: " + str(Dim) + "\n")
g.write("The estimated neighborhood is: " + str(K)+ "\n")
g.write("The estimated cluster number is: " +  str(nClust)+ "\n")
#g.write("The ARI index is: " +  str(ARI)+ "\n")
g.close()

# output PCA table with clustering labels
#PCA_Data["RealLabel"] = truelabel
PCA_Data_sub["PeakLabel"] = Cluster_Peaks_sub.Label.values
PCA_Data_sub["AllLabel"] = Cluster_All_sub.Label.values
PCA_Data_sub.to_csv(name + "/PCA.csv")

# output normalized gene expression table with clustering labels
S4_sub = S4_sub.T
#S4["RealLabel"] = truelabel
S4_sub["PeakLabel"] = Cluster_Peaks_sub.Label.values
S4_sub["AllLabel"]  = Cluster_All_sub.Label.values

S4_sub.to_csv(name + "/S4.csv")

# output other important results
tSNE_sub.to_csv(name + "/tSNE.csv")

# output clustering results
PH_sub.to_csv(name + "/PersistenceDiagram.csv")
PH_Peaks_sub.to_csv(name + "/PersistenceDiagramPeaks.csv")
Cluster_All_sub.to_csv(name + "/ClusteringLabel.csv")
Cluster_Peaks_sub.to_csv(name + "/ClusteringLabelPeaks.csv")

### Step-11 Further analysis
# intra-cluster stability(pair-wise pearson correlation)
PCA_Data_sub = pd.read_csv(name + "/PCA.csv", header=0, index_col=0, sep=",")
IntraClusterCorrelation(PCA_Data_sub, "AllLabel", name)
IntraClusterCorrelation(PCA_Data_sub, "PeakLabel", name)
peak_tpm = RepresentativeTranscriptome2(S4_sub,"PeakLabel")
all_tpm = RepresentativeTranscriptome2(S4_sub,"AllLabel")

Markers = ClusterSpecificGeneSet(S4_sub, name, top=1)








