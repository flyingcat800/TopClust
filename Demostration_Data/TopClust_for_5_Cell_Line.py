#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 17:24:28 2022

@author: cshu
"""

### Step-0 Import packages
import os
os.chdir("E:/Project/tda_clust/")
from script.functions import *


### Step-1 Specify dataset
if not os.path.exists("Summary"):
    os.makedirs("Summary")

starttime = datetime.datetime.now()
file =  "Input/GSM3618014_gene_count.csv"
label = "Input/GSM3618014_gene_count_metadata.csv"
name = "5_Cell_Line"
Dir = "Layer-1/" + name
if not os.path.exists(Dir):
    os.makedirs(Dir)
    
### Step-2 Data processing
Truelabel = pd.read_csv(label, header=0, index_col=0, sep=",")
truelabel = Truelabel["cell_line_demuxlet"]
kept_cells = Truelabel.index.tolist()

S1 = pd.read_csv(file, header=0, index_col=0, sep=",")
S1 = S1.loc[:, kept_cells]
S2 = Quality_Control(S1, cellNum=0, UMINumLowLimit= 500,UMINumUpLimit=3000,mtPct=0.1)
S3 = Normalize_Rawdata(S2)
S4_1 = Select_HEGenes(S3,ngenes=2000)
# dubstepR = pd.read_csv("GSM3618014_DUPStepR_Genes.txt", header=0, index_col= None, sep="\t")
# dubstepR = dubstepR.DUPStepR_Genes.values
# S4_1 = S3.loc[dubstepR,:]
S5 = Scaling_Data(S4_1)
S6 = S5.T   



### Step-3 PCA dimension reduction for density peaks finding
Dimension = 100
pca = PCA(n_components=Dimension,  random_state = 1) 
PCA_s = pd.DataFrame(pca.fit_transform(S6), index=S6.index, columns=["PC%d" % k for k in range(1,Dimension+1)]).iloc[:,:Dimension]
ElbowPlot(S6,Dir)
Dim = PCA_Permutation(S6, Dir, nperm = 100, nPC = 40)

PCA_Data_1 = np.round(PCA_s.iloc[:,0:Dim],4)
PCA_Data_Prime = deepcopy(PCA_Data_1)
CellID = pd.DataFrame(PCA_Data_1.index, columns=["ID"])
nCells = PCA_Data_1.shape[0]

for K in range(10,101):
    print(K)
    PCA_Data = deepcopy(PCA_Data_1)
    S4 = deepcopy(S4_1)
    Dir = "Layer-1/" + name + "/" + str(K) 
    if not os.path.exists(Dir):
        os.makedirs(Dir)
    pca_components = pd.DataFrame(pca.components_, index=["PC%d" % k for k in range(1,Dimension+1)], columns = S6.columns )
    pca_components.to_csv(Dir+"/Layer1_HVGs_Contribution_to_Each_PCs.csv")
    
    ### Step-4 tSNE dimension reduction for visualization
    umap_coor = UMAP(random_state=1).fit(PCA_Data)
    umap_coor = pd.DataFrame(umap_coor.embedding_,columns=["UMAP-1","UMAP-2"])
     
    tsne = TSNE(n_components=2, random_state=1).fit_transform(PCA_Data)
    tSNE = pd.DataFrame(tsne,columns=["tSNE-1","tSNE-2"]) 
    # tSNE = pd.read_csv("tSNE_Ref.csv", header=0, index_col=0, sep=",")
    tSNE["TrueLabel"] = truelabel.values
    Scatter_Plot(tSNE,Dir+"/tSNE","tSNE-1","tSNE-2")
    
    for i in range(Dim):
        for j in range(i+1,Dim):
            a = i + 1
            b = j + 1
            Scatter_Plot_Label(PCA_Data.loc[:,["PC"+str(a),"PC"+str(b)]], Dir + "/Label_of_reference_pca_"+"_PC"+str(a)+"_PC"+str(b), truelabel, "PC-"+str(a), "PC-"+str(b), legend={"loc":'center right', "bbox_to_anchor":(1.4,0.5), "fontsize":5, "ncol":1})
    
    Scatter_Plot_Label(umap_coor, Dir + "/Label_of_reference_UMAP", truelabel, "UMAP-1", "UMAP-2",legend={"loc":'center right', "bbox_to_anchor":(1.4,0.5), "fontsize":5, "ncol":1})
    
    endtime = datetime.datetime.now()
    print("\nPreprocessing Time:")
    print((endtime-starttime))  
    tSNEPlot(tSNE,Dir+"/Label_of_Reference",truelabel,legend={"loc":'center right', "bbox_to_anchor":(1.3,0.5),"fontsize":5, "ncol":1})
    ### Step-5 Node-weighted SNN Graph
    starttime = datetime.datetime.now()
    #K = K_Estimation(PCA_Data,5)
    
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
    PH, Cluster, Cluster_U = PersistentHomology(PCA_Data, tSNE, Edges, Den, Dir, pd_name = "raw", Stratified = True, iter_plot = False, cluster_plot = False)   
    endtime = datetime.datetime.now()
    print("\nPeaks Identification Time:")
    print((endtime-starttime))  
    
    
    ### Step-7 Report peaks_fraction of each cluster
    starttime = datetime.datetime.now()
    PH_Peaks, Cluster_Peaks = Peaks_Cells_Label(PH, Cluster, Cluster_U)
    uniq_label = list(np.unique(Cluster_Peaks.Label.values))
    if 0 not in uniq_label:
        uniq_label.append(0)
        uniq_label.sort()
    endtime = datetime.datetime.now()
    print("\nPeaks Identification Time:")
    print((endtime-starttime))  
    
    ### Step-8 Rename peaks
    rename_dict = dict(zip(uniq_label, np.arange(len(uniq_label))))
    Cluster_Peaks.loc[:,"Label"] =  [rename_dict[i] for i in Cluster_Peaks.Label.values]
    
    
    ### Step-9 Iterative KNN to call all cells back
    Cluster_All = Iterative_KNN_Calling(Cluster_Peaks, PCA_Data, tSNE, Dir, iter_plot=False)
    nClust = len(list(np.unique(Cluster_All.Label.values)))
    
    predictlabel = Cluster_All.Label.values
    #ARI = np.round(adjusted_rand_score(truelabel,predictlabel),3)
    endtime = datetime.datetime.now()
    print("\nAssign All Other Cells Time:")
    print((endtime-starttime))
    
    ### Step-10 Output results
    starttime = datetime.datetime.now()
    # report parameters
    g=open(Dir+"/summary.txt","w")
    g.write("The estimated number of PCs is: " + str(Dim) + "\n")
    g.write("The estimated neighborhood is: " + str(K)+ "\n")
    g.write("The estimated cluster number is: " +  str(nClust)+ "\n")
    #g.write("The ARI index is: " +  str(ARI)+ "\n")
    g.close()
    
    
    # output density
    Den = pd.DataFrame(Den)
    Den.index = S6.index
    Den.columns = ["Density"]
    Den.to_csv(Dir + "/Density.csv",index=True)  
    
    
    # output PCA table with clustering labels
    #PCA_Data["RealLabel"] = truelabel
    PCA_Data["L1P"] = Cluster_Peaks.Label.values
    PCA_Data["L1A"] = Cluster_All.Label.values
    #PCA_Data["TrueLabel"] = truelabel.values
    PCA_Data.to_csv(Dir+"/PCA.csv")
    
    # output normalized gene expression table with clustering labels
    S4 = S4.T
    #S4["RealLabel"] = truelabel
    S4["L1P"] = Cluster_Peaks.Label.values
    S4["L1A"]  = Cluster_All.Label.values
    #S4["TrueLabel"] = truelabel.values
    S4.to_csv(Dir+"/S4.csv")
    
    # output other important results
    tSNE["L1P"] = Cluster_Peaks.Label.values
    tSNE["L1A"]  = Cluster_All.Label.values
    tSNE.to_csv(Dir+"/tSNE.csv")
    
    # output clustering results
    PH.to_csv(Dir+"/PersistenceDiagram.csv")
    PH_Peaks.to_csv(Dir+"/PersistenceDiagramPeaks.csv")
    Cluster_All.to_csv(Dir+"/ClusteringLabel.csv")
    Cluster_Peaks.to_csv(Dir+"/ClusteringLabelPeaks.csv")
    
    # ### Step-11 Further analysis
    # # intra-cluster stability(pair-wise pearson correlation)
    # #PCA_Data = pd.read_csv(Dir+"/PCA.csv", header=0, index_col=0, sep=",")
    # IntraClusterCorrelation(PCA_Data, "L1A",Dir)
    # IntraClusterCorrelation(PCA_Data, "L1P", Dir)
    # IntraClusterCorrelation(PCA_Data, "TrueLabel", Dir)
    # PCA_Data["MergedLabel"] = PCA_Data["TrueLabel"]
    # PCA_Data.loc[PCA_Data.MergedLabel == "CD14+ Mono","MergedLabel"] = "Merged_7"
    # PCA_Data.loc[PCA_Data.MergedLabel == "DC","MergedLabel"] = "Merged_7"
    # IntraClusterCorrelation(PCA_Data, "MergedLabel", Dir)
    # Markers = ClusterSpecificGeneSet(S4, Dir, "L1P", top=10)
    # PH.to_csv(Dir+"/PersistenceDiagram.csv")
    # PH_Peaks.to_csv(Dir+"/PersistenceDiagramPeaks.csv")
    # Cluster_All.to_csv(Dir+"/ClusteringLabel.csv")
    # PredictLabel = Cluster_All.Label.values
    # LabelID = np.unique(PredictLabel)
    
    
    
    # find pair-wise markers
    Labels = np.unique(S4.L1A.values)
    # for i in range(len(Labels)):
    #     for j in range(i+1,len(Labels)):
    #         print(Labels[i],Labels[j])
    #         try:
    #             Markers = ClusterSpecificGeneSet(S4, Dir, "L1P", top=10, group=(Labels[i],Labels[j]))
    #         except ValueError:
    #             print('No DEGs')
    
    
    ## Peak cells
    S4["L1D"] = Den.values
    S4["L1Max"] = 0
    
    for i in Labels:
        S4_Peak = S4.loc[S4.L1P==i,:]
        ind = S4_Peak.L1D.argmax()
        index = S4_Peak.index[ind]
        S4_Peak.loc[index,"L1Max"] = 1
        S4.loc[index,"L1Max"] = 1
    
    PCA_Data["L1Max"] = S4.L1Max
    Wilcox = pd.DataFrame(columns=["ClusterA","ClusterB","Wilcox.Pvalue"])
    for i in Labels:
        for j in Labels:
            if j > i:
                S4_Peak_i = S4.loc[S4.L1P==i,:]
                S4_Peak_j = S4.loc[S4.L1P==j,:]
                S4_Peak_i_max = S4_Peak_i.loc[(S4_Peak_i.L1Max == 1),:]
                S4_Peak_j_max = S4_Peak_j.loc[(S4_Peak_j.L1Max == 1),:]
                Genes = S4_Peak_i_max.columns.to_list()[0:2000]
                S4_Peak_i_max = S4_Peak_i_max.iloc[:,0:2000].values
                S4_Peak_j_max = S4_Peak_j_max.iloc[:,0:2000].values
                
                PCA_Data_Peak_i = PCA_Data.loc[PCA_Data.L1P==i,:]
                PCA_Data_Peak_j = PCA_Data.loc[PCA_Data.L1P==j,:]
                PCA_Data_Peak_i_max = PCA_Data_Peak_i.loc[(PCA_Data_Peak_i.L1Max == 1),:]
                PCA_Data_Peak_j_max = PCA_Data_Peak_j.loc[(PCA_Data_Peak_j.L1Max == 1),:]
                pc_ids = PCA_Data_Peak_i_max.columns.to_list()[0:Dim]
                PCA_Data_Peak_i_max = PCA_Data_Peak_i_max.iloc[:,0:Dim].values
                PCA_Data_Peak_j_max = PCA_Data_Peak_j_max.iloc[:,0:Dim].values
                
                Direct_Vec_HVG = S4_Peak_j_max - S4_Peak_i_max 
                Direct_Vec_HVG = Direct_Vec_HVG/np.linalg.norm(Direct_Vec_HVG)
                Direct_Vec_HVG = pd.DataFrame(Direct_Vec_HVG.T,index=Genes,columns=["Direction_Vector"])
                Direct_Vec_HVG.to_csv(Dir+"/DirectionVector_HVG_form_Cluster_" + str(i) + "_to_" + str(j) + ".csv")
                
                Direct_Vec_PCA = PCA_Data_Peak_j_max - PCA_Data_Peak_i_max 
                Direct_Vec_PCA = Direct_Vec_PCA/np.linalg.norm(Direct_Vec_PCA)
                Direct_Vec_PCA_T = pd.DataFrame(Direct_Vec_PCA.T,index=pc_ids,columns=["Direction_Vector"])
                Direct_Vec_PCA_T.to_csv(Dir+"/DirectionVector_PCA_form_Cluster_" + str(i) + "_to_" + str(j) + ".csv")
                
                Peak_Label_i = S4_Peak_i.index.to_list()
                Peak_Label_j = S4_Peak_j.index.to_list()
                PCA_Peak_i = PCA_Data_Prime.loc[Peak_Label_i,:]
                PCA_Peak_j = PCA_Data_Prime.loc[Peak_Label_j,:]
                pcs = PCA_Peak_i.columns.to_list()
                
                
                proj_i = pd.DataFrame((PCA_Peak_i * Direct_Vec_PCA).sum(axis=1),columns=["Projection"])
                proj_j = pd.DataFrame((PCA_Peak_j * Direct_Vec_PCA).sum(axis=1),columns=["Projection"])
                proj_i.to_csv(Dir+"/Cluster" +str(i) +"_Projected_on_DirectionVector_PCA_form_Cluster_" + str(i) + "_to_" + str(j) + ".csv")
                proj_j.to_csv(Dir+"/Cluster" +str(j) +"_Projected_on_DirectionVector_PCA_form_Cluster_" + str(i) + "_to_" + str(j) + ".csv")
                
                Wilcox_stat, Wilcox_pvalue = stats.ranksums(proj_i.Projection, proj_j.Projection, alternative='two-sided')
                Wilcox_line = [i, j, Wilcox_pvalue]
                Wilcox_line = pd.DataFrame(Wilcox_line).T
                Wilcox_line.columns = Wilcox.columns
                Wilcox = pd.concat([Wilcox,Wilcox_line],axis=0)
        PCA_Peak_i.to_csv(Dir+"/PCA_of_Cluster-" + str(i) +"_peak_cells.csv")
    Wilcox.to_csv(Dir+"/Layer-1_Direction_between_summit_Wilcox_rank_sum_test_on_peak_cells.csv")            
    
    Wilcox = pd.DataFrame(columns=["ClusterA","ClusterB","Wilcox.Pvalue"])
    for i in Labels:
        for j in Labels:
            if j > i:
                S4_Peak_i = S4.loc[S4.L1A==i,:]
                S4_Peak_j = S4.loc[S4.L1A==j,:]
                S4_Peak_i_max = S4_Peak_i.loc[(S4_Peak_i.L1Max == 1),:]
                S4_Peak_j_max = S4_Peak_j.loc[(S4_Peak_j.L1Max == 1),:]
                Genes = S4_Peak_i_max.columns.to_list()[0:2000]
                S4_Peak_i_max = S4_Peak_i_max.iloc[:,0:2000].values
                S4_Peak_j_max = S4_Peak_j_max.iloc[:,0:2000].values
                
                PCA_Data_Peak_i = PCA_Data.loc[PCA_Data.L1P==i,:]
                PCA_Data_Peak_j = PCA_Data.loc[PCA_Data.L1P==j,:]
                PCA_Data_Peak_i_max = PCA_Data_Peak_i.loc[(PCA_Data_Peak_i.L1Max == 1),:]
                PCA_Data_Peak_j_max = PCA_Data_Peak_j.loc[(PCA_Data_Peak_j.L1Max == 1),:]
                pc_ids = PCA_Data_Peak_i_max.columns.to_list()[0:Dim]
                PCA_Data_Peak_i_max = PCA_Data_Peak_i_max.iloc[:,0:Dim].values
                PCA_Data_Peak_j_max = PCA_Data_Peak_j_max.iloc[:,0:Dim].values
                
                Direct_Vec_HVG = S4_Peak_j_max - S4_Peak_i_max 
                Direct_Vec_HVG = Direct_Vec_HVG/np.linalg.norm(Direct_Vec_HVG)
                Direct_Vec_HVG = pd.DataFrame(Direct_Vec_HVG.T,index=Genes,columns=["Direction_Vector"])
                Direct_Vec_HVG.to_csv(Dir+"/DirectionVector_HVG_form_Cluster_" + str(i) + "_to_" + str(j) + ".csv")
                
                Direct_Vec_PCA = PCA_Data_Peak_j_max - PCA_Data_Peak_i_max 
                Direct_Vec_PCA = Direct_Vec_PCA/np.linalg.norm(Direct_Vec_PCA)
                Direct_Vec_PCA_T = pd.DataFrame(Direct_Vec_PCA.T,index=pc_ids,columns=["Direction_Vector"])
                Direct_Vec_PCA_T.to_csv(Dir+"/DirectionVector_PCA_form_Cluster_" + str(i) + "_to_" + str(j) + ".csv")
                
                Peak_Label_i = S4_Peak_i.index.to_list()
                Peak_Label_j = S4_Peak_j.index.to_list()
                PCA_Peak_i = PCA_Data_Prime.loc[Peak_Label_i,:]
                PCA_Peak_j = PCA_Data_Prime.loc[Peak_Label_j,:]
                pcs = PCA_Peak_i.columns.to_list()
                
                
                proj_i = pd.DataFrame((PCA_Peak_i * Direct_Vec_PCA).sum(axis=1),columns=["Projection"])
                proj_j = pd.DataFrame((PCA_Peak_j * Direct_Vec_PCA).sum(axis=1),columns=["Projection"])
                proj_i.to_csv(Dir+"/Cluster" +str(i) +"_Projected_on_DirectionVector_PCA_form_Cluster_" + str(i) + "_to_" + str(j) + ".csv")
                proj_j.to_csv(Dir+"/Cluster" +str(j) +"_Projected_on_DirectionVector_PCA_form_Cluster_" + str(i) + "_to_" + str(j) + ".csv")
                
                Wilcox_stat, Wilcox_pvalue = stats.ranksums(proj_i.Projection, proj_j.Projection, alternative='two-sided')
                Wilcox_line = [i, j, Wilcox_pvalue]
                Wilcox_line = pd.DataFrame(Wilcox_line).T
                Wilcox_line.columns = Wilcox.columns
                Wilcox = pd.concat([Wilcox,Wilcox_line],axis=0)
        PCA_Peak_i.to_csv(Dir+"/PCA_of_Cluster-" + str(i) +"_peak_cells.csv")
    Wilcox.to_csv(Dir+"/Layer-1_Direction_between_summit_Wilcox_rank_sum_test_on_all_cells.csv")            
    
    
    
    
                                                               