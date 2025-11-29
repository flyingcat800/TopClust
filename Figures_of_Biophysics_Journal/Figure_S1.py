# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 17:53:21 2022

@author: huchu
"""

### Import Library
import os
import gc
import datetime
import networkx as nx
import numpy as np
import pandas as pd
import scienceplots
from umap import UMAP
from math import pi
from copy import deepcopy
from numpy.random import multivariate_normal
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib import gridspec
from matplotlib.colors import ListedColormap
from mpl_toolkits import axisartist
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import gamma
from scipy.spatial.distance import squareform
from scipy.spatial.distance import pdist
from sklearn.metrics import adjusted_mutual_info_score 
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import NearestNeighbors
from sklearn.neighbors import kneighbors_graph


### Working Dir
path = "D:/Project/tda_paper/5_Cell_Line/Layer_1\\"
plt.style.use(['science','nature','no-latex'])
os.chdir(path)


dirs = os.listdir()
dirs = [r for r in dirs if not r.endswith(".png")]
dirs = [r for r in dirs if not r.endswith(".jpg")]
dirs = [r for r in dirs if not r.endswith(".pdf")]
dirs = [r for r in dirs if not r.endswith(".tiff")]
dirs = [r for r in dirs if not r.endswith(".R")]
dirs = [r for r in dirs if not r.endswith(".tsv")]
dirs = [r for r in dirs if not r.endswith(".csv")]
dirs = [r for r in dirs if not r.endswith(".txt")]

result = []
for p in dirs:
    summary = os.path.join(path,p,"summary.txt") 
    if not summary:
        print(p)
    summary_text = pd.read_csv(summary, header = None , index_col= None, sep = "\t" )
    cell = summary_text.iloc[2,0]
    cellNum = cell.split(":")[1].split(" ")[1]
    result.append( [int(p), int(cellNum)] )

result = pd.DataFrame(result, columns=["K", "Cluster_Number"])
result = result.sort_values(by="K", ascending=True)
result.to_csv(path + "K_ClusterNum.tsv", sep="\t", index=None)



Dirs = [int(i) for i in dirs]
Dirs.sort()
mat = pd.DataFrame()
Dirs.remove(min(Dirs))


result = []
for p in Dirs:
    #print(p)
    tSNE_before = str(p-1) + "/tSNE.csv"
    tSNE_before = pd.read_csv(tSNE_before, header = 0 , index_col= 0, sep = "," )
    tSNE_current = str(p) + "/tSNE.csv"
    tSNE_current = pd.read_csv(tSNE_current, header = 0 , index_col= 0, sep = "," )
    ARI_Neighbor = np.round(adjusted_rand_score(tSNE_before["L1A"], tSNE_current["L1A"]),3) 
    result.append( [int(p), ARI_Neighbor] )
    #g.write(str(p) + "\t" + str(ARI_Neighbor)  + "\n")
result = pd.DataFrame(result, columns=["K", "ARI_Neighbor"])  
result = result.sort_values(by="K", ascending=True) 
result.to_csv(path + "K_ARI2.tsv", sep="\t", index=None)



name = path
dat1 = pd.read_csv(name + "K_ARI2.tsv", header = 0 , index_col= None, sep = "\t" )
dat2 = pd.read_csv(name + "K_ClusterNum.tsv", header = 0 , index_col= None, sep = "\t" )
#dat1 = dat1.sort_values(by=["K"],ascending=True)
data_ari = pd.concat([dat1,dat2],axis=1)
data_ari = data_ari.iloc[:,[0,1,3]]






### Plotting 
### Image 1
fig = plt.figure(constrained_layout=True, figsize=(6,6))
gs = gridspec.GridSpec(2, 2, figure=fig)
ax1 = fig.add_subplot(gs[0, :])
ax3 = fig.add_subplot(gs[1, 0])
ax4 = fig.add_subplot(gs[1, 1])

# subplot_71
# 左侧Y轴
ax1.set_xlabel('K')
ax1.set_ylabel('ClustNum', color='darkred')
ax1.plot(data_ari["K"], data_ari["Cluster_Number"], color='darkred',linestyle='--', marker='o', label='ClustNum' ,markersize=2,linewidth=0.5)
ax1.tick_params(axis='y', labelcolor='darkred')

# 右侧Y轴
ax2 = ax1.twinx()
ax2.set_ylabel('ARI', color='darkcyan')
ax2.plot(data_ari["K"], data_ari["ARI_Neighbor"], color='darkcyan', linestyle='--', marker='o',label='ARI', markersize=2,linewidth=0.5)
ax2.tick_params(axis='y', labelcolor='darkcyan')
ax2.set_ylim([0.8, 1.02])

# 添加图例
lines, labels = ax7.get_legend_handles_labels()
lines2, labels2 = ax8.get_legend_handles_labels()
ax7.legend(lines + lines2, labels + labels2, loc="center right")

# 添加标题
ax7.annotate("g",xy=(-0.046,1.1),  xycoords = 'axes fraction', weight = "bold", fontsize = 14)
ax7.set_title("Clustering stability across different Ks")


# subplot_8 
# peak cell label
# the chosen K = 57
tSNE = pd.read_csv("57/tSNE.csv", header = 0)
tSNE_L1P = tSNE.loc[tSNE.L1P == 1] 
tSNE_L2P = tSNE.loc[tSNE.L1P == 2] 
tSNE_L3P = tSNE.loc[tSNE.L1P == 3] 
tSNE_ambigous = tSNE.loc[(tSNE.L1P == 0)  , :] 
ax9.plot(tSNE_L1P["tSNE-2"], tSNE_L1P["tSNE-1"], linestyle='None', color="crimson", marker='o', markersize=2, linewidth=0.5, label = "1")
ax9.plot(tSNE_L2P["tSNE-2"], tSNE_L2P["tSNE-1"], linestyle='None', color="coral", marker='o', markersize=2, linewidth=0.5, label = "2")
ax9.plot(tSNE_L3P["tSNE-2"], tSNE_L3P["tSNE-1"], linestyle='None', color="dodgerblue", marker='o', markersize=2, linewidth=0.5, label = "3")
ax9.plot(tSNE_ambigous["tSNE-2"], tSNE_ambigous["tSNE-1"], linestyle='None', color="grey", marker='o', markersize=2, linewidth=0.5, label = "3")

ax9.invert_yaxis() 
ax9.set_title("Clustering label of peak cells")
ax9.annotate("h",xy=(-0.10,1.1),  xycoords = 'axes fraction', weight = "bold", fontsize = 14)

# subplot_9 
# peak cell label
# the chosen K = 57
tSNE_L1A = tSNE.loc[tSNE.L1A == 1] 
tSNE_L2A = tSNE.loc[tSNE.L1A == 2] 
tSNE_L3A = tSNE.loc[tSNE.L1A == 3] 
ax10.plot(tSNE_L1A["tSNE-2"], tSNE_L1A["tSNE-1"], linestyle='None', color="crimson", marker='o', markersize=2, linewidth=0.5, label = "1")
ax10.plot(tSNE_L2A["tSNE-2"], tSNE_L2A["tSNE-1"], linestyle='None', color="coral", marker='o', markersize=2, linewidth=0.5, label = "2")
ax10.plot(tSNE_L3A["tSNE-2"], tSNE_L3A["tSNE-1"], linestyle='None',color="dodgerblue", marker='o', markersize=2, linewidth=0.5, label = "3")
ax10.invert_yaxis() 
ax10.set_title("Clustering label of all cells")
ax10.annotate("i",xy=(-0.10,1.1),  xycoords = 'axes fraction', weight = "bold", fontsize = 14)


### Overall settings
fig.tight_layout()#调整整体空白
plt.subplots_adjust(wspace = 0.8, hspace = 0.6)  #调整子图间距
plt.savefig("Figure1.pdf", dpi=1000)
plt.savefig("Figure1.tiff", dpi=1000)
plt.savefig("Figure1.png", dpi=1000)
plt.close("all")






