### Step-0 Import packages
# pip install SciencePlots
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
from matplotlib import gridspec
from matplotlib.colors import ListedColormap
from mpl_toolkits import axisartist
from mpl_toolkits.mplot3d import Axes3D
from sklearn.neighbors import kneighbors_graph
from scipy.special import gamma
from scipy.spatial.distance import squareform
from scipy.spatial.distance import pdist
from sklearn.metrics import adjusted_rand_score

### Working Dir
plt.style.use(['science','nature','no-latex'])
plt.rcParams['font.sans-serif'] = "Helvetica"
plt.rcParams['font.family'] = "sans-serif"


path = "E:/Project/tda_paper/3_Cell_Line/Layer_1/"
os.chdir(path)



### Functions
def KNN_Density_Ball_SNN(PCAData,K):
    starttime = datetime.datetime.now()
    N = PCAData.shape[0]
    d = PCAData.shape[1]
    KNN_Dist = kneighbors_graph(PCAData, n_neighbors=K, mode='distance', metric='minkowski', p=2, include_self=False)
    Distance = np.round(squareform(pdist(np.array(PCAData), metric='euclidean')),4)
    R = np.array([])
    for i in range(Distance.shape[0]):
        line = Distance[i]
        line.sort()
        R = np.append(R,line[K])
        
    v = (pi**(d/2))/gamma(d/2+1)*(R**d)
    del Distance,R
    Den = (K/float(N))/v
    Den = Den/Den.sum()
    
    Den = pd.Series(Den)
    KNNetwork = kneighbors_graph(PCAData, n_neighbors=K, include_self=False)
    KNNeighbor =[set(KNNetwork[i].nonzero()[1]) for i in range(len(PCAData))]
    starttime = datetime.datetime.now()
    SNN= np.array([[len(KNNeighbor[i].intersection(KNNeighbor[j]))/K for j in range(len(KNNeighbor))] for i in range(len(KNNeighbor))])
    SNN = pd.DataFrame(SNN)    
    SNN = np.triu(SNN,1)
    SNN = pd.DataFrame(SNN)
    Edges = SNN.stack()
    Edges.index.rename(['NodeA', 'NodeB'], inplace=True)
    Edges = Edges.to_frame('Weight').reset_index()
    Edges = Edges.loc[Edges.NodeA != Edges.NodeB,:]   
    prune = 6/(K+K-6)
    #prune = 0.6
    Edges = Edges.loc[Edges.Weight >= prune,:]
    Edges["NodeADen"] = Den[Edges.NodeA].values
    Edges["NodeBDen"] = Den[Edges.NodeB].values
    ## To define the higher vertice as the "From" vertice, the lower vertice as the "To" vertice 
    Diff = (Edges["NodeADen"] - Edges["NodeBDen"]) < 0
    Edges.loc[Diff,['NodeA','NodeB']] = Edges.loc[Diff,['NodeB','NodeA']].values
    Edges.loc[Diff,['NodeADen','NodeBDen']] = Edges.loc[Diff,['NodeBDen','NodeADen']].values
    # Edges = Edges.sort_values(by="NodeADen" , ascending = False)
    # Edges.index = np.arange(Edges.shape[0])
    
    # lowest = Edges.NodeBDen.argmin()
    # lowestNodeB = int(Edges.values[lowest,1])
    # lowestNodeBDen = Edges.values[lowest,4]
    # nodes = np.unique(Edges.NodeA)
    # np.random.shuffle(nodes)
    # nodes = nodes[0:50]
    # append = pd.DataFrame(columns=Edges.columns)
    # append.NodeA = nodes
    # append.NodeB = lowestNodeB
    # append.Weight = prune
    # append.NodeADen = Den[append.NodeA].values
    # append.NodeBDen = lowestNodeBDen
    # Edges = pd.concat([Edges,append],axis=0,ignore_index=True)
    Edges["EdgeDen"] = (Edges["NodeADen"] + Edges["NodeBDen"])/2
    Edges = Edges.sort_values(by=["EdgeDen","NodeADen"] , ascending = False)
    Edges.index = np.arange(Edges.shape[0])    
    gc.collect()
    endtime = datetime.datetime.now()
    print("\nShared Nearest Networt Construction Time:")
    print((endtime-starttime))  
    return(Edges,Den)

def root(i,idx):
    while i != idx[i]:
        idx[i] = idx[ID[i]]
        i = idx[i]
    return i

def union(f,t,idx):
    F = root(f,idx)
    T = root(t,idx)
    idx[T] = F




### Data 
# Generate synthetic data with 3 clusters
mean_1 = [0, 4]
mean_2 = [-3, 0]
mean_3 = [3, 0]
cov = [[1, 0], [0, 1]]  # diagonal covariance

# Each cluster follows a 2-d Gaussian distribution
np.random.seed(0)
data_1 = np.random.multivariate_normal(mean_1, cov, 500)
data_2 = np.random.multivariate_normal(mean_2, cov, 300)
data_3 = np.random.multivariate_normal(mean_3, cov, 200)

# Cluster boundry
print(data_1.min(axis=0),data_1.max(axis=0))
print(data_2.min(axis=0),data_2.max(axis=0))
print(data_3.min(axis=0),data_3.max(axis=0))

# Merge 3 clusters
data = np.concatenate((data_1, data_2, data_3), axis=0)
plt.plot(data[:,0], data[:,1], 'p')
plt.axis('equal')
plt.show()

# Build Node-weighted SNN network
ID = np.array(list(range(data.shape[0])))
SNN, Den = KNN_Density_Ball_SNN(data, K = 30)
edge_list = [(SNN.iloc[i,0],SNN.iloc[i,1]) for i in range(SNN.shape[0])]
G = nx.from_edgelist(edge_list)

# Generate coordiante of nodes and edges
# 3D coordinates of nodes
node_xyz = np.concatenate([data, Den.values.reshape(-1,1)],axis=1)
# 3D coordinates of edges
edge_xyz = np.array([(node_xyz[u], node_xyz[v]) for u, v in edge_list])
# 2D coordinates of nodes
node_xy = deepcopy(data)
# 2D coordinates of edges
edge_xy = np.array([(node_xy[u], node_xy[v]) for u, v in edge_list])




dirs = os.listdir()
dirs = [r for r in dirs if not r.endswith(".png")]
dirs = [r for r in dirs if not r.endswith(".pdf")]
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

g = open(path + "K_ARI2.tsv","a")
g.write("K\tARI_Neighbor\n")

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
fig = plt.figure(constrained_layout=True, figsize=(7,7))
gs = gridspec.GridSpec(3, 3, figure=fig)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[0, 2])
ax4 = fig.add_subplot(gs[1, 0], projection="3d")
ax5 = fig.add_subplot(gs[1, 1], projection="3d")
ax6 = fig.add_subplot(gs[1, 2])
ax7 = fig.add_subplot(gs[2, :])



# subplot_1
m1 = np.ones((200,150))
#ax1.imshow(m1, cmap=ListedColormap("deepskyblue"), alpha=0.7)
ax1.imshow(m1, cmap=ListedColormap("lightseagreen"), alpha=0.7)
grid_size=10
ax1.grid(which='both', color='black', linestyle='--', linewidth=0.3)
ax1.set_xticks(np.arange(0, m1.shape[1], grid_size) - 0.5, minor=True)
ax1.set_yticks(np.arange(0, m1.shape[0], grid_size) - 0.5, minor=True)
ax1.grid(which='minor', color='black', linestyle='--', linewidth=0.3)
ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_xlabel("c Cells", labelpad=10)
ax1.set_ylabel("r Genes", labelpad=10)
ax1.set_title("Raw data", pad=10)
ax1.annotate("a",xy=(-0.54,1.1),  xycoords = 'axes fraction', weight = "bold", fontsize = 14)


# subplot_2
m2 = np.ones((160,120))
ax2.imshow(m2, cmap=ListedColormap("coral"), alpha=0.7)
grid_size=10
ax2.grid(which='both', color='black', linestyle='--', linewidth=0.3)
ax2.set_xticks(np.arange(0, m1.shape[1], grid_size) - 0.5, minor=True)
ax2.set_yticks(np.arange(0, m1.shape[0], grid_size) - 0.5, minor=True)
ax2.grid(which='minor', color='black', linestyle='--', linewidth=0.3)
ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_xlabel("n Cells", labelpad=10, color='coral')
ax2.set_ylabel("m Genes", labelpad=10, color='coral')
ax2.xaxis.label.set_x(0.45) 
ax2.yaxis.label.set_y(0.6) 
ax2.set_title("Quality Control", pad=10)
ax2.annotate("b",xy=(-0.54,1.1),  xycoords = 'axes fraction', weight = "bold", fontsize = 14)


# subplot_3
m3 = np.ones((60,120))
ax3.imshow(m3, cmap=ListedColormap("crimson"), alpha=0.7)
grid_size=10
ax3.grid(which='both', color='black', linestyle='--', linewidth=0.3)
ax3.set_xticks(np.arange(0, m1.shape[1], grid_size) - 0.5, minor=True)
ax3.set_yticks(np.arange(0, m1.shape[0], grid_size) - 0.5, minor=True)
ax3.grid(which='minor', color='black', linestyle='--', linewidth=0.3)
ax3.set_xticks([])
ax3.set_yticks([])
ax3.set_xlabel("n Cells",                labelpad=10, color='crimson')
ax3.set_ylabel("d PCs", labelpad=10, color='crimson')
ax3.set_title("PCA Dimension Reduction", pad=10)
ax3.xaxis.label.set_x(0.45) 
ax3.yaxis.label.set_y(0.83) 
ax3.annotate("c",xy=(-0.54,1.1),  xycoords = 'axes fraction', weight = "bold", fontsize = 14)


# subplot_4
#ax4 = fig.add_subplot(2,3,4, projection="3d")
N = 1
idx = np.array(list(range(data.shape[0])))
edge_highlight = edge_xyz[0:N]
raw_shape = edge_highlight.shape
npoints = raw_shape[0] * raw_shape[1]
node_highlight = edge_highlight.reshape(1,npoints,3)[0]
# Plot the nodes - alpha is scaled by "depth" automatically
ax4.plot(*node_xyz.T, "o", markersize = 3, mec = "black", mew = 0.1, mfc="lightgrey")
# Plot water level 
xticks =  np.linspace(data.min(axis=0)[0], data.max(axis=0)[0], 10)
yticks =  np.linspace(data.min(axis=0)[1], data.max(axis=0)[1], 10)
xx, yy = np.meshgrid(xticks, yticks)
den = node_highlight.min(axis=0)[2]
zz = np.ones_like(xx) * den

# Formatting axis
ax4.zaxis.set_rotate_label(False) 
ax4.set_xlabel("Dimension-1",labelpad=-8)
ax4.set_ylabel("Dimension-2",labelpad=-8)
ax4.set_zlabel("Density",    labelpad=-6, rotation=90)
ax4.tick_params(axis="x",pad=-4)
ax4.tick_params(axis="y",pad=-4)
ax4.tick_params(axis="z",pad= -2)
ax4.view_init(30,30)
ax4.set_title("Node-weighted SNN Graph")
ax4.annotate("d",xy=(-0.32,1.13),  xycoords = 'axes fraction', weight = "bold", fontsize = 14)


# subplot_5
## the 3D network
# the first n edges of highest density
N = 17000
idx = np.array(list(range(data.shape[0])))
#ax5 = fig.add_subplot(2,3,5, projection="3d")
edge_highlight = edge_xyz[0:N]
raw_shape = edge_highlight.shape
npoints = raw_shape[0] * raw_shape[1]
node_highlight = edge_highlight.reshape(1,npoints,3)[0]
# Union Find
edges = edge_list[0:N]
for i in range(N):
    union(edges[i][0],edges[i][1],idx)
idx = [root(n,idx) for n in ID ]
seq, count = np.unique(idx, return_counts=True)
codes = (np.where(count>1))[0]
if len(codes) == 1:
    code1 = seq[codes[0]]
    code2 = ""
    code3 = ""
if len(codes) == 2:
    code1 = seq[codes[0]]
    code2 = seq[codes[1]]
    code3 = ""
if len(codes) == 3:
    code1 = seq[codes[0]]
    code2 = seq[codes[1]]
    code3 = seq[codes[2]]
# Plot the nodes - alpha is scaled by "depth" automatically
ax5.plot(*node_xyz.T, "o", markersize = 3, mec = "black", mew = 0.1, mfc="lightgrey")
#ax5.plot(node_highlight[:,0], node_highlight[:,1], node_highlight[:,2], "o", markersize = 6, mec = "black", mew = 0.3, mfc="lightgrey")
# Labeling highlight nodes
n1 = ID[idx == code1].tolist()
n2 = ID[idx == code2].tolist()
n3 = ID[idx == code3].tolist()
ax5.plot(node_xyz[n1,0], node_xyz[n1,1], node_xyz[n1,2], "o", markersize = 3, mec = "black", mew = 0.1, mfc="crimson")
ax5.plot(node_xyz[n2,0], node_xyz[n2,1], node_xyz[n2,2], "o", markersize = 3, mec = "black", mew = 0.1, mfc="coral")
ax5.plot(node_xyz[n3,0], node_xyz[n3,1], node_xyz[n3,2], "o", markersize = 3, mec = "black", mew = 0.1, mfc="dodgerblue")  
# Plot the edges
# for edge in edge_highlight:
#     ax5.plot(*edge.T, linestyle="dashed", color="lightgrey", linewidth=0.2)
# Plot water level 
xticks =  np.linspace(data.min(axis=0)[0], data.max(axis=0)[0], 10)
yticks =  np.linspace(data.min(axis=0)[1], data.max(axis=0)[1], 10)
xx, yy = np.meshgrid(xticks, yticks)
den = node_highlight.min(axis=0)[2]
zz = np.ones_like(xx) * den
ax5.plot_surface(xx, yy, zz, alpha=0.2, color = "springgreen")
# Formatting axis
ax5.zaxis.set_rotate_label(False) 
ax5.set_xlabel("Dimension-1",labelpad=-8)
ax5.set_ylabel("Dimension-2",labelpad=-8)
ax5.set_zlabel("Density",    labelpad=-6, rotation=90)
ax5.tick_params(axis="x",pad=-4)
ax5.tick_params(axis="y",pad=-4)
ax5.tick_params(axis="z",pad=-2)
ax5.set_title("Superlevel Set Filtration")
ax5.view_init(30,30)
ax5.annotate("e",xy=(-0.3,1.15),  xycoords = 'axes fraction', weight = "bold", fontsize = 14)


# subplot_6
## Persistence Barcode
#ax6 = fig.add_subplot(2,3,6)
d2max = Den[n2].max()
d2min = Den[n2].min()
ax6.barh(["1","2","3"], [Den[n1].min(), Den[n2].min(), Den[n3].min()], 0.2, color="w")
ax6.barh(["1","2","3"], [Den[n1].max()-Den[n1].min(), Den[n2].max()-Den[n2].min(), Den[n3].max()-Den[n3].min()], 0.2, left=[Den[n1].min(), Den[n2].min(), Den[n3].min()], color=["crimson","coral","dodgerblue"])
ax6.invert_yaxis() 
ax6.invert_xaxis() 
#ax6.ticklabel_format(style='scientific', scilimits=(0,0), useMathText=True, axis="x")
ax6.set_xlabel("Density Level")
ax6.set_ylabel("Persistence Barcodes ($H_{0}$)")
ax6.set_xlim(Den.max(),Den.min())
ax6.set_title("Persistence Barcodes")
ax6.annotate("f",xy=(-0.3,1.1),  xycoords = 'axes fraction', weight = "bold", fontsize = 14)

# subplot_7
# 左侧Y轴
ax7.set_xlabel('K')
ax7.set_ylabel('ClustNum', color='darkred')
ax7.plot(data_ari["K"], data_ari["Cluster_Number"], color='darkred',linestyle='--', marker='o', label='ClustNum' ,markersize=2,linewidth=0.5)
ax7.tick_params(axis='y', labelcolor='darkred')
ax7.set_yticks([0,5,10,15])

# 右侧Y轴
ax8 = ax7.twinx()
ax8.set_ylabel('ARI', color='darkcyan')
ax8.plot(data_ari["K"], data_ari["ARI_Neighbor"], color='darkcyan', linestyle='--', marker='o',label='ARI', markersize=2,linewidth=0.5)
ax8.tick_params(axis='y', labelcolor='darkcyan')
ax8.set_ylim([0.8, 1.02])

# 添加图例
lines, labels = ax7.get_legend_handles_labels()
lines2, labels2 = ax8.get_legend_handles_labels()
ax7.legend(lines + lines2, labels + labels2, loc="center right")

# 添加标题
ax7.annotate("g",xy=(-0.075,1.1),  xycoords = 'axes fraction', weight = "bold", fontsize = 14)



### Overall settings
fig.tight_layout()#调整整体空白
plt.subplots_adjust(wspace = 0.6, hspace = 0.6)  #调整子图间距
plt.savefig("Figure1.pdf", dpi=1000)
plt.savefig("Figure1.tiff", dpi=1000)
plt.savefig("Figure1.png", dpi=1000)
plt.close("all")




