library(Seurat)
library(patchwork)

setwd("C:/Users/Administrator/Desktop/pbmc3k/")
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "pbmc3k_filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

hvgs <- pbmc@assays$RNA@var.features
hvgs <- data.frame(hvgs)
write.table(hvgs,"Seurat_HVGs.csv", row.names = F, col.names = T,sep=",",quote=F)