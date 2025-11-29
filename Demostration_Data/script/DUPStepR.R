#install.packages("DUBStepR")
setwd("C:/Users/Administrator/Desktop/pbmc3k/")
library(DUBStepR)
library(Seurat)
library(dplyr)

pbmc = read.table("pbmc3k_raw.txt", header = T, sep="\t" , check.names = F, stringsAsFactors = F)
seuratObj <- CreateSeuratObject(counts = pbmc, assay = "RNA", project = "PBMC")
seuratObj <- NormalizeData(object = seuratObj, normalization.method = "LogNormalize")

dubstepR.out <- DUBStepR(input.data = seuratObj@assays$RNA@data, min.cells = 0.05*ncol(seuratObj), optimise.features = TRUE, k = 10, num.pcs = 20, error = 0)
seuratObj@assays$RNA@var.features <- dubstepR.out$optimal.feature.genes
DUPStepR_Genes <- dubstepR.out$optimal.feature.genes
DUPStepR_Genes <- data.frame(DUPStepR_Genes)
write.table(DUPStepR_Genes,"DUPStepR_Genes.txt",col.names = T, row.names = F, sep="\t", quote=F)





setwd("C:/Users/Administrator/Desktop/pbmc3k/")
library(DUBStepR)
library(Seurat)
library(dplyr)

pbmc = read.table("GSM3618014_gene_count.csv", header = T, sep="," , check.names = F, stringsAsFactors = F, row.names = 1)
seuratObj <- CreateSeuratObject(counts = pbmc, assay = "RNA", project = "PBMC")
seuratObj <- NormalizeData(object = seuratObj, normalization.method = "LogNormalize")

dubstepR.out <- DUBStepR(input.data = seuratObj@assays$RNA@data, min.cells = 0.05*ncol(seuratObj), optimise.features = TRUE, k = 10, num.pcs = 20, error = 0)
seuratObj@assays$RNA@var.features <- dubstepR.out$optimal.feature.genes
DUPStepR_Genes <- dubstepR.out$optimal.feature.genes
DUPStepR_Genes <- data.frame(DUPStepR_Genes)
write.table(DUPStepR_Genes,"DUPStepR_Genes.txt",col.names = T, row.names = F, sep="\t", quote=F)



