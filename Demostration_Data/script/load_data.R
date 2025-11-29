setwd("C:/Users/Administrator/Desktop/pbmc3k/")

Seurat_object <- readRDS("GSE197017_scRNAseq_gene_expression_Seurat.rds")
labels = data.frame(Batch=Seurat_object$batch, Treat=Seurat_object$Treatment, ID=Seurat_object$orig.ident, SPiDer=Seurat_object$SPiDER, CD45=Seurat_object$CD45, Seurat_Cluster=Seurat_object$seurat_clusters)
write.table(labels,"GSE197017_Label.tsv",sep="\t",row.names = F,col.names = T,quote = F)

counts = Seurat_object@assays$RNA@counts
count = data.frame(as.matrix(counts))
write.table(count,"GSE197017_Count.tsv",sep="\t",row.names = T,col.names = T,quote = F)


umap = Seurat_object@reductions$umap@cell.embeddings
write.table(umap,"GSE197017_Umap.tsv",sep="\t",row.names = T,col.names = T,quote = F)