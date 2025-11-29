# TopClust

TopClust is a clustering algorithm for scRNA-seq data clustering.

Single-cell RNA sequencing (scRNA-seq) clustering methods often require subjective parameter selection, leading to biased or inconsistent results. Here, we present TopClust, a data-driven method for scRNA-seq analysis based on topological data analysis (TDA). TopClust constructs a high-dimensional density map via a Shared Nearest Neighbor (SNN) model and uses TDA to identify the local maximum (summit cell) as well as a core group of cells of each cluster (peak cells) in the map. We validate our method on benchmark datasets (lung cancer cell lines and FACS-sorted immune cells), demonstrating excellent agreement with the known cell types. We show that analyses of the peak cells alone reduces transcriptional variation that characterizes each cluster to improve accuracy in cell type annotation. We calculate a directional vector (DV) between pairs of summit cells and then project the genes on this vector to identify the most cluster-discriminating genes, supplementing conventional differential expression analysis. Overall, TopClust enables unbiased, reproducible scRNA-seq clustering, offering a powerful tool for examining cellular heterogeneity.

The Demostration_Data filefolder has the scripts to run TopClust with a benchmark dataset.
