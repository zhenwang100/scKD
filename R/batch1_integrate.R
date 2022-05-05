# This script is used to integrate samples collected in Batch 1.
# Seurat 3.0.2
library(Seurat)

# Samples and UMI cutoffs
samples<-c("S1_1224", "S1_1229", "S2_0317", "S2_0320", "S3_0414", "S3_0416", "S4_0601", "S4_0603", "H1_0601", "H2_0601", "H3_0601")
low_cutoff<-c(1000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000)
high_cutoff<-c(60000, 60000, 60000, 60000, 60000, 60000, 60000, 60000, 60000, 60000, 60000)

names(low_cutoff)<-samples
names(high_cutoff)<-samples

# Read data and data QC
data_list<-lapply(samples, FUN = function(sample) {
	sample_data<-Read10X(data.dir = sample)
	# project should be different for each sample
	sample_data<-CreateSeuratObject(counts = sample_data, project = sample, assay = "RNA")

	# data QC
	sample_data[["percent.mt"]]<-PercentageFeatureSet(sample_data, pattern = "^MT-")
	sample_data<-SubsetData(sample_data, subset.name = "nCount_RNA", 
		low.threshold = low_cutoff[sample], high.threshold = high_cutoff[sample])
	sample_data<-SubsetData(sample_data, subset.name = "percent.mt", high.threshold = 5)
})

# Normalization and variable feature selection
data_list<-lapply(data_list, FUN = function(sample_data) {
	sample_data<-NormalizeData(sample_data, normalization.method = "LogNormalize")
	sample_data<-FindVariableFeatures(sample_data, selection.method = "vst", nfeatures = 2000)
})

# Integrate data
# Integrated data has a cell-level meta data (comb_data@meta.data) named "orig.ident" storing the sample where the cell is from.
# Cell barcode will also be modified to add sample index.
anchors<-FindIntegrationAnchors(object.list = data_list, dims = 1:20)
comb_data<-IntegrateData(anchorset = anchors, dims = 1:20)

DefaultAssay(comb_data)<-"integrated"

# Dimensional reduction
comb_data<-ScaleData(comb_data)
comb_data<-RunPCA(comb_data)
comb_data<-RunUMAP(comb_data, reduction = "pca", dims = 1:30)

# Cell cluster is also the default identity (comb_data@active.ident), which can be get or re-set by Idents(comb_data). 
comb_data<-FindNeighbors(comb_data, reduction = "pca", dims = 1:30)
# Smaller resolution, fewer clusters
comb_data<-FindClusters(comb_data, resolution = 0.1)

# Visualization
pdf("umap_0.1.pdf")
DimPlot(comb_data, reduction = "umap", label = TRUE)
dev.off()

# Plot marker genes
DefaultAssay(comb_data)<-"RNA"
pdf("markers.pdf")
DotPlot(comb_data, features = rev(c("CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "NCAM1", "KLRB1", "NKG7", "CD19", "MS4A1", "CD38", "CD14", "CD68", "FCGR3A", "LILRA4", "PPBP", "HBB"))) + RotatedAxis()
dev.off()

# Get finer clusters
comb_data[["clusters_0.1"]]<-comb_data[["seurat_clusters"]]
comb_data<-FindClusters(comb_data, resolution = 1.2)
comb_data[["clusters_1.2"]]<-comb_data[["seurat_clusters"]]

pdf("umap_1.2.pdf")
DimPlot(comb_data, reduction = "umap", label = TRUE)
dev.off()

save(list = ls(all.names = TRUE), file = "batch1.Rdata")

# Sum cell count of each sample and cluster for differential abundance analysis
sample_celltype<-paste(comb_data[["orig.ident"]][,1], comb_data[["clusters_1.2"]][,1], sep = "_")
sample_celltype_sum<-table(sample_celltype)
write.table(x = data.frame(sample_celltype_sum), file = "batch1_cellcount.txt", row.names = F, quote = F)


# Reset cell identity by combining sample and cell cluster
sample_celltype<-paste(comb_data[["orig.ident"]][,1], comb_data[["clusters_0.1"]][,1], sep = "_")
Idents(comb_data)<-sample_celltype

# Get raw UMI counts (not normalized) aggregated for each cell type (pseudo-bulk for DEG detection)
avg_count<-AverageExpression(comb_data, assays = "RNA", use.counts = T)$RNA
cell_num<-table(Idents(comb_data))
total_count<-t(t(as.matrix(avg_count)) * as.numeric(cell_num))
write.table(file = "batch1_UMIcount.txt", x = total_count, sep = "\t", quote = F)

