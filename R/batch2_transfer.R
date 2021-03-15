# Transfer cell labels from Batch 1 to Batch 2 samples
# Seurat 3.0.2
library(Seurat)

load("batch1.Rdata")
DefaultAssay(comb_data)<-"integrated"

# Samples of Batch 2
samples<-c("P5_1112", "P5_1114", "P6_1120", "P6_1123", "P7_1123", "P7_1125")
low_cutoff<-c(2000, 2000, 2000, 2000, 2000, 2000)
high_cutoff<-c(60000, 60000, 60000, 60000, 60000, 60000)

names(low_cutoff)<-samples
names(high_cutoff)<-samples

# Read data and data QC
data_list<-lapply(samples, FUN = function(sample) {
	sample_data<-Read10X(data.dir = sample)
	# project should be different for each sample
	sample_data<-CreateSeuratObject(counts = sample_data, project = sample, assay = "RNA")

	# Quality control
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

# Transfer cell labels from Batch 1
data_list<-lapply(data_list, FUN = function(sample_data) {
	trans_anchors<-FindTransferAnchors(reference = comb_data, query = sample_data)
	pred_1<-TransferData(anchorset = trans_anchors, refdata = comb_data$clusters_0.1)
	pred_2<-TransferData(anchorset = trans_anchors, refdata = comb_data$clusters_1.2)

	sample_data$clusters_0.1<-pred_1$predicted.id
	sample_data$clusters_1.2<-pred_2$predicted.id

	return(sample_data)
})

# Merge Batch 2 samples
merge_data<-merge(x = data_list[[1]], y = data_list[2:length(data_list)])

save(list = c("data_list", "merge_data"), file = "batch2.Rdata")

# Sum cell count for each sample and cluster
sample_celltype<-paste(merge_data[["orig.ident"]][,1], merge_data[["clusters_1.2"]][,1], sep = "_")
sample_celltype_sum<-table(sample_celltype)
write.table(x = data.frame(sample_celltype_sum), file = "batch2_cellcount.txt", row.names = F, quote = F)

# Reset cell identity by combining sample and cell cluster
sample_celltype<-paste(merge_data[["orig.ident"]][,1], merge_data[["clusters_0.1"]][,1], sep = "_")
Idents(merge_data)<-sample_celltype

# Get raw UMI counts (not normalized) aggregated for each cell type (pseudo-bulk for DEG detection)
avg_count<-AverageExpression(merge_data, assays = "RNA", use.counts = T)$RNA
cell_num<-table(Idents(merge_data))
total_count<-t(t(as.matrix(avg_count)) * as.numeric(cell_num))
write.table(file = "batch2_UMIcount.txt", x = total_count, sep = "\t", quote = F)

