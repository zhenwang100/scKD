# Perform cell type annotation for Seurat clusters using SingleR
# https://github.com/LTLA/SingleR
library(SingleR)
library(Seurat)
library(ggplot2)

# Load Seurat object
load("batch1.Rdata")

# Extract log-count expression matrix
exp<-GetAssayData(comb_data, slot = "data", assay = "RNA")

# Load internal reference of human immune cells for SingleR
# https://bioconductor.org/packages/3.12/data/experiment/vignettes/celldex/inst/doc/userguide.html
refs<-list(ref1 = HumanPrimaryCellAtlasData(), 
	ref2 = BlueprintEncodeData(),
	ref3 = DatabaseImmuneCellExpressionData(),
	ref4 = NovershternHematopoieticData(), 
	ref5 = MonacoImmuneData())

# Perform annotation using each reference, respectively
sapply(names(refs), FUN = function(ref_idx) {
	ref<-refs[[ref_idx]]
	ref_name<-paste("SingleR", ref_idx, sep = ".")

	# Annotation is performed on cluster-level rather than default single-cell level
	rst<-SingleR(test = exp, ref = ref, method = "cluster", clusters = comb_data[["clusters_1.2"]], labels = ref$label.fine)

	# Assign predicted cell labels to seurat object
	# SingleR assigns labels for all clusters, even for those not in the reference
	# So use pruned.labels to remove uncertain labels 
	comb_data[[ref_name]]<-rst$pruned.labels[match(comb_data[["clusters_1.2"]], rownames(rst))]
	
	# Visualize cell labels
	p<-DimPlot(comb_data, reduction = "umap", group.by = ref_name, label = TRUE)
	ggsave(filename = paste(ref_name, "pdf", sep = "."), plot = p, width = 10, height = 7)	
})

