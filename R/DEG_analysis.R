# Performing DEG for each cell compartment with DESeq2
library(DESeq2)

# Combine UMI count across batches
batch1<-read.table("batch1_UMIcount.txt", header = T, row.names = 1)
batch2<-read.table("batch2_UMIcount.txt", header = T, row.names = 1)
data<-cbind(batch1, batch2)

# Cluster ids for analysis
cluster<-c(4)		# monocyte
#cluster<-c(1, 10)	# B cell
#cluster<-c(0, 3, 6)	# CD4 cell 
#cluster<-c(2, 7, 8)	# CD8 cell
#cluster<-c(5)		# NK
#cluster<-c(0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 13)	# All PBMC

##### PART I. Compare before and after #####
sample1<-c("S1_1224", "S1_1229", "S2_0317", "S2_0320", "S3_0414", "S3_0416", "S4_0601", "S4_0603", 
	   "P5_1112", "P5_1114", "P6_1120", "P6_1123", "P7_1123", "P7_1125")
meta1<-data.frame(patient = c("P1", "P1", "P2", "P2", "P3", "P3", "P4", "P4", "P5", "P5", "P6", "P6", "P7", "P7"), 
		  treatment = c("before", "after", "before", "after", "before", "after", "before", "after", 
				"before", "after", "before", "after", "before", "after"))

# Retrive cluster specific expression profile
col_name1<-paste(sample1, cluster[1], sep = "_")
subdata1<-data[,col_name1]

# If more than two cluster ids are assigned, add their expression profile
if (length(cluster) > 1) {	
	for (i in 2:length(cluster)) {
		col_name1<-paste(sample1, cluster[i], sep = "_")
		subdata1<-subdata1 + data[,col_name1]
	}
}

# Construct DESeq dataset
# paired design
dds1<-DESeqDataSetFromMatrix(countData = subdata1, colData = meta1, design = ~patient + treatment)

# Normalized by CPM/FPM
norm1<-fpm(dds1)
write.table(x = cbind(NAME = rownames(norm1), DESCRIPTION = "na", norm1), file = "CPM1.txt", 
	    sep = "\t", quote = F, row.names = F)	# GSEA convention

# Pre-filtering low expression genes
# Only minimal pre-filtering is needed here because DESeq2 will automatically apply strict filtering in results function
keep1<-rowSums(counts(dds1)) >= 10
dds1<-dds1[keep1,]

# Perform DE analysis
dds1<-DESeq(dds1)

# Get result (default for p value of the last variable in the design formula)
res1<-results(dds1)

# Summary result with FDR < 0.05
summary(res1, alpha = 0.05)

# Log fold change shrinkage (useful for ranking of fold change, but do not change p or FDR)
#resLFC1<-lfcShrink(dds1, coef = "treatment_before_vs_after")
#plotMA(res1)
#plotMA(resLFC1)

# Output results
write.table(as.data.frame(res1), file = "stat1.txt", sep = "\t", quote = F)	# for clusterProfiler analysis

##### PART II. Compare before and health #####
sample2<-c("S1_1224", "S2_0317", "S3_0414", "S4_0601", "P5_1112", "P6_1120", "P7_1123", 
	   "H1_0601", "H2_0601", "H3_0601")
meta2<-data.frame(patient = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "H1", "H2", "H3"), 
		  treatment = c("before", "before", "before", "before", "before", "before", "before",
				"health", "health", "health"))

# Retrive cluster specific expression profile
col_name2<-paste(sample2, cluster[1], sep = "_")
subdata2<-data[,col_name2]
if (length(cluster) > 1) {	
	for (i in 2:length(cluster)) {
		col_name2<-paste(sample2, cluster[i], sep = "_")
		subdata2<-subdata2 + data[,col_name2]
	}
}

dds2<-DESeqDataSetFromMatrix(countData = subdata2, colData = meta2, design = ~treatment)

# Normalized by CPM/FPM
norm2<-fpm(dds2)
write.table(x = cbind(NAME = rownames(norm2), DESCRIPTION = "na", norm2), file = "CPM2.txt", 
	    sep = "\t", quote = F, row.names = F)	# GSEA convention

# Pre-filtering low expression genes
# Only minimal pre-filtering is needed here because DESeq2 will automatically apply strict filtering in results function
keep2<-rowSums(counts(dds2)) >= 10
dds2<-dds2[keep2,]

# Perform DE analysis
dds2<-DESeq(dds2)

# Get result (default for p value of the last variable in the design formula)
res2<-results(dds2)

# Summary result with FDR < 0.05
summary(res2, alpha = 0.05)

# Output results
write.table(as.data.frame(res2), file = "stat2.txt", sep = "\t", quote = F)
