# Summarize cell abundance by types for visiualization and statistical test
library(ggplot2)

# Extract cell compartments for analysis
cell<-"monocyte"
#cell<-"B"
#cell<-"CD4 T"
#cell<-"CD8 T"
#cell<-"NK"

# Combine cell counts across batches
batch1<-read.table("batch1_cellcount.txt", header = T, sep = "\t")
batch2<-read.table("batch2_cellcount.txt", header = T, sep = "\t")
data<-rbind(batch1, batch2)

# Manually refined cell cluster annotations
ann<-read.table("cell_anno.txt", header = T, sep = "\t")
data<-merge(data, ann, by = "Cluster")

# Filtering cell types
data<-subset(data, grepl(cell, Type))

# Drop filtered levels
data<-droplevels(data)

# Sum count for a cell type for each sample
type.count<-aggregate(x = data['Count'], by = data[c('Patient', 'Treatment', 'Type')], sum)
# Total count for a sample
total.count<-aggregate(x = data['Count'], by = data[c('Patient', 'Treatment')], sum)

# Merge cell type count and total count to calculate cell type proportion
merge.count<-merge(type.count, total.count, by = c("Patient", "Treatment"), suffixes = c(".type",".total"))
merge.count$Prop<-merge.count$Count.type / merge.count$Count.total

# Boxplot cell type proportion for each sample group
plot<-ggplot(merge.count, aes(x = Treatment, y = Prop * 100)) +
	scale_x_discrete(limits = c("before", "after", "health"), name = NULL) +
	ylab(paste("% of", cell, "cells", sep = " ")) +
	theme_bw() +	
	facet_wrap(~Type, scales = "free") +
	geom_boxplot()+
	geom_point(aes(colour = Patient)) +
	geom_line(aes(group = Patient, colour = Patient))

# Save plot
ggsave(paste(cell, "abundance.pdf", sep = "_"), plot)

# Perform statistical test
# Convert the data format to matrix (cell type * sample)
merge.count$Sample<-paste(merge.count$Patient, merge.count$Treatment, sep = "_")
prop.matrix<-t(tapply(merge.count$Prop, INDEX = list(merge.count$Sample, merge.count$Type), sum))

# Get p value for each cell type between groups
p.value<-apply(prop.matrix, 1, FUN = function(x) {
	# Group index
	health<-grepl("health", names(x))
	before<-grepl("before", names(x))
	after<-grepl("after", names(x))

	# Pairwise p value
	p1<-t.test(x = x[health], y = x[before])$p.value
	p2<-t.test(x = x[health], y = x[after])$p.value
	p3<-t.test(x = x[before], y = x[after], paired = T)$p.value
	
	p<-c(p1, p2, p3)
	names(p)<-c("health_before", "health_after", "before_after")
	return(p)
})

