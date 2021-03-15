# GO and KEGG enrichment with clusterProfiler
library(clusterProfiler)
library(ggplot2)

# Input file (before vs after and health vs before, respectively)
stat1<-read.table("D:/KD/Integrate/DEG/NK/stat1.txt", header = T, sep = "\t")
stat2<-read.table("D:/KD/Integrate/DEG/NK/stat2.txt", header = T, sep = "\t")

genecut<-0.05	# DEG cut off
enrichcut<-0.05	# enrichment cut off

# Significant DEGs
sig1<-subset(stat1, padj < genecut)
sig2<-subset(stat2, padj < genecut)

symbols1<-rownames(sig1)
symbols2<-rownames(sig2)

map1<-bitr(symbols1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
map2<-bitr(symbols2, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

id1<-map1$ENTREZID
id2<-map2$ENTREZID

##### PART I. GO enrichment #####
# GO enrichment for one group
#ego1<-enrichGO(gene = id1, OrgDb = org.Hs.eg.db, ont = "MF", pvalueCutoff = enrichcut)
#ego1<-simplify(ego1)	# Remove redundancy
#ego1<-setReadable(ego1, 'org.Hs.eg.db', 'ENTREZID')	# convert ID to symbol
#as.data.frame(ego1)

# Compare two groups (MF)
ego.mf<-compareCluster(geneCluster = list("before vs after" = id1, "before vs health" = id2), fun = "enrichGO", 
	OrgDb = org.Hs.eg.db, ont = "MF", pvalueCutoff = enrichcut)
ego.mf<-simplify(ego.mf)
ego.mf<-setReadable(ego.mf, 'org.Hs.eg.db', 'ENTREZID')

write.table(x = as.data.frame(ego.mf), file = "ego_MF.txt", row.names = F, quote = F, sep = "\t")

# Compare two groups (BP)
ego.bp<-compareCluster(geneCluster = list("before vs after" = id1, "before vs health" = id2), fun = "enrichGO", 
	OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = enrichcut)
ego.bp<-simplify(ego.bp)
ego.bp<-setReadable(ego.bp, 'org.Hs.eg.db', 'ENTREZID')

write.table(x = as.data.frame(ego.bp), file = "ego_BP.txt", row.names = F, quote = F, sep = "\t")

##### PART II. KEGG enrichment #####
ekegg1<-enrichKEGG(gene = id1, organism = 'hsa', pvalueCutoff = enrichcut)
ekegg1<-setReadable(ekegg1, 'org.Hs.eg.db', 'ENTREZID')
ekegg1<-as.data.frame(ekegg1)
if (nrow(ekegg1) > 0) ekegg1<-data.frame(Cluster = "before vs after", ekegg1)

ekegg2<-enrichKEGG(gene = id2, organism = 'hsa', pvalueCutoff = enrichcut)
ekegg2<-setReadable(ekegg2, 'org.Hs.eg.db', 'ENTREZID')
ekegg2<-as.data.frame(ekegg2)
if (nrow(ekegg2) > 0) ekegg2<-data.frame(Cluster = "before vs health", ekegg2)

write.table(x = rbind(ekegg1, ekegg2), 
	file = "ekegg.txt", row.names = F, quote = F, sep = "\t")

##### PART III. Hallmark gene sets enrichment ##### 
# Read hallmark gmt file
gmt<-read.gmt("D:/GSEA_Database/h.all.v7.0.symbols.gmt")

egmt1<-enricher(symbols1, pvalueCutoff = enrichcut, TERM2GENE = gmt)
egmt2<-enricher(symbols2, pvalueCutoff = enrichcut, TERM2GENE = gmt)

egmt1<-as.data.frame(egmt1)
if (nrow(egmt1) > 0) egmt1<-data.frame(Cluster = "before vs after", egmt1)

egmt2<-as.data.frame(egmt2)
if (nrow(egmt2) > 0) egmt2<-data.frame(Cluster = "before vs health", egmt2)

write.table(x = rbind(egmt1, egmt2), file = "egmt.txt", row.names = F, quote = F, sep = "\t")

