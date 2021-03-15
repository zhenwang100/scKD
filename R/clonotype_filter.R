# This script extract productive and paired IG chain
# Input from cellranger: consensus_annotations.csv, filtered_contig_annotations.csv

# Clonotype consensus
consensus_tab<-read.csv("consensus_annotations.csv")

# Select productive consensus
consensus_tab<-consensus_tab[consensus_tab$productive == "True",]

# Consider IGH and IGL/IGK, respectively 
consensus_IGH<-consensus_tab[consensus_tab$chain == "IGH",]
consensus_IGL<-consensus_tab[consensus_tab$chain == "IGL" | consensus_tab$chain == "IGK",]

# If a clonotype has more than one IGH chains, select the one with most UMI supported
rep_IGH<-tapply(consensus_IGH$umis, INDEX = consensus_IGH$clonotype_id, FUN = max)
consensus_IGH<-merge(consensus_IGH, data.frame(clonotype_id = names(rep_IGH), umis = rep_IGH))

# If a clonotype has more than one IGL/IGK chains, select the one with most UMI supported
rep_IGL<-tapply(consensus_IGL$umis, INDEX = consensus_IGL$clonotype_id, FUN = max)
consensus_IGL<-merge(consensus_IGL, data.frame(clonotype_id = names(rep_IGL), umis = rep_IGL))

# Only preserve clonotypes with both IGH and IGK/IGL chains
pair_IG<-intersect(consensus_IGH$clonotype_id, consensus_IGL$clonotype_id)
consensus_IGH<-consensus_IGH[consensus_IGH$clonotype_id %in% pair_IG,]
consensus_IGL<-consensus_IGL[consensus_IGL$clonotype_id %in% pair_IG,]

# Preserve only cells with productive and paired IG chain
contig_tab<-read.csv("filtered_contig_annotations.csv")
contig_tab<-contig_tab[contig_tab$raw_consensus_id %in% c(consensus_IGH$consensus_id, consensus_IGL$consensus_id),]

write.csv(contig_tab, file = "filtered_contig_annotations_productive_pair.csv", quote = F, row.names = F)
 
