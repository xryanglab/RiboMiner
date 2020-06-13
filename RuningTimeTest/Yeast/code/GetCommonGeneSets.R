up_genes <- read.table("../results/MetageneAnalysis/eIF5A_up_CDS_normed_transcript_id.txt",sep="\t",stringsAsFactors = F,header = T,check.names = F)
up_genes_common <- as.data.frame(Reduce(intersect,up_genes))

down_genes <- read.table("../results/MetageneAnalysis/eIF5A_down_CDS_normed_transcript_id.txt",sep="\t",stringsAsFactors = F,header = T,check.names = F)
down_genes_common <- as.data.frame(Reduce(intersect,down_genes))

unblocked_genes <- read.table("../results/MetageneAnalysis/eIF5A_unblocked_CDS_normed_transcript_id.txt",sep="\t",stringsAsFactors = F,header = T,check.names = F)
unblocked_genes_common <- as.data.frame(Reduce(intersect,unblocked_genes))

colnames(up_genes_common) <- "trans_id"
colnames(down_genes_common) <- "trans_id"
colnames(unblocked_genes_common) <- "trans_id"

write.table(up_genes_common,"../results/up_common_trans.txt",quote = F,sep="\t",row.names = F,col.names = T)
write.table(down_genes_common,"../results/down_common_trans.txt",quote = F,sep="\t",row.names = F,col.names = T)
write.table(unblocked_genes_common,"../results/unblocked_common_trans.txt",quote = F,sep="\t",row.names = F,col.names = T)
