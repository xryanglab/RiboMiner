all_mean_density <- read.table("../results/MetageneAnalysis/eIF5A_all_local_mean_density.txt",sep="\t",stringsAsFactors = F,header = T,check.names = F)
all_mean_density_non_zeros <- all_mean_density[Reduce(intersect,lapply(all_mean_density,function(x){which(x>0)})),]
#all_mean_density_non_zeros['si-Ctrl'] <- (all_mean_density_non_zeros$`si-Ctrl-1`+all_mean_density_non_zeros$`si-Ctrl-2`)/2
all_mean_density_non_zeros['si-Ctrl'] <- all_mean_density_non_zeros$`si-Ctrl-2`
all_mean_density_non_zeros['si-eIF5A'] <- all_mean_density_non_zeros$`si-eIF5A-2`
#all_mean_density_non_zeros['si-eIF5A'] <- (all_mean_density_non_zeros$`si-eIF5A-1`+all_mean_density_non_zeros$`si-eIF5A-2`)/2
all_mean_density_non_zeros['FC(si-eIF5A/si-Ctrl)'] <- all_mean_density_non_zeros$`si-eIF5A`/all_mean_density_non_zeros$`si-Ctrl`
colnames(all_mean_density_non_zeros)[1] <- "trans_id"

FC_le_0.67 <- all_mean_density_non_zeros[all_mean_density_non_zeros$`FC(si-eIF5A/si-Ctrl)`<=0.67,]
FC_ge_1.5 <- all_mean_density_non_zeros[all_mean_density_non_zeros$`FC(si-eIF5A/si-Ctrl)`>=1.5,]
FC_lt1.5_gt0.67 <- all_mean_density_non_zeros[all_mean_density_non_zeros$`FC(si-eIF5A/si-Ctrl)`<1.5&all_mean_density_non_zeros$`FC(si-eIF5A/si-Ctrl)`>0.67,]

write.table(FC_le_0.67,"../results/MetageneAnalysis/down_trans.txt",quote = F,sep="\t",row.names = F,col.names = T)
write.table(FC_ge_1.5,"../results/MetageneAnalysis/up_trans.txt",quote = F,sep="\t",row.names = F,col.names = T)
write.table(FC_lt1.5_gt0.67,"../results/MetageneAnalysis/unblocked_trans.txt",quote = F,sep="\t",row.names = F,col.names = T)

