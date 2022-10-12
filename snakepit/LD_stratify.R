#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

lds_seg = read.table(args[1],header=T,colClasses=c("character",rep("numeric",8)))
quartiles=summary(lds_seg$ldscore_SNP)

lb1 = which(lds_seg$ldscore_SNP <= quartiles[2])
lb2 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP <= quartiles[3])
lb3 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP <= quartiles[5])
lb4 = which(lds_seg$ldscore_SNP > quartiles[5])

lb1_snp = lds_seg$SNP[lb1]
lb2_snp = lds_seg$SNP[lb2]
lb3_snp = lds_seg$SNP[lb3]
lb4_snp = lds_seg$SNP[lb4]

write.table(lb1_snp, paste0(args[2],"LD_group1.ids"), row.names=F, quote=F, col.names=F)
write.table(lb2_snp, paste0(args[2],"LD_group2.ids"), row.names=F, quote=F, col.names=F)
write.table(lb3_snp, paste0(args[2],"LD_group3.ids"), row.names=F, quote=F, col.names=F)
write.table(lb4_snp, paste0(args[2],"LD_group4.ids"), row.names=F, quote=F, col.names=F)
