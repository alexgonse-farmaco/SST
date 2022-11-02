bmikg2=read.table("sumstat_QCed.txt",header=T)

mergd0=data.frame()
mergd=data.frame()
for (i in 1:22){
  loadname=paste0('target_clumped_1_chr',i,'.clumped')
  clumped=read.table(loadname,header=T)
  bmikg3=merge(bmikg2,clumped,by="SNP",sort=F)
  bmikg3a1=subset(bmikg3,select=c("SNP","V5","W"))
  mergd0=rbind.data.frame(mergd0,bmikg3)
  mergd=rbind.data.frame(mergd,bmikg3a1)
}

colnames(mergd)<-c("SNP","A1","W")
write.table(mergd,"target_clumped_1_nothreshold.raw",col.names=T,row.names=F,quote=F,sep='\t')
mergd2=subset(mergd0,select=c("SNP","P.x"))
write.table(mergd2,"target_clumped_1_nothreshold.snppvalues",col.names=T,row.names=F,quote=F,sep='\t')

q()
n
