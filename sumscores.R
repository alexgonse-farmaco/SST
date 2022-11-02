####Sum chr scores

library(data.table)
library(dplyr)

profilepath="/external/rprshnas01/kcni/asegura/SST_Fernanda/results"
profile01=grep(list.files(path=profilepath),pattern = ".qv0.1.profile", value = T)
profile1=grep(list.files(path=profilepath),pattern = ".qv1.profile", value = T)

tempID=fread('SSTPRS_cmcimp_chr22.qv0.1.profile')
tempcbind=tempID%>%dplyr::select(.,IID)
tempID=tempID%>%dplyr::select(.,IID)
tempsum=data.frame()

for (chrom in profile01){
  chr01=fread(chrom)
  chr01=chr01%>%dplyr::select(.,SCORESUM)
  tempcbind=cbind.data.frame(tempcbind,chr01)
}
tempsum01=rowSums(tempcbind[,-1])
prof01=cbind.data.frame(tempID,tempsum01)

tempID=fread('SSTPRS_cmcimp_chr22.qv1.profile')
tempcbind=tempID%>%dplyr::select(.,IID)
tempID=tempID%>%dplyr::select(.,IID)
tempsum=data.frame()

for (chrom in profile1){
  chr1=fread(chrom)
  chr1=chr01%>%dplyr::select(.,SCORESUM)
  tempcbind=cbind.data.frame(tempcbind,chr1)
}

tempsum1=rowSums(tempcbind[,-1])
prof1=cbind.data.frame(tempID,tempsum1)

prof=merge(prof01,prof1,by='IID')
colnames(prof)=c('ID','SST_PRS_0.1','SST_PRS_1')
fwrite(prof,'CMC_PRS_SST.txt',sep='\t')
