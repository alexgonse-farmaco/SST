
kg=read.table('target.bim',header=F, fill = T)
bmi=read.table('sumstat.txt',header=T)

head(kg)
nrow(kg)

head(bmi)
nrow(bmi)
## Attribute new column names to the sumstat file

colnames(bmi)<-c("SNP", "A1", "A2", "Weight", "P")

bmikg=merge(bmi,kg,by.x="SNP",by.y="V2")      # the merge step

nrow(bmikg)
head(bmikg)
#### See column names to confirm


### To mark ambiguous G/C or A/T SNPs for removal

## let’s create a column flagging the ambiguous SNPs

bmikg$ambiguousSNPs[bmikg$A1 == "A" & bmikg$A2 == "T"]<-1
bmikg$ambiguousSNPs[bmikg$A1 == "C" & bmikg$A2 == "G"]<-1
bmikg$ambiguousSNPs[bmikg$A1 == "G" & bmikg$A2 == "C"]<-1
bmikg$ambiguousSNPs[bmikg$A1 == "T" & bmikg$A2 == "A"]<-1

head(bmikg)


### To remove ambiguous SNPs

bmikg1=bmikg[is.na(bmikg$ambiguousSNPs),] # to remove ambiguous SNPs

nrow(bmikg1)      # number of non-unambiguous SNPs remained


### To find perfect allele matches (beta will be the same) 

bmikg1aa=bmikg1[as.character(bmikg1$A1) == bmikg1$V5 & as.character(bmikg1$A2) == bmikg1$V6,]

nrow(bmikg1aa)


### To find allele matches for flipped strand -- beta will be the same

bmikg1ab=bmikg1[bmikg1$A1 == "A" & bmikg1$A2 == "C" & bmikg1$V5 == "T" & bmikg1$V6 == "G",]
bmikg1ac=bmikg1[bmikg1$A1 == "A" & bmikg1$A2 == "G" & bmikg1$V5 == "T" & bmikg1$V6 == "C",]
bmikg1ad=bmikg1[bmikg1$A1 == "C" & bmikg1$A2 == "A" & bmikg1$V5 == "G" & bmikg1$V6 == "T",]
bmikg1ae=bmikg1[bmikg1$A1 == "C" & bmikg1$A2 == "T" & bmikg1$V5 == "G" & bmikg1$V6 == "A",]
bmikg1af=bmikg1[bmikg1$A1 == "G" & bmikg1$A2 == "A" & bmikg1$V5 == "C" & bmikg1$V6 == "T",]
bmikg1ag=bmikg1[bmikg1$A1 == "G" & bmikg1$A2 == "T" & bmikg1$V5 == "C" & bmikg1$V6 == "A",]
bmikg1ah=bmikg1[bmikg1$A1 == "T" & bmikg1$A2 == "C" & bmikg1$V5 == "A" & bmikg1$V6 == "G",]
bmikg1ai=bmikg1[bmikg1$A1 == "T" & bmikg1$A2 == "G" & bmikg1$V5 == "A" & bmikg1$V6 == "C",]
bmikg1a=rbind(bmikg1aa,bmikg1ab,bmikg1ac,bmikg1ad,bmikg1ae,bmikg1af,bmikg1ag,bmikg1ah,bmikg1ai)   # to combine all the datasets created above + bmikg1aa with non-ambiguous SNPs

nrow(bmikg1a)


### To add column W for Weight to be used in PRS; W=Weight for matching SNPs and W=-Weight for non-matching SNPs

bmikg1a$W=bmikg1a$Weight

### To find perfect allele switches between A1 and A2

bmikg1ba=bmikg1[as.character(bmikg1$A1) == bmikg1$V6 & as.character(bmikg1$A2) == bmikg1$V5,]

nrow(bmikg1ba)


### To find flipped strands but perfect allele switches between A1 and A2

bmikg1bb=bmikg1[bmikg1$A1 == "A" & bmikg1$A2 == "C" & bmikg1$V5 == "G" & bmikg1$V6 == "T",]
bmikg1bc=bmikg1[bmikg1$A1 == "A" & bmikg1$A2 == "G" & bmikg1$V5 == "C" & bmikg1$V6 == "T",]
bmikg1bd=bmikg1[bmikg1$A1 == "C" & bmikg1$A2 == "A" & bmikg1$V5 == "T" & bmikg1$V6 == "G",]
bmikg1be=bmikg1[bmikg1$A1 == "C" & bmikg1$A2 == "T" & bmikg1$V5 == "A" & bmikg1$V6 == "G",]
bmikg1bf=bmikg1[bmikg1$A1 == "G" & bmikg1$A2 == "A" & bmikg1$V5 == "T" & bmikg1$V6 == "C",]
bmikg1bg=bmikg1[bmikg1$A1 == "G" & bmikg1$A2 == "T" & bmikg1$V5 == "A" & bmikg1$V6 == "C",]
bmikg1bh=bmikg1[bmikg1$A1 == "T" & bmikg1$A2 == "C" & bmikg1$V5 == "G" & bmikg1$V6 == "A",]
bmikg1bi=bmikg1[bmikg1$A1 == "T" & bmikg1$A2 == "G" & bmikg1$V5 == "C" & bmikg1$V6 == "A",]
bmikg1b=rbind(bmikg1ba,bmikg1bb,bmikg1bc,bmikg1bd,bmikg1be,bmikg1bf,bmikg1bg,bmikg1bh,bmikg1bi)

nrow(bmikg1b)


### NOTE: there are 11 SNPs with mismatching A1 and A2 between Base GWAS summary stat and target GWAS data

### NOTE: there are 0 INDELs with mismatching A1 and A2 between Base GWAS summary stat and target GWAS data

nrow(bmikg1b)+nrow(bmikg1a)

### W = -Weight for SNPs with switched A1 and A2
bmikg1b$W=0-bmikg1b$Weight
### ***NOTE: Use V5 as A1 & V6 as A2 & W as the effect estimate for each A1 in genotype*** 

bmikg2=rbind(bmikg1a,bmikg1b)

nrow(bmikg2) 

head(bmikg2)

write.table(bmikg2,"sumstat_QCed.txt",quote=F,sep='\t',col.names=T,row.names=F)

bmikg2s=subset(bmikg2,select=c("SNP","P"))
colnames(bmikg2s)<-c("SNP","P")
write.table(bmikg2s,"SNPs_P_forclumping.txt",quote=F,sep='\t',col.names=T,row.names=F)
q()
