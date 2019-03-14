r = read.table("norm_tum_zyg.txt", header=T, sep = "\t",stringsAsFactors=FALSE)
r1 <- subset(r[,1:11])
r2 <- subset(r[,12:21])
R1 <- r1[,11] #ctrl tot
R2 <- r2[,10] #case tot
N <- R1 +R2
N11 <- r1[,8]
N12 <- r1[,9]
N13 <- r1[,10]
N21 <- r2[,7]
N22 <- r2[,8]
N23 <- r2[,9]

C1 <- N11+N21
C2 <- N12 + N22
C3 <- N13 + N23
T1 <- 0*(N11*R2-N21*R1) + 1*(N12*R2-N22*R1)+2*(N13*R2-N23*R1)
VT <- (R1*R2/N)*(((C2*(N-C2))+(4*C3*(N-C3))) - (4*C2*C3))
CT <- T1^2/(VT)
pval <-  1-pchisq(CT,df=1)
padj <- p.adjust(pval, method = "bonferroni" , n = length(pval))
gwas <- cbind(r[,1:6],R1,R2,N11,N12,N13,N21,N22,N23,C1,C2,C3,T1,VT,CT,pval,padj)


write.table(gwas,"Tum_gwas.txt",sep="\t",quote=F,row.names=F)
