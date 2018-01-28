#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
a1<-read.delim(args[1],header=T,check.names=F)
a2<-read.delim(args[2],header=T,check.names=F)
g1dfgi1<-a1[,c(2:length(a1))]
g1dfgi1_t<-t(g1dfgi1)
n<-ncol(g1dfgi1_t)
set3<-numeric(0)
for (x in 1:n){
    normcol1<-t((g1dfgi1_t[,x])/(max(t(g1dfgi1)[,x])))
    set3<-rbind(set3,normcol1)
}
colg1<-as.vector(a1$genename)
set3t<-data.frame(t(set3))
colnames(set3t)<-a1$genename
g1dfgi1<-a2[,c(2:length(a2))]
g1dfgi1_t<-t(g1dfgi1)
n<-ncol(g1dfgi1_t)
set4<-numeric(0)
for (x in 1:n){
    normcol1<-t((g1dfgi1_t[,x])/(max(t(g1dfgi1)[,x])))
    set4<-rbind(set4,normcol1)
}
colg1<-as.vector(a2$genename)
set4t<-data.frame(t(set4))
colnames(set4t)<-a2$genename
write.csv(set3t,quote=F,row.names=F,"set3t-full.csv")
write.csv(set4t,quote=F,row.names=F,"set4t-full.csv")