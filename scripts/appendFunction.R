args = commandArgs(trailingOnly=TRUE)
df1<-read.delim(args[1],header=F)
df2<-read.delim(args[2],header=T)
df3<-read.csv(args[3],header=T)
colnames(df1) = c("genename","GOTerm","probability")
df3df1 = merge(df3,df1,by="genename")
dfcom = df3df1[ which(df3df1$probability > as.numeric(args[4])),]
dfcom1 = merge(dfcom,df2,by="GOTerm")
dfcom2 = unique(dfcom1)
dfcom3 = dfcom2[,c(1,2,3,5,6)]
colnames(dfcom3) = c("GOTerm","genename","probability","function","type")
write.table(dfcom3,quote=F,row.names=F,sep="\t",args[5])