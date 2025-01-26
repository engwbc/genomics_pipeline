args<-commandArgs(trailingOnly=TRUE)
imiss <- args[1]
genome <- args[2]
d1<-read.table(imiss, header=T)
d2<-read.table(genome, header=T)
rownames(d1) = d1$IID 
out = NULL
for (i in 1:nrow(d2)){
 id1 = as.character(d2[i, "IID1"])
 id2 = as.character(d2[i, "IID2"])
 
 if (d1[id1, "F_MISS"] > d1[id2, "F_MISS"]){
  out =c(out, id1)
 }
 else{
  out =c(out, id2)
 }  

}
data<-as.data.frame(out)
data[,2]<-data[,1]
write.table(data,file= "", quote=F, col.names = F, row.names = F)
