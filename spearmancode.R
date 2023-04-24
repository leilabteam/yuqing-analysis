#dia second
data.dia.second.mean <-read.csv("TMT4mean.csv",row.names = 1,stringsAsFactors = F)
data.dia.second.mean<-data.dia.second.mean[!nchar(rownames(data.dia.second.mean)) > 20,]
#matrix.rowname<-rownames(data.dia.second.mean)
#matrix.rowname<- t(as.data.frame(strsplit(matrix.rowname,split = ".",fixed = T)))
#rownames(data.dia.second.mean)<- as.vector(matrix.rowname[,1])

#dia second
data.TMT.mean <-read.csv("mRNA4mean.csv",row.names = 1,stringsAsFactors = F)
data.TMT.mean<-data.TMT.mean[!nchar(rownames(data.TMT.mean)) > 20,]
#matrix.rowname<-rownames(data.TMT.mean)
#matrix.rowname<- t(as.data.frame(strsplit(matrix.rowname,split = ".",fixed = T)))
#rownames(data.TMT.mean)<- as.vector(matrix.rowname[,1])

gene.cor.cal.from.two.matix<-function(matrix1,matrix2,outname){
  tmp.c<- c()
  tmp.c.pvalue<-c()
  #计算交集基因：
  intersect.gene.name<- intersect(rownames(matrix1),rownames(matrix2))
  for (i in 1:length(intersect.gene.name)){
    #tmp.matrix<-rbind(as.vector(intersect.data.dia.second.mean[i,]),as.vector(intersect.ddata.TMT.mean[i,]),stringsAsFactors=F) 
    tmp.matrix<-rbind(matrix1[intersect.gene.name[i],],matrix2[intersect.gene.name[i],]) 
    if (sum(is.na(tmp.matrix))<2 ){
      #print(intersect.gene.name[i])
      print(i)
      tmp.matrix<-t( na.omit(t(tmp.matrix)) )
      tmp.cor<-cor.test(x=as.numeric(as.character(as.vector(tmp.matrix[1,]))),y=as.numeric(as.character(as.vector(tmp.matrix[2,]))),method="pearson")
      tmp.cor.value<- tmp.cor$estimate
      tmp.cor.pvalue<- tmp.cor$p.value
      names(tmp.cor.value)<- intersect.gene.name[i]
      names(tmp.cor.pvalue)<- intersect.gene.name[i]
      tmp.c<-c(tmp.c,tmp.cor.value)
      tmp.c.pvalue<-c(tmp.c.pvalue,tmp.cor.pvalue)
    }
  }
  dia.second.and.TMT.gene.cor<- as.data.frame(cbind(tmp.c,tmp.c.pvalue))
  outname
  pic.name<- paste(outname,".png",sep = "")
  p<-hist(dia.second.and.TMT.gene.cor[,1])
  png(pic.name)
  print(p)
  dev.off()
  csv.name<- paste(outname,".csv",sep = "")
  colnames(dia.second.and.TMT.gene.cor)<-c("cor","pvalue")
  write.csv(dia.second.and.TMT.gene.cor,csv.name)
  
} 

#matrix1<-data.TMT.mean
##matrix2<-data.dia.second.mean
#outname<-"test.3.10"



gene.cor.cal.from.two.matix(data.TMT.mean,data.dia.second.mean,"TMTmRNA4mean")

#dev.off()

