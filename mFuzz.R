data<- read.csv("spectronautmatrix.csv",row.names = 1)

#每四个重复取平局值：
mean_data<- data.frame(c(1:nrow(data)))
for( i in seq(1,length(data),by=4) ) {
  ii<- i+3
  s_data<-data[,c(i:ii)]
  cc <- data.frame(apply(s_data,1,mean,na.rm=TRUE))
  
  mean_data<-cbind(mean_data,cc)
}
mean_data<- mean_data[,-1]




#挑选出需要的gene list
use_gene<- read.csv("pvalue0.05.csv",header = F)
use_mean_gene<- mean_data[match(as.matrix(use_gene),rownames(mean_data)),]





#colnames(mean_data)<- c("0h","1h","3h","6h","9h","12h","1d","1.5d","2d","2.5d","3d","4d","5d","6d","7d","10d","14d")
#scale 数据集
#scale_data<- scale(t(mean_data), center = TRUE, scale = TRUE)
#scale_data<-t(scale_data)

############################################
#2.进行mfuzz分析：
library(Mfuzz)
count <- data.matrix(use_mean_gene)
count<- na.omit(count)
dim(count)
eset <- new("ExpressionSet",exprs = count)
# 根据标准差去除样本间差异太小的基因
#这里暂时不去除
eset <- filter.std(eset,min.std=0)

#标准化
eset <- standardise(eset)

#聚类
# 聚类个数
c <- 16
#  评估出最佳的m值
m <- mestimate(eset)


dim(eset)
# 聚类
cl <- mfuzz(eset, c = c, m = m)

# 查看每个cluster中的基因个数
cl$size
# 提取某个cluster下的基因
cl$cluster[cl$cluster == 1]
# 查看基因和cluster之间的membership
cl$membership
gene_membership<- as.data.frame(cl$membership)
write.csv(gene_membership,"p0.01NA-mean-gene_membership.csv")

gene_cluster<-as.data.frame(cl$cluster)
write.csv(gene_cluster,"p0.01NA-mean-gene_cluster.csv")

#可视化：
mfuzz.plot(
  eset,
  cl,
  mfrow=c(3,3),
  new.window= FALSE)


png("p0.01NA-mean-16-cluster.png",width = 1600,height = 1600, res=72*3)
mfuzz.plot(
  eset,
  cl,
  mfrow=c(2,4),
  new.window= FALSE)
dev.off()

#保存相关数据
data1<-t(as.data.frame(eset))
write.csv(data1,"use_mean_12time_data_log_normlization.csv")

