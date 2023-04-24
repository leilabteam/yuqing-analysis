#install.packages("psych")
library(psych)
data <- as.data.frame(read.csv("mRNAmatrix.csv",row.names = 1))

dim(data)
data<- t(data)

#group<- c(rep("A",3),rep("B",3),rep("C",3),rep("D",3),rep("E",3),rep("F",3),rep("G",3),rep("H",3),rep("I",3),rep("J",3),rep("K",3)) #分组信息根据课题设定
#快速分组方法
group_num<- nrow(data)/4
group<-paste("A",c(rep(1:group_num,each=4)),sep="")
rownames(data)<- group 
dim(data)

#开始进行anova分析：
anova_matrix<- data
#anova_matrix[is.na(anova_matrix)]<- 1 

mode(anova_matrix)
####测试某个基因，同时收集行名
anova_matrix[,2]
fit <- aov(anova_matrix[,2]~group)
A_test<-summary(fit)
Anova_Pvalue <-A_test[[1]][1,5] 
multi_p<-TukeyHSD(fit)
p_data<-as.matrix(multi_p$group[,4])
#如果表达矩阵中NA值过多，则row_names不能加入p_data名称，同时下面的Anova_test函数也不能加入p_data的数据，因为每个基因的p_data的数量不一样，导致最后无法data.frame
#row_names<- c("Anova_Pvalue",rownames(p_data),paste(rownames(p_data),"fold"))
row_names<- c("Anova_Pvalue")
group_names<-rownames(p_data)

#建立方差分析函数，同时用循环得到每两组之间的差异倍数:
#针对蛋白数据做调整，只输出p值
Anova_test <- function(x){
  fit <- aov(x~group)
  A_test<-summary(fit)
  Anova_Pvalue <-A_test[[1]][1,5] 
  multi_p<-TukeyHSD(fit)
  p_data<-as.matrix(multi_p$group[,4])
  Test_summary<- rbind(Anova_Pvalue,p_data)
  
  des_data<-describeBy(x, list(group))
  length_des<- length(des_data)
  
  #开始循环计算foldchange值
  #  mean_data<- c(mean_data,as.data.frame(des_data[j])[3]/as.data.frame(des_data[i])[3])
  #mean_data<-t(as.data.frame(mean_data))
  #Test_summary2 <- rbind(Anova_Pvalue,p_data, mean_data) #如果表达矩阵中NA值过多，则不能加入p_data
  Test_summary2 <<- Anova_Pvalue
  #  Test_summary2 <<- rbind(Anova_Pvalue)
  #Test_summary2[is.na(Test_summary2)] <- "N"
  count<<- count+1
  print(count)
  return(Test_summary2)
}

#check 某些基因的NA值是否过多，如果过多则需要去掉这些基因，比如一个基因只在一个组内有表达，则无法求组件差异，就会报错退出：
#比如有66个样本，3个样本一组，这里设置为60个NA以上才会去掉，保证至少两个组之间可以求差异不报错；
#此处去掉部分NA值数量多的基因；
ve<-c()
for  ( i in 1:ncol(anova_matrix))  {
  if ( sum(is.na(anova_matrix[,i])) >= nrow(anova_matrix)-6 ) {
    ve<-c(ve,i)
  }
}
anova_matrix<-anova_matrix[,-ve]
mode(anova_matrix)
#使用apply函数循环做方差检验
count<-0
#anova_matrix<-anova_matrix[,c(10:200)]
#apply(anova_matrix, 2, Anova_test)
rm(test_res)
test_res <- apply(anova_matrix, 2, Anova_test)
test_res
#做data.frame可能还是会有报错
test_res2 <- as.data.frame(test_res,check.rows=F)
colnames(test_res2) <- row_names

#FDR矫正
#p.adjust(p, method = p.adjust.methods, n = length(p))
dim(test_res)
padjust<- as.data.frame(p.adjust(test_res2[,1],method = p.adjust.methods,n=length(test_res2[,1])))
test_res_padjust<- cbind(padjust,test_res2)

#输出最终文件
write.csv(test_res_padjust, file="20200411-anova_test.csv")
