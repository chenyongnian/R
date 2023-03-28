library(GEOquery)
library(limma)
gse <- getGEO("GSE75214", GSEMatrix = TRUE)
expr <- exprs(gse[[1]]);
expr<-data.frame(rownames(expr),expr);
names(expr)[1]<-"PROBEID"
pheno <- pData(gse[[1]])

#probe ID to gene ID R
library(hugene10sttranscriptcluster.db)
ann_table <- select(hugene10sttranscriptcluster.db, keys = rownames(expr), columns = "SYMBOL")
probID<-na.omit(ann_table)
expr<-merge(probID,expr,by="PROBEID")
#合并同一基因的多个探针：求平均
EXPR<-aggregate(expr[,-c(1,2)],list(expr$SYMBOL),"mean")

#患者分组：UC,CD
CD_GSM<-rownames(pheno)[grep("CD",pheno$title)]
CD<-EXPR[,c(1,which(names(SUB) %in% CD_GSM))]
UC_GSM<-rownames(pheno)[grep("UC",pheno$title)]
UC<-EXPR[,c(1,which(names(SUB) %in% UC_GSM))]

#result<-cor.test(as.numeric(EXPR[which(EXPR$Group.1=="FABP5"),-1]),as.numeric(EXPR[which(EXPR$Group.1=="TNF"),-1]),method='pearson')
#把result做成数据框，每次合并到外部的结果数据框中。使用rbind（）
p_value<-result$p.value
cor<-as.numeric(result$estimate)
result<-data.frame(p_value,cor)
relevanceAnalysis<-function(gene1,gene2){
  gene2Name<-gene2$Group.1
  result<-cor.test(as.numeric(gene1[,-1]),as.numeric(gene2[,-1]),method='pearson')
  p_value<-result$p.value
  cor<-as.numeric(result$estimate)
  result<-data.frame(gene2Name,p_value,cor)
  return(result)
}
FABP5_CD<-CD[which(EXPR$Group.1=="FABP5"),]
FABP5_UC<-UC[which(EXPR$Group.1=="FABP5"),]
#用于存储各基因与FABP5相关性的数据框
cor_CD<-data.frame()
cor_UC<-data.frame()
for (name in CD$Group.1) {
  gene2<-CD[which(CD$Group.1==name),]
  result<-relevanceAnalysis(FABP5_CD,gene2)
  cor_CD<-rbind(cor_CD,result)
}

for (name in UC$Group.1) {
  gene2<-UC[which(UC$Group.1==name),]
  result<-relevanceAnalysis(FABP5_UC,gene2)
  cor_UC<-rbind(cor_UC,result)
}
write.csv(cor_CD,"cor_CD.csv")
write.csv(cor_UC,"cor_UC.csv")
