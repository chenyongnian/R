library(GEOquery)
library(limma)
gse <- getGEO("GSE75214", GSEMatrix = TRUE)
expr <- exprs(gse[[1]]);
expr<-data.frame(rownames(expr),expr);
names(expr)[1]<-"PROBEID"
pheno <- pData(gse[[1]])

#probe ID to gene ID R
library(hugene10sttranscriptcluster.db)
ann_table <- select(hugene10sttranscriptcluster.db, keys = rownames(expr), columns = c("SYMBOL", "ENTREZID"))

# Select genes of interest
genes <- c("IL6", "IL1B", "TNF", "FABP5")
probID<-ann_table[ann_table$SYMBOL %in% genes,]
probID<-probID[,-3]
expr_sub<- expr[rownames(expr) %in% probID$PROBEID, ]
sub<-merge(probID,expr_sub,by="PROBEID")
#合并同一基因的多个探针：求平均
SUB<-aggregate(sub[,-c(1,2)],list(sub$SYMBOL),"mean")

#患者分组：UC,CD
CD_GSM<-rownames(pheno)[grep("CD",pheno$title)]
CD<-SUB[,c(1,which(names(SUB) %in% CD_GSM))]
UC_GSM<-rownames(pheno)[grep("UC",pheno$title)]
UC<-SUB[,c(1,which(names(SUB) %in% UC_GSM))]
#export dataframe 
write.table(CD,"CD.csv", sep = ",", quote = FALSE, row.names = FALSE)
write.table(UC,"UC.csv", sep = ",", quote = FALSE, row.names = FALSE)
