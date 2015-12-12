source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
setwd("C:/Users/malmwamu/Documents/Epid600/Epid600-project/")
directory <- "C:/Users/malmwamu/Documents/Epid600/Epid600-project/Normal-AMD-PR/"
sampleFiles <- grep("HTseq-reverse", list.files(directory), value = TRUE)
sampleCondition <- c("NPR", "NPR", "NPR", "NPR", "AMDPR", "AMDPR", "AMDPR","AMDPR", "AMDPR", "NPR", "NPR", "NPR")
sampleTable <- data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)
nrow(ddsHTSeq)
ddsHTSeq<- ddsHTSeq[rowSums(counts(ddsHTSeq)) > 1, ]
nrow(ddsHTSeq)
ddsHTSeq <-DESeq(ddsHTSeq)
res <- results(ddsHTSeq)
res
resOrdered <- res[order(res$padj),]
resOrdered
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)
res05 <- results(ddsHTSeq, alpha = 0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm = TRUE)
plotMA(res, main="Peripheral Retina Normal versus AMD eyes DESeq2")
abline(h = c(-2, 2), col='green')
abline(h = c(-1, 1), col='blue')
write.csv(as.data.frame(resOrdered), file = "Results_Retina-Normal-vs-AMD_DESEq2.csv")
Results.unpaired <- as.data.frame(resOrdered)
resSig <- subset(resOrdered, padj < 0.05)
resSig
write.csv(as.data.frame(resSig), file = "DE_Retina-Normal-vs-AMD_DESeq2.csv")
DE.unpaired <- as.data.frame(resSig)
#volcano plot
Results.unpaired$threshold = as.factor(abs(Results.unpaired$log2FoldChange) > 1 & Results.unpaired$padj < 0.05)

g <-ggplot(data = Results.unpaired, aes(log2FoldChange, -log10(padj), colour=threshold)) + geom_point(alpha=0.7, size=3) + coord_cartesian(xlim = c(-5,5), ylim = c(0,50)) + xlab("log2 fold change") + ylab("-log10 FDR adjusted p-value")+ theme(legend.position = "none")+ ggtitle("Volcano Plot Normal vs AMD Retina \n DESeq2") + theme(plot.title=element_text(vjust = -4.5, face = "italic", size = 24), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(colour= "blue", size = 16), axis.text.y = element_text(colour= "blue",size = 16))
g

SigQ <- Results.unpaired[Results.unpaired$padj < 0.05 & abs(Results.unpaired$log2FoldChange) > 1,]          
#SigQ= 2878