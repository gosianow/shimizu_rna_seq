library(DESeq2)


dds <- makeExampleDESeqDataSet(n=100,m=12)
dds$genotype <- factor(rep(rep(c("I","II"),each=3),2))
design(dds) <- ~ genotype + condition + genotype:condition
dds <- DESeq(dds)
resultsNames(dds)

model.matrix(~ genotype + condition + genotype:condition, data=as.data.frame(colData(dds)))

res <- results(dds, contrast=c("condition","B","A"))


## Example 3: two conditions, three genotypes
# ~~~ Using interaction terms ~~~


dds <- makeExampleDESeqDataSet(n=100,m=18)
dds$genotype <- factor(rep(rep(c("I","II","III"),each=3),2))
design(dds) <- ~ genotype + condition + genotype:condition
dds <- DESeq(dds)
resultsNames(dds)

model.matrix(~ genotype + condition + genotype:condition, data=as.data.frame(colData(dds)))


################
# Time course experiments
################

library("ggplot2")

library("fission")
data("fission")



ddsTC <- DESeqDataSet(fission, ~ strain + minute + strain:minute)

coldata  <- data.frame(colData(ddsTC))

modelMatrix <- model.matrix(~ strain + minute + strain:minute, coldata)

qr(modelMatrix)$rank
ncol(modelMatrix)



ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ strain + minute)
resTC <- results(ddsTC)
resTC$symbol <- mcols(ddsTC)$symbol
head(resTC[order(resTC$padj),],4)


data <- plotCounts(ddsTC, which.min(resTC$padj), intgroup=c("minute","strain"), returnData=TRUE)


ggp <- ggplot(data, aes(x=minute, y=count, color=strain, group=strain)) + 
  geom_point() + 
stat_smooth(se=FALSE,method="loess") +  
scale_y_log10()

pdf("deseq2_ts1.pdf")
print(ggp)
dev.off()


resultsNames(ddsTC)



### Wald tests for the log2 fold changes at individual time points can be investigated using the test argument to results:


res30 <- results(ddsTC, name="strainmut.minute30", test="Wald")
res30[which.min(resTC$padj),]

### Clustering

betas <- coef(ddsTC)
colnames(betas)


library("pheatmap")


topGenes <- head(order(resTC$padj),20)
mat <- betas[topGenes, -c(1,2)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr

pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101), cluster_col=FALSE, filename = "deseq2_hm.pdf")







### Switch to minute:strain


ddsTC <- DESeqDataSet(fission, ~ strain + minute + minute:strain)

ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ strain + minute)
resTC <- results(ddsTC)
resTC$symbol <- mcols(ddsTC)$symbol
head(resTC[order(resTC$padj),],4)


data <- plotCounts(ddsTC, which.min(resTC$padj), intgroup=c("minute","strain"), returnData=TRUE)


ggp <- ggplot(data, aes(x=minute, y=count, color=strain, group=strain)) + 
  geom_point() + 
stat_smooth(se=FALSE,method="loess") +  
scale_y_log10()

pdf("deseq2_ts2.pdf")
print(ggp)
dev.off()



















