setwd("/home/gosia/Analysis/Shimizu_RNA_seq/")

# parse sample information
samps <- read.table("Data/orig/sample_list.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE)

#write.table(samps, "sample_list.xls", sep="\t", row.names = F)

# get table of counts
ch <- as.character(read.table("Data/orig/raw_data.csv", sep=",",nrow=1,stringsAsFactors=FALSE)[1,-c(1,2,57)])
x <- read.table("Data/orig/raw_data.csv", sep=",",skip=3, header=F, row.names=1)
nc <- ncol(x)
x <- x[,-c(1, nc)]
stopifnot( ncol(x)==length(ch) )
colnames(x) <- ch

head(x)
names(x)


# time format change 
tmO <- strptime(samps$year.month.day, "%Y.%m.%d")
tm <- as.numeric(tmO)

library(stringr)
# colors by tree
colors.trees <- as.factor(samps$tree_ID)
levels(colors.trees) <- c("orange", "cyan", "green", "blue", "red", "magenta")
colors.trees <- as.character(colors.trees)
colors.trees[samps$developmental_stage=="flower_bud"] <- "darkmagenta"
colors.trees <- data.frame(as.character(samps$sample_name), as.character(samps$tree_ID), colors.trees)
colors.trees$samps.short.date <- paste(str_sub(string= colors.trees[,1], start = -6L, end = -5L), str_sub(string= colors.trees[,1], start = -4L, end = -3L), sep=".")
rownames(colors.trees) <- samps$sample_name


library(edgeR)
d <- DGEList(x, group=samps$tree_ID)
d <- calcNormFactors(d)

samps$flowered <- samps$tree_ID=="8266"
samps$tree_ID <- as.factor(samps$tree_ID)

# make sure a gene is expressed (CPM > 1) in more than 2 samples
cps <- cpm(d, normalized.lib.sizes=TRUE)
d <- d[ rowSums(cps>1) > 2, ]

# reorder by tree and time
o <- order(d$samples$group, tm)
samps <- samps[o,]
d <- d[,o]
tmO <- tmO[o]
tm <- tm[o]


# create d$genes -> TAIR annotation / AT ID / functional description

blast <- read.table("Data/genes_descr_control/best_hit_blast_result_Sl_predicted_exons.txt",sep=",")
m <- match( rownames(d$counts), blast$V1)
d$genes <- data.frame(assembl_id=rownames(d$counts), at_id=as.character(blast$V2[m]), at_symbol="", stringsAsFactors=FALSE)


library(org.At.tair.db)
nna <- !is.na(d$genes$at_id)
syms <- mget(d$genes$at_id[nna], org.At.tairSYMBOL)
syms <- sapply(syms, function(u) if(is.na(u[1])) "" else paste(u,collapse=";"))
d$genes$at_symbol[nna] <- syms

tfd <- read.table("Data/genes_descr_control/TAIR10_functional_descriptions",sep="\t", header=TRUE, quote="", comment.char="", stringsAsFactors=FALSE)
t10id <- gsub("\\.[1-9]","",tfd$Model_name)
m <- match(d$genes$at_id[nna],t10id)
d$genes$description[nna] <- as.character(tfd$Short_description[m])

cps <- cpm(d, normalized.lib.sizes=TRUE)
d$genes <- cbind(d$genes, round(cps,1))


# plots: MDS, 
library(DESeq)
library(ggplot2)

dir.create("Plots_Crude_analysis")


pdf("Plots_Crude_analysis/mds_500.pdf",w=7,h=7)
mds <- plotMDS(d, col=as.character(colors.trees[rownames(d$samples), 3]), top=500, main="method=logFC, prior.count=2", labels=colors.trees[rownames(d$samples), "samps.short.date"])
mds10 <- plotMDS(d,col=as.character(colors.trees[rownames(d$samples), 3]),main="method=logFC, prior.count=10", prior.count=10, top=500, labels=colors.trees[rownames(d$samples), "samps.short.date"])
mdsb <- plotMDS(d, col=as.character(colors.trees[rownames(d$samples), 3]), method="bcv", top=500, labels=colors.trees[rownames(d$samples), "samps.short.date"])
dev.off()


pdf("Plots_Crude_analysis/mds_1000.pdf",w=7,h=7)
mds <- plotMDS(d, col=as.character(colors.trees[rownames(d$samples), 3]), top=1000, main="method=logFC, prior.count=2", labels=colors.trees[rownames(d$samples), "samps.short.date"])
mds10 <- plotMDS(d,col=as.character(colors.trees[rownames(d$samples), 3]),main="method=logFC, prior.count=10", prior.count=10, top=1000, labels=colors.trees[rownames(d$samples), "samps.short.date"])
mdsb <- plotMDS(d, col=as.character(colors.trees[rownames(d$samples), 3]), method="bcv", top=1000, labels=colors.trees[rownames(d$samples), "samps.short.date"])
dev.off()


