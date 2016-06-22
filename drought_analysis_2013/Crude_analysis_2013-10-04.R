setwd("/Users/gosia/Analysis/Shimizu_RNA_seq/")

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
#> all( samps$sample_name == colnames(x) )
#[1] TRUE

# time format change 

tmO <- strptime(samps$year.month.day, "%Y.%m.%d")
tm <- as.numeric(tmO)

full.time <- read.table("Data/Unique_days_short2.csv")
full.time <- strptime(full.time[1:374,], "%d.%m.%y")
full.time <- as.numeric(full.time)

full.time <- data.frame(full.time, as.character(colorRampPalette(c("red","blue"))(374)))

colors.mds.time <- merge(data.frame(samps$sample_name, tm), full.time, by.x=2, by.y=1, all.x=T)
rownames(colors.mds.time) <- colors.mds.time[,2]

plot(1:length(colors.mds.time[,1]),colors.mds.time[,1], col=as.character(colors.mds.time[,3]))

new.samps <- read.table("Samples_out/new_samps_interpolation_November.csv", header=TRUE, sep=";", stringsAsFactors=FALSE)

colors.palette <- as.character(colorRampPalette(c("red","blue"))(50))
colors.wp <- new.samps[, c("sample_name", "Water.Potential")]
colors.wp$Water.Potential <- 100*round(colors.wp$Water.Potential, digits=2)
colors.wp$colors <- colors.palette[colors.wp$Water.Potential - min(colors.wp$Water.Potential[!is.na(colors.wp$Water.Potential)]) +1 ]
colors.wp$colors[is.na(colors.wp$colors)] <- colors.palette[50]
rownames(colors.wp) <- colors.wp$sample_name

colors.trees <- as.factor(samps$tree_ID)
levels(colors.trees) <- 1:length(levels(colors.trees))
levels(colors.trees) <- apply(col2rgb(levels(colors.trees)), 2, function(c2r){ rgb(c2r[1], c2r[2], c2r[3], maxColorValue=255)})
colors.trees <- as.character(colors.trees)
colors.trees[samps$developmental_stage=="flower_bud"] <- "orange"
colors.trees <- data.frame(as.character(samps$sample_name), as.character(samps$tree_ID), colors.trees)
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

blast <- read.table("Data/new_data_18Oct/trinity_tair10_obh_clean.csv",sep=",")
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

#all( samps$sample_name == colnames(d) )
#[1] TRUE

# plots: MDS, 
library(DESeq)
library(ggplot2)
library(ggdendro)

cols <- c("black","red","green","blue","cyan","magenta")


pdf("Analysis_till_Nov06/Crude_analysis/mds3_500.pdf",w=10,h=10)
mds <- plotMDS(d, col=as.character(colors.trees[rownames(d$samples), 3]), top=500, main="method=logFC, prior.count=2")
#mds <- plotMDS(d, col=as.character(colors.mds.time[rownames(d$samples),3]), top=500, main="method=logFC, prior.count=2")
#mds <- plotMDS(d, col=as.character(colors.wp[rownames(d$samples),3]), top=500, main="method=logFC, prior.count=2")
mds10 <- plotMDS(d,col=as.character(colors.trees[rownames(d$samples), 3]),main="method=logFC, prior.count=10", prior.count=10, top=500)
#mds10 <- plotMDS(d,col=as.character(colors.mds.time[rownames(d$samples),3]),main="method=logFC, prior.count=10", prior.count=10, top=500)
mdsb <- plotMDS(d, col=as.character(colors.trees[rownames(d$samples), 3]), method="bcv", top=500)
#mdsb <- plotMDS(d, col=as.character(colors.mds.time[rownames(d$samples),3]), method="bcv", top=500)

# hc <- hclust(as.dist(mdsb$distance.matrix))
# #plot(hc, col=as.numeric(d$samples$group), hang=-1)
# ggdendrogram(hc, rotate=TRUE, size=4, theme_dendro=FALSE, 
#              color=(cols[as.numeric(d$samples$group)])[hc$order])
# 
# cds <- newCountDataSet(countData=d$counts, conditions=samps$tree_ID)
# cds <- estimateSizeFactors( cds )
# cds <- estimateDispersions( cds, method="blind" )
# vsd <- varianceStabilizingTransformation( cds )
# plotPCA( vsd, ntop=10000)
# plotPCA( vsd, ntop=5000)
# plotPCA( vsd, ntop=1000)
dev.off()

pdf("Analysis_till_Nov06/Crude_analysis/mds3_1000.pdf",w=10,h=10)
mds <- plotMDS(d, col=as.character(colors.trees[rownames(d$samples), 3]), top=1000, main="method=logFC, prior.count=2")
mds10 <- plotMDS(d,col=as.character(colors.trees[rownames(d$samples), 3]),main="method=logFC, prior.count=10", prior.count=10, top=1000)
mdsb <- plotMDS(d, col=as.character(colors.trees[rownames(d$samples), 3]), method="bcv", top=1000)
dev.off()


#p <- function(...) smoothScatter(..., nrpoints=0, nbin=64, add=TRUE, transformation = function(x) x^.1)
#g <- unique(d$samples$group)
#pdf("all.pdf",w=20,h=20)
#for(i in g)
#pairs( log2(1+d$counts[,as.character(d$samples$group)==i]), pch=".", panel=p, lower.panel=NULL)
#dev.off()



dbak <- d

# take subset
keep <- samps$tree_ID != "990" & samps$developmental_stage != "flower_bud"
d <- d[,keep]

# design model matrix

design <- model.matrix(~-1 + tree_ID, data=samps[keep,])
design <- design[, colSums(design)>0]

# estimate dispersion

d <- estimateGLMCommonDisp(d,design)
d <- estimateGLMTrendedDisp(d,design)
d <- estimateGLMTagwiseDisp(d,design)

# glmFit fits genewise negative binomial glms, all with the same design matrix but possibly different dispersions, offsets and weights

fit <- glmFit(d,design)

# make contrast to test any difference from control
mc <- makeContrasts(contrasts=c("(tree_ID970+tree_ID8212)/2-(tree_ID1099+tree_ID1377)/2", "tree_ID8266-(tree_ID1099+tree_ID1377)/2", "tree_ID8266-(tree_ID1099+tree_ID1377+tree_ID970+tree_ID8212)/4"),levels=colnames(design))
colnames(mc) <- c("sheet_vs_control","flower_vs_control","flower_vs_rest")

# glmLRT conducts likelihood ratio tests for one or more coefficients in the linear model

lrt <- glmLRT(fit,contrast=mc[,"flower_vs_rest"])



plotGene <- function(d,id="GID000013_4400", tm, cols=c("black","blue","blue","blue","orange","black"), 
                     FUN=sqrt,main=NULL,...) 
{
  cps <- cpm(d, normalized.lib.sizes=TRUE)
  x <- cps[id,]
  g <- d$samples$group
  ug <- unique(as.character(g))
  # main title for plot
  if (is.null(main)) {
    main <- id
    m <- match(id, d$genes$assembl_id)
    if( !is.na(d$genes$at_id[m]) ) main <- paste0(main," // ",d$genes$at_id[m])
    if( nchar(d$genes$at_symbol[m]) > 0 ) main <- paste0(main," // ",d$genes$at_symbol[m])
  }
  for(tr in 1:length(ug)) {
    k <- g==ug[tr]
    if(tr==1) {
      plot(tm[k], FUN(x[k]), pch=19, lwd=3, type="b", ylim=range(sqrt(x)),col=cols[tr], 
           xlab="Date", ylab="expression (counts per million)",main=main,...)
    } else {
      points(tm[k], FUN(x[k]), pch=19, lwd=3, type="b",col=cols[tr])
    }
  }
  
  dt <- as.numeric(as.POSIXlt("2009-01-01"))
  abline(v=dt,lwd=4,col="grey")
  abline(v=as.numeric(as.POSIXlt(paste0("2009-",1:12,"-01"))),lwd=1,col="grey")
  abline(v=as.numeric(as.POSIXlt(paste0("2008-",11:12,"-01"))),lwd=1,col="grey")
  
  legend("topright",as.character(ug),col=cols,pch=19,lwd=3,bg="white")
}



rn <- topTags(lrt,n=100)$table$assembl_id # top genes

pdf("Crude_analysis/differential_flower_vs_rest.pdf",w=18,h=12)
par(mfrow=c(2,3))
for(ii in rn)
  plotGene(d,ii,tm[keep],cols=cols)
dev.off()


# tt <- topTags(lrt,n=nrow(d))$table$assembl_id
# tt <- tt[(nrow(d)-100) : nrow(d)]
# 
# pdf("Crude_analysis/NOTdifferential_drought_vs_control.pdf",w=18,h=12)
# par(mfrow=c(2,3))
# for(ii in tt)
#   plotGene(d,ii,tm[keep],cols=c("grey30","darkblue","blue","grey50","orange"))
# dev.off()


# tt <- topTags(lrt,n=nrow(d))$table
# write.table(tt,"crude_global_analysis_12june2013.xls", row.names=FALSE, sep="\t")

