#####################################################################################################
# BioC 3.0

### crude analysis - MDS plots
### plots of expression for flowering genes

#####################################################################################################

### load data

# RPath <- "/home/gosia/R/R_Shimizu_RNA_seq/Shimizu_analysis_2014-12-03/"
# dataPath <- "/home/Shared/data/seq/Shimizu_RNA_seq/Data/"
# analysisPath <- "Analysis_2014-12-03"
# analysisPath <- paste0("/home/Shared/data/seq/Shimizu_RNA_seq/", analysisPath)

RPath <- "/Users/gosia/Dropbox/Shimizu_time_course_RNA_seq/R_Shimizu_RNA_seq/Shimizu_analysis_2014-12-03/"
dataPath <- "/Users/gosia/Dropbox/Shimizu_time_course_RNA_seq/Shimizu_RNA_seq/Data/"
analysisPath <- "Analysis_2014-12-03"
analysisPath <- paste0("/Users/gosia/Dropbox/Shimizu_time_course_RNA_seq/Shimizu_RNA_seq/", analysisPath)

dir.create(analysisPath, showWarnings = FALSE)
setwd(analysisPath)

load("Shimizu_workspace.Rdata")


all(colnames(x)==rownames(new.samps))

#####################################################################################################
### MDS plots
#####################################################################################################

dim(x)


library(edgeR)
d <- DGEList(x, group=new.samps$tree_ID)
d <- calcNormFactors(d)

# make sure a gene is expressed (CPM > 1) in more samples
cps <- cpm(d, normalized.lib.sizes=TRUE)
d <- d[ rowSums( cps > 1 ) > 27, ]
dim(d$counts)


# plots: MDS, 
library(ggplot2)
library(EDASeq)
library(ggdendro)

dir.create("Plots_MDS", showWarnings = FALSE)


pdf("Plots_MDS/mds_500.pdf",w=10,h=10)
# mds <- plotMDS(d, col=new.samps$tree_col, top=500, main="method=logFC, prior.count=2", labels=new.samps$short.name)
# plotMDS(mds, col=new.samps$colors.time, main="method=logFC, prior.count=2", labels=new.samps$short.name) 
# mds10 <- plotMDS(d,col=new.samps$tree_col, main="method=logFC, prior.count=10", prior.count=10, top=500, labels=new.samps$short.name)
mdsb <- plotMDS(d, col=new.samps$tree_col, method="bcv", top=500, labels=new.samps$short.name)
plotMDS(mdsb, col=new.samps$tree_col, labels=new.samps$short.name, cex.axis = 2, xlim=c(-1, 1), ylim=c(-1, 1), cex = 1.2)
dev.off()



pdf("Plots_MDS/dendrogram_500.pdf", w=10, h=10)
hc <- hclust(as.dist(mdsb$distance.matrix))
ggdendrogram(hc, size=4, theme_dendro = FALSE, rotate=TRUE) + theme(axis.text.y = element_text(colour = new.samps[hc$order,"tree_col"]))
ggdendrogram(hc, size=4, theme_dendro = FALSE, rotate=TRUE) + theme(axis.text.y = element_text(colour = new.samps[hc$order,"colors.time"]))
dev.off()


pdf("Plots_MDS/mds_1000.pdf",w=10,h=10)
mds <- plotMDS(d, col=new.samps$tree_col, top=1000, main="method=logFC, prior.count=2", labels=new.samps$short.name)
mds10 <- plotMDS(d,col=new.samps$tree_col,main="method=logFC, prior.count=10", prior.count=10, top=1000, labels=new.samps$short.name)
mdsb <- plotMDS(d, col=new.samps$tree_col, method="bcv", top=1000, labels=new.samps$short.name)
dev.off()


pdf("Plots_MDS/dendrogram_1000.pdf", w=10, h=10)
hc <- hclust(as.dist(mdsb$distance.matrix))
ggdendrogram(hc, size=4, theme_dendro = FALSE, rotate=TRUE) + theme(axis.text.y = element_text(colour = new.samps[hc$order,"tree_col"]))
dev.off()



pdf("Plots_MDS/pca.pdf",w=10,h=10)
plotPCA(d$counts, col = new.samps$tree_col)
plotPCA(cpm(d, normalized.lib.sizes=TRUE), col = new.samps$tree_col)
dev.off()


pdf("Plots_MDS/rle.pdf",w=15,h=7)
plotRLE(d$counts, col = new.samps$tree_col, outline = FALSE, main="Raw counts")
plotRLE(cpm(d, normalized.lib.sizes=TRUE), col = new.samps$tree_col, outline = FALSE, main = "CPM normalization")
dev.off()



######### MDS plots of housekeeping genes 

Housekeeping_Genes <- read.table(paste0(dataPath, "GeneControlSets/Housekeeping_Genes/Table4_Arabidopsis_Housekeepieng_Gene_ID.txt"), head=F, skip=1, stringsAsFactors = FALSE)
head(Housekeeping_Genes)
Housekeeping_Genes <- toupper(Housekeeping_Genes[,1])
length(Housekeeping_Genes)

write.table(Housekeeping_Genes, paste0(dataPath, "GeneControlSets/Housekeeping_Genes/Table4_Arabidopsis_Housekeepieng_Gene_ID_upper.txt"), quote=FALSE, col.names=FALSE, row.names=FALSE)


library(edgeR)
d <- DGEList(x, group=new.samps$tree_ID)
d <- calcNormFactors(d)

# make sure a gene is expressed (CPM > 1) in more than 2 samples
cps <- cpm(d, normalized.lib.sizes=TRUE)
d <- d[ rowSums( cps > 1 ) > 27, ]
dim(d$counts)

d <- d[rownames(d$counts) %in% genes.full.description[genes.full.description$AT_ID %in% Housekeeping_Genes, "ID"], ]
dim(d$counts)



pdf("Plots_MDS/Housekeeping_mds_500.pdf",w=10,h=10)

mdsb <- plotMDS(d, col=new.samps$tree_col, method="bcv", top=500, labels=new.samps$short.name, cex.axis = 2, xlim=c(-1, 1), ylim=c(-1, 1), cex = 1.2)

dev.off()





pdf("Plots_MDS/Housekeeping_mds_500.pdf",w=10,h=10)
mds <- plotMDS(d, col=new.samps$tree_col, top=500, main="method=logFC, prior.count=2", labels=new.samps$short.name)
mds10 <- plotMDS(d,col=new.samps$tree_col, main="method=logFC, prior.count=10", prior.count=10, top=500, labels=new.samps$short.name)
mdsb <- plotMDS(d, col=new.samps$tree_col, method="bcv", top=500, labels=new.samps$short.name)

r <- c(min(mdsb$x, mdsb$y), max(mdsb$x, mdsb$y))
plotMDS(mdsb, col=new.samps$tree_col, labels=new.samps$short.name, xlim=r, ylim=r)
plotMDS(mdsb, col=new.samps$colors.time, labels=new.samps$short.name, xlim=r, ylim=r)

dev.off()



pdf("Plots_MDS/Housekeeping_dendrogram_500.pdf", w=10, h=10)
hc <- hclust(as.dist(mdsb$distance.matrix))
ggdendrogram(hc, size=4, theme_dendro = FALSE, rotate=TRUE) + theme(axis.text.y = element_text(colour = new.samps[hc$order,"tree_col"]))
ggdendrogram(hc, size=4, theme_dendro = FALSE, rotate=TRUE) + theme(axis.text.y = element_text(colour = new.samps[hc$order,"colors.time"]))
dev.off()


pdf("Plots_MDS/Housekeeping_mds_1000.pdf",w=10,h=10)
mds <- plotMDS(d, col=new.samps$tree_col, top=1000, main="method=logFC, prior.count=2", labels=new.samps$short.name)
mds10 <- plotMDS(d,col=new.samps$tree_col,main="method=logFC, prior.count=10", prior.count=10, top=1000, labels=new.samps$short.name)
mdsb <- plotMDS(d, col=new.samps$tree_col, method="bcv", top=1000, labels=new.samps$short.name)
dev.off()


pdf("Plots_MDS/Housekeeping_dendrogram_1000.pdf", w=10, h=10)
hc <- hclust(as.dist(mdsb$distance.matrix))
ggdendrogram(hc, size=4, theme_dendro = FALSE, rotate=TRUE) + theme(axis.text.y = element_text(colour = new.samps[hc$order,"tree_col"]))
dev.off()



pdf("Plots_MDS/Housekeeping_pca.pdf",w=10,h=10)
plotPCA(d$counts, col = new.samps$tree_col)
plotPCA(cpm(d, normalized.lib.sizes=TRUE), col = new.samps$tree_col)
dev.off()


pdf("Plots_MDS/Housekeeping_rle.pdf",w=15,h=7)
plotRLE(d$counts, col = new.samps$tree_col, outline = FALSE, main="Raw counts")
plotRLE(cpm(d, normalized.lib.sizes=TRUE), col = new.samps$tree_col, outline = FALSE, main = "CPM normalization")
dev.off()





#####################################################################################################
### plots of expression for flowering genes from raw data
#####################################################################################################


AT.genes <- unique(genes.full.description[!is.na(genes.full.description$Flowering), c("AT_ID", "AT_symbol")])
AT.id <- genes.full.description[!is.na(genes.full.description$AT_ID), c("ID","AT_ID")]

AT.id <- AT.id[AT.id$ID %in% c("GID031739_3527284", "GID037469_3529277"),]


elim.samps=c(flowered.samps, new.samps$sample_name[!grepl("8266", new.samps$sample_name) ])
trees.order <- trees.order[trees.order$legend %in% c("8266-drought-leaf_bud"), ]


# elim.samps=NULL
x <- x[,!names(x) %in% elim.samps]
new.samps <- new.samps[!new.samps$sample_name %in% elim.samps, ]

all(colnames(x)==rownames(new.samps))


library(edgeR)
d.org <- DGEList(x, group=new.samps$tree_ID)
d.org <- calcNormFactors(d.org)

d.cpm <- cpm(d.org, normalized.lib.sizes=TRUE)
d.cpm.l <- log(d.cpm + min(d.cpm[d.cpm != 0]))

d.cpm[AT.id$ID, "E7_8266_20090416"]



out.dir <- "Plots_of_flowering_genes/"
dir.create(out.dir, showWarnings=F, recursive=T)


pdf(paste0(out.dir , "Flowering_genes_from_S4_table_Nov_fl" ,".pdf"), h=5, w=10)

for(g in 1:nrow(AT.genes)){
  # g=1
  cat(paste(g, ", "))
  genes <- AT.id[AT.id[,2] == AT.genes[g, 1], 1]
  
  if(length(genes!=0)){
    for(j in 1:length(genes)){
      # j=1
 
      plot(0, type="n", main=paste0(AT.genes[g, 1] ,"\n", AT.genes[g, 2] , "\n",  genes[j]) ,xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(d.cpm[genes[j], ])), max(na.omit(d.cpm[genes[j], ]))), xlab="Time", ylab="Gene Expression in cpm", xaxt = "n")
      axis(side=1, at=month.days[,2], labels=month.days[,1])
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[246])
      for(t in trees.order$legend){
        # t=trees.order$legend[1]
        lines(new.samps$time_nr[new.samps$tree_legend == t], d.cpm[genes[j], new.samps$tree_legend == t] , col=trees.order$color[trees.order$legend==t], type="b", pch=trees.order$pch[trees.order$legend==t], cex=trees.order$cex[trees.order$legend==t], lwd=3) 
      }
      legend("topright", legend = trees.order$legend, col=trees.order$color, cex=0.5, text.col=trees.order$color)
      
    }
    
  }
}

dev.off()



pdf(paste0(out.dir , "Variables" ,".pdf"), h=5, w=10)

plot.vars <- c("Water.Potential", "Soil.Moisture")

for(v in plot.vars){
  plot(0, type="n", xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(new.samps[,v])), max(na.omit(new.samps[,v]))), xlab="Time", ylab=v, xaxt = "n")
  for(t in trees.order$legend){
    lines(new.samps$time_nr[new.samps$tree_legend==t], new.samps[new.samps$tree_legend==t, v], col=trees.order$color[trees.order$legend==t], type="b", pch=trees.order$pch[trees.order$legend==t], cex=trees.order$cex[trees.order$legend==t], lwd=4) 
  }
  axis(side=1, at=month.days[,2], labels=month.days[,1])
  legend("bottomright", legend = trees.order$legend, col=trees.order$color, cex=0.7, text.col=trees.order$color)
  
}
dev.off()


### log cpm

pdf(paste0(out.dir , "Flowering_genes_from_S4_table_Nov_fl_log" ,".pdf"), h=5, w=10)

for(g in 1:nrow(AT.genes)){
  # g=1
  cat(paste(g, ", "))
  genes <- AT.id[AT.id[,2] == AT.genes[g, 1], 1]
  
  if(length(genes!=0)){
    for(j in 1:length(genes)){
      # j=1
      
      plot(0, type="n", main=paste0(AT.genes[g, 1] ,"\n", AT.genes[g, 2] , "\n",  genes[j]) ,xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(d.cpm.l[genes[j], ])), max(na.omit(d.cpm.l[genes[j], ]))), xlab="Time", ylab="Gene Expression in cpm", xaxt = "n")
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[246])
      for(t in trees.order$legend){
        # t=trees.order$legend[1]
        lines(new.samps$time_nr[new.samps$tree_legend ==t], d.cpm.l[genes[j], new.samps$tree_legend ==t] , col=trees.order$color[trees.order$legend==t], type="b", pch=trees.order$pch[trees.order$legend==t], cex=trees.order$cex[trees.order$legend==t], lwd=3) 
      }
      axis(side=1, at=month.days[,2], labels=month.days[,1])
      legend("topright", legend = trees.order$legend, col=trees.order$color, cex=0.5, text.col=trees.order$color)
      
    }
    
  }
}


dev.off()




#####################################################################################################
### plots of expression for flowering genes from x.red data 
#####################################################################################################

load(paste0("Plots_MeanForReplicates/" , "data.reduced" ,".RData"))


AT.genes <- unique(genes.full.description[!is.na(genes.full.description$Flowering), c("AT_ID", "AT_symbol")])
head(AT.genes)

AT.genes <- AT.genes[AT.genes$AT_ID %in% c("AT1G65480", "AT2G22540"),]

AT.id <- genes.full.description[!is.na(genes.full.description$AT_ID), c("ID","AT_ID")]
AT.id <- AT.id[AT.id$ID %in% c("GID031739_3527284", "GID037469_3529277"),]

x <- x.red
new.samps <- new.samps.red

all(colnames(x)==rownames(new.samps))


# trees.order <- trees.order[trees.order$legend %in% c("8266-drought-leaf_bud"), ]



library(edgeR)
d.org <- DGEList(x, group=new.samps$tree_ID)
d.org <- calcNormFactors(d.org)

d.cpm <- cpm(d.org, normalized.lib.sizes=TRUE)
# d.cpm.l <- log(d.cpm + min(d.cpm[d.cpm != 0]))
d.cpm.l <- log2(d.cpm + 1)

out.dir <- "Plots_of_flowering_genes/Reduced_data/"
dir.create(out.dir, showWarnings=F, recursive=T)


pdf(paste0(out.dir , "Flowering_genes_from_S4_table" ,".pdf"), h=5, w=10)

for(g in 1:nrow(AT.genes)){
  # g=1
  cat(paste(g, ", "))
  genes <- AT.id[AT.id[,2] == AT.genes[g, 1], 1]
  
  if(length(genes!=0)){
    for(j in 1:length(genes)){
      # j=1
      
      plot(0, type="n", main=paste0(AT.genes[g, 1] ,"\n", AT.genes[g, 2] , "\n",  genes[j]) ,xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(d.cpm[genes[j], ])), max(na.omit(d.cpm[genes[j], ]))), xlab="Time", ylab="Gene Expression in cpm", xaxt = "n")
      axis(side=1, at=month.days[,2], labels=month.days[,1])
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[246])
      for(t in trees.order$legend){
        # t=trees.order$legend[1]
        lines(new.samps$time_nr[new.samps$tree_legend == t], d.cpm[genes[j], new.samps$tree_legend == t] , col=trees.order$color[trees.order$legend==t], type="b", pch=trees.order$pch[trees.order$legend==t], cex=trees.order$cex[trees.order$legend==t], lwd=3) 
      }
      legend("topright", legend = trees.order$legend, col=trees.order$color, cex=0.5, text.col=trees.order$color)
      
    }
    
  }
}

dev.off()



### log cpm

pdf(paste0(out.dir , "Flowering_genes_from_S4_table_log" ,".pdf"), h=5, w=10)

for(g in 1:nrow(AT.genes)){
  # g=1
  cat(paste(g, ", "))
  genes <- AT.id[AT.id[,2] == AT.genes[g, 1], 1]
  
  if(length(genes!=0)){
    for(j in 1:length(genes)){
      # j=1
      
#       plot(0, type="n", main=paste0(AT.genes[g, 1] ,"\n", AT.genes[g, 2] , "\n",  genes[j]) ,xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(d.cpm.l[genes[j], ])), max(na.omit(d.cpm.l[genes[j], ]))), xlab="Time", ylab="Gene Expression in cpm", xaxt = "n")
      
      plot(0, type="n", main=paste0(AT.genes[g, 1] ,"\n", AT.genes[g, 2] , "\n",  genes[j]) ,xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(0, 10), xlab="Time", ylab="Gene Expression in cpm", xaxt = "n")
      
      
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[246])
      for(t in trees.order$legend){
        # t=trees.order$legend[1]
        lines(new.samps$time_nr[new.samps$tree_legend ==t], d.cpm.l[genes[j], new.samps$tree_legend ==t] , col=trees.order$color[trees.order$legend==t], type="b", pch=trees.order$pch[trees.order$legend==t], cex=trees.order$cex[trees.order$legend==t], lwd=3) 
      }
      axis(side=1, at=month.days[,2], labels=month.days[,1])
      legend("topright", legend = trees.order$legend, col=trees.order$color, cex=0.5, text.col=trees.order$color)
      
    }
    
  }
}


dev.off()


