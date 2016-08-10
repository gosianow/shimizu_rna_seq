#####################################################################################################
### BioC 3.0
### clustering the genes from controll sets (flowering genes, drought responce, Mol Ecol clusters)
### Clustering spline fitted data 
#####################################################################################################


### load data

RPath <- "/home/gosia/R/R_Shimizu_RNA_seq/Shimizu_analysis_2014-12-03/"
dataPath <- "/home/Shared/data/seq/Shimizu_RNA_seq/Data/"

analysisPath <- "Analysis_2014-12-03"
analysisPath <- paste0("/home/Shared/data/seq/Shimizu_RNA_seq/", analysisPath)
dir.create(analysisPath, showWarnings = FALSE)
setwd(analysisPath)

load("Shimizu_workspace.Rdata")

dir.create("Plots_Clustering_GeneControlSets/", showWarnings=F, recursive=T)


source(paste0(RPath, "Kmeans_clustering_splines.R"))



##########################################################################################################
## cluster flowering genes 
##########################################################################################################


genes.fl <- as.character(genes.full.description[!is.na(genes.full.description$Flowering), "ID"])


### do not consider 990 and 8212 because they were observed too few times
elim.samps <- c(new.samps$sample_name[grepl(pattern = "990|8212", new.samps$sample_name)], new.samps$sample_name[grepl(pattern = "flower_bud", new.samps$developmental_stage)])

x <- x[,!names(x) %in% elim.samps]
new.samps <- new.samps[!new.samps$sample_name %in% elim.samps, ]

all(colnames(x)==rownames(new.samps))

trees.order <- trees.order[!grepl(pattern = "990|8212", trees.order$tree_ID), ]
trees.order <- trees.order[trees.order$legend != "8266-drought-flower_bud", ]


library(edgeR)
d.org <- DGEList(x, group=new.samps$tree_ID)
d.org <- calcNormFactors(d.org)

d.cpm <- cpm(d.org, normalized.lib.sizes=TRUE)


d.cpm <- d.cpm[genes.fl, ]
dim(d.cpm)

d.cpm <- d.cpm[ rowSums( d.cpm > 1 ) > 5, ]
dim(d.cpm)

"GID031739_3527284" %in% rownames(d.cpm) ### FT gene

d.cpm <- normalize.counts(d.cpm, norm.method= "01")


########################## fit splines 

dspl <- vector("list", 4)
dspl <- lapply(c(12, 6, 12, 11), function(n) matrix(0, nrow(d.cpm), n) )
names(dspl) <-  trees.order$legend 

genes <- rownames(d.cpm)
  
genes.full.description$ID <- as.character(genes.full.description$ID )
  
pdf(paste0("Plots_Clustering_GeneControlSets/FloweringGenes/" , "Spline_fitting_01" ,".pdf"), h=5, w=10)

for(j in 1:length(genes)){
  # j=3
  
  plot(0, type="n", main = paste0(genes.full.description[genes[j], c("ID", "AT_ID", "AT_symbol")], collapse = " \n ") , xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(d.cpm[genes[j], ])), max(na.omit(d.cpm[genes[j], ]))), xlab="Time", ylab="Gene Expression", xaxt = "n")
  axis(side=1, at=month.days[,2], labels=month.days[,1])
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[246])
  
  for(t in trees.order$legend){
    # t=trees.order$legend[2]
    
    time <- new.samps$time_nr[new.samps$tree_legend == t]
    expr <- d.cpm[genes[j], new.samps$tree_legend == t]
    
    ### plot raw expression
    lines(time, expr , col=trees.order$color[trees.order$legend==t], type="b", pch=trees.order$pch[trees.order$legend==t], cex=trees.order$cex[trees.order$legend==t]/2, lwd=1, lty = 3) 
    
    ### smooth with splines
    sm.spl <- smooth.spline(time, expr , spar=0.3)
    
    lines(sm.spl, col=trees.order$color[trees.order$legend==t], type="l", pch=trees.order$pch[trees.order$legend==t], cex=trees.order$cex[trees.order$legend==t], lwd=5) 
    
    dspl[[t]][j,] <- sm.spl$y

    
    if(j == 1)
      colnames(dspl[[t]]) <- paste0(unique(new.samps$short.name[new.samps$tree_legend == t]))  
    
  }
  legend("topright", legend = trees.order$legend, col=trees.order$color, cex=0.5, text.col=trees.order$color)
  
  
}

dev.off()


dspl <- do.call(cbind, dspl)

rownames(dspl) <- genes

new.samps.spl <- unique(new.samps[-which(colnames(new.samps) %in% c("sample_num", "sample_name", "sample_ID"))])
rownames(new.samps.spl) <- new.samps.spl$short.name
new.samps.spl <- new.samps.spl[colnames(dspl),]
new.samps.spl$sample_name <- new.samps.spl$short.name

save(dspl, new.samps.spl, file = paste0("Plots_Clustering_GeneControlSets/FloweringGenes/" , "Spline_fitting_01" ,".RData"))


########################## do clustering


# x <- dspl
# samps <- new.samps.spl
# out.path="Plots_Clustering_GeneControlSets/FloweringGenes/"
# out.name="Kmeans"
# prior.nr.cl=5:50
# mc.cores=5
# iter.max = 1000
# ylim=c(0,1)
# 
# allClusters <- kmeans.clustering.splines(x, samps, genes.full.description, out.path, out.name, prior.nr.cl, mc.cores, iter.max)  
# 
# load(paste0(out.path ,"Kmeans_allClusters.RData"))
# 
# plot.kmeans.clustering.splines(allClusters, x, samps, genes.full.description, out.path, out.name, prior.nr.cl, ylim)
# 
# plot.kmeans.clustering.splines.Exact_Expr(allClusters, x, samps, genes.full.description, out.path, out.name, prior.nr.cl, ylim, cl2plot = "all")




##########################################################################################################
## recluster MolEcol clusters
##########################################################################################################



molecol <- read.table(paste0(dataPath, "GeneControlSets/MolEcolClusters_Genes/gene_list_in_all_clusters_unigenes_DEgenes_Mol_Ecol.csv"), header = TRUE, sep=",", stringsAsFactors = FALSE) 
colnames(molecol) <- c("SB_ID", "AT_ID", "MolEcolClust")
dim(molecol) ### 1128


Sb2Sl <- read.table(paste0(dataPath, "Data_2013-07-01/reciprocal_best_hit_Sb_vs_Sl.txt"), header = FALSE, stringsAsFactors = FALSE) 
colnames(Sb2Sl) <- c("SB_ID", "ID")

dim(Sb2Sl)
length(unique(Sb2Sl$SB_ID)) ### 6829
length(unique(Sb2Sl$ID))

molecol <- merge(molecol, Sb2Sl, by="SB_ID", all.x=TRUE )

sum(!is.na(molecol$ID)) ### 986

molecol <- molecol[!is.na(molecol$ID), ]
molecol <- molecol[order(molecol$MolEcolClust, decreasing = FALSE), ]
rownames(molecol) <- molecol$ID


genes.molecol <- molecol$ID

genes.full.description <- merge(genes.full.description, molecol[,c("SB_ID", "ID", "MolEcolClust")], by="ID", all.x=TRUE)
rownames(genes.full.description) <- genes.full.description$ID



########################## prepare expression data


### do not consider 990 and 8212 because they were observed too few times
elim.samps <- c(new.samps$sample_name[grepl(pattern = "990|8212", new.samps$sample_name)], new.samps$sample_name[grepl(pattern = "flower_bud", new.samps$developmental_stage)])

x <- x[,!names(x) %in% elim.samps]
new.samps <- new.samps[!new.samps$sample_name %in% elim.samps, ]

all(colnames(x)==rownames(new.samps))

trees.order <- trees.order[!grepl(pattern = "990|8212", trees.order$tree_ID), ]
trees.order <- trees.order[trees.order$legend != "8266-drought-flower_bud", ]


library(edgeR)
d.org <- DGEList(x, group=new.samps$tree_ID)
d.org <- calcNormFactors(d.org)

d.cpm <- cpm(d.org, normalized.lib.sizes=TRUE)


d.cpm <- d.cpm[genes.molecol, ]
dim(d.cpm)

d.cpm <- d.cpm[ rowSums( d.cpm > 1 ) > 5, ]
dim(d.cpm)

"GID031739_3527284" %in% rownames(d.cpm) ### FT gene

d.cpm <- normalize.counts(d.cpm, norm.method= "01")


########################## fit splines 
dir.create("Plots_Clustering_GeneControlSets/MolEcolClusters/", showWarnings=F, recursive=T)

dspl <- vector("list", 4)
dspl <- lapply(c(12, 6, 12, 11), function(n) matrix(0, nrow(d.cpm), n) )
names(dspl) <-  trees.order$legend 

genes <- rownames(d.cpm)

genes.full.description$ID <- as.character(genes.full.description$ID )



pdf(paste0("Plots_Clustering_GeneControlSets/MolEcolClusters/" , "Spline_fitting_01" ,".pdf"), h=5, w=10)

for(j in 1:length(genes)){
  # j=3
  
  plot(0, type="n", main = paste0(genes.full.description[genes[j], c("ID", "AT_ID", "AT_symbol", "MolEcolClust")], collapse = " \n ") , xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(d.cpm[genes[j], ])), max(na.omit(d.cpm[genes[j], ]))), xlab="Time", ylab="Gene Expression", xaxt = "n", cex.main = 0.7)
  axis(side=1, at=month.days[,2], labels=month.days[,1])
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[246])
  
  for(t in trees.order$legend){
    # t=trees.order$legend[2]
    
    time <- new.samps$time_nr[new.samps$tree_legend == t]
    expr <- d.cpm[genes[j], new.samps$tree_legend == t]
    
    ### plot raw expression
    lines(time, expr , col=trees.order$color[trees.order$legend==t], type="b", pch=trees.order$pch[trees.order$legend==t], cex=trees.order$cex[trees.order$legend==t]/2, lwd=1, lty = 3) 
    
    ### smooth with splines
    sm.spl <- smooth.spline(time, expr , spar=0.3)
    
    lines(sm.spl, col=trees.order$color[trees.order$legend==t], type="l", pch=trees.order$pch[trees.order$legend==t], cex=trees.order$cex[trees.order$legend==t], lwd=5) 
    
    dspl[[t]][j,] <- sm.spl$y
    
    
    if(j == 1)
      colnames(dspl[[t]]) <- paste0(unique(new.samps$short.name[new.samps$tree_legend == t]))  
    
  }
  legend("topright", legend = trees.order$legend, col=trees.order$color, cex=0.5, text.col=trees.order$color)
  
  
}

dev.off()


dspl <- do.call(cbind, dspl)

rownames(dspl) <- genes

new.samps.spl <- unique(new.samps[-which(colnames(new.samps) %in% c("sample_num", "sample_name", "sample_ID"))])
rownames(new.samps.spl) <- new.samps.spl$short.name
new.samps.spl <- new.samps.spl[colnames(dspl),]
new.samps.spl$sample_name <- new.samps.spl$short.name

save(dspl, new.samps.spl, file = paste0("Plots_Clustering_GeneControlSets/FloweringGenes/" , "Spline_fitting_01" ,".RData"))



########################## check the expression of MolEcol clusters

x <- dspl
samps <- new.samps.spl
out.path="Plots_Clustering_GeneControlSets/MolEcolClusters/MolEcolClustering"
out.name="Kmeans"
prior.nr.cl=7
ylim=c(0,1)

molecol <- molecol[rownames(x), ]
all(rownames(x) == molecol$ID)


### create allClusters object
CL7 <- list()
CL7$clusters <- molecol$MolEcolClust

CL7$centers <- matrix(0, 7, ncol(x))
colnames(CL7$centers) <- colnames(x)
rownames(CL7$centers) <- 1:7
CL7$size <- rep(0, 7)
names(CL7$size) <- 1:7

for(i in 1:7){
  
  CL7$centers[i,] <- colMeans(x[CL7$clusters == i, ])
  CL7$size[i] <- nrow(x[CL7$clusters == i, ])
  
}

allClusters <- list()
allClusters[["CL7"]] <- CL7


plot.kmeans.clustering.splines(allClusters, x, samps, genes.full.description, out.path, out.name, prior.nr.cl, ylim)



########################## recluster MolEcol clusters


x <- dspl
samps <- new.samps.spl
out.path="Plots_Clustering_GeneControlSets/MolEcolClusters/"
out.name="Kmeans"
prior.nr.cl=3:20
mc.cores=10
iter.max = 1000
ylim=c(0,1)

allClusters <- kmeans.clustering.splines(x, samps, genes.full.description, out.path, out.name, prior.nr.cl, mc.cores, iter.max)  

load(paste0(out.path ,"Kmeans_allClusters.RData"))

plot.kmeans.clustering.splines(allClusters, x, samps, genes.full.description, out.path, out.name, prior.nr.cl, ylim)























