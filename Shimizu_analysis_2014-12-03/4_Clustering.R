
#####################################################################################################

### clustering

#####################################################################################################

setwd("/home/Shared/data/seq/Shimizu_RNA_seq/")
load("Shimizu_workspace.Rdata")


############################################
### select genes for clustering
############################################

x.orig <- x
new.samps.orig <- new.samps

# elim.samps=c(flowered.samps, "K4_990_20081205", "K5_990_20090511")
elim.samps=NULL
x.orig <- x.orig[,new.samps.orig$sample_name]
x <- x.orig[,!colnames(x.orig) %in% elim.samps]
new.samps <- new.samps.orig[!rownames(new.samps.orig) %in% elim.samps, ]

all(colnames(x)==rownames(new.samps))

library(limma)
library(edgeR)

d.org <- DGEList(x, group=new.samps$tree_ID)

# normalization between samples, default = "TMM"
d.org <- calcNormFactors(d.org)

### make sure a gene is expressed (CPM > 1) in more than 2 samples
d.cpm.org <- cpm(d.org, normalized.lib.sizes=TRUE)
dim(d.cpm.org)

# sum(d.cpm.org["GID031739_3527284",] > 1)
# sum( rowSums(d.cpm.org > 1) > 2 )

d.org <- d.org[ rowSums(d.cpm.org > 1) > 2, ]
d.cpm.org <- cpm(d.org, normalized.lib.sizes=TRUE)
d.cpm.org.l <- log(d.cpm.org  + min(d.cpm.org [d.cpm.org  != 0]))

count.data <- list()
count.data$d.cpm.org <- d.cpm.org
count.data$d.cpm.org.l <- d.cpm.org.l
selected.genes <- list()
selected.genes$all.genes <- rownames(d.cpm.org)
selected.genes$AT.genes <- selected.genes$all.genes[selected.genes$all.genes %in% AT.id[,1]]


write.table(data.frame(contigs=rownames(d.org$counts) , d.org$counts), "Clustering/Data_original.xls", quote=FALSE, sep="\t", row.names=FALSE)


write.table(data.frame(contigs=rownames(d.cpm.org) , d.cpm.org), "Clustering/Data_normalized_cpm.xls", quote=FALSE, sep="\t", row.names=FALSE)



############################################
# Kmeans run
############################################

dir.create("Clustering", showWarnings=F, recursive=T)

# install.packages("/home/gosia/R/packages/", "/home/gosia/R/libraries/3.0.2/")
library(cluster)
library(clusterSim)
library(lattice)
library(parallel)
library(gtools)
library(gplots) 
library(RColorBrewer)
source("/home/gosia/R/R_Shimizu_RNA_seq/heatmap2.r")



source("/home/gosia/R/R_Shimizu_RNA_seq/Kmeans_clustering.R")
# source("/home/gosia/R/R_Shimizu_RNA_seq/Kmeans_clustering_tmp.R")

pam.x <- count.data$d.cpm.org
#pam.x <- pam.x[1:200,]
dim(pam.x)



# all
kmeans.clustering(pam.x, new.samps, genes.full.description, colors.wp, norm.method="norm", out.path="Clustering/Kmeans_all2/", out.name="CPM_norm_all", prior.nr.cl=10:50, mc.cores=15, ylim=c(-2, 5), clustering.samps="all", colors=c("yellow","blue"))



# drought trees (no flowering)
drought.samp <- new.samps[new.samps$tree_ID %in% c("970", "8266"), "sample_name"]
# no flowering sample E7_8266_20090416
drought.samp <- drought.samp[!drought.samp %in% c("E7_8266_20090416")]

kmeans.clustering(pam.x[,drought.samp], new.samps[drought.samp,], genes.full.description, colors.wp, norm.method="norm", out.path="Clustering/Kmeans_drought_trees/", out.name="CPM_norm_drought", prior.nr.cl=5:40, mc.cores=15, ylim=c(-2, 5), clustering.samps="all", colors=c("yellow","blue"))




# per tree

kmeans.clustering(pam.x, new.samps, genes.full.description, colors.wp, norm.method="norm", out.path="Clustering/Kmeans_per_tree/", out.name="8266", prior.nr.cl=10:20, mc.cores=5, ylim=c(-2, 5), clustering.samps=new.samps$sample_name[new.samps$tree_ID=="8266"], colors=c("yellow","blue"))


kmeans.clustering(pam.x, new.samps, genes.full.description, colors.wp, norm.method="norm", out.path="Clustering/Kmeans_per_tree/", out.name="970", prior.nr.cl=5:20, mc.cores=5, ylim=c(-2, 5), clustering.samps=new.samps$sample_name[new.samps$tree_ID=="970"], colors=c("yellow","blue"))


kmeans.clustering(pam.x, new.samps, genes.full.description, colors.wp, norm.method="norm", out.path="Clustering/Kmeans_per_tree/", out.name="1099", prior.nr.cl=5:20, mc.cores=5, ylim=c(-2, 5), clustering.samps=new.samps$sample_name[new.samps$tree_ID=="1099"], colors=c("yellow","blue"))


kmeans.clustering(pam.x, new.samps, genes.full.description, colors.wp, norm.method="norm", out.path="Clustering/Kmeans_per_tree/", out.name="1377", prior.nr.cl=5:20, mc.cores=5, ylim=c(-2, 5), clustering.samps=new.samps$sample_name[new.samps$tree_ID=="1377"], colors=c("yellow","blue"))




