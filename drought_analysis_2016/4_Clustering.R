
#####################################################################################################

### clustering

#####################################################################################################


### load data

RPath <- "/home/gosia/R/R_Shimizu_RNA_seq/Shimizu_analysis_2014-12-03/"
dataPath <- "/home/Shared/data/seq/Shimizu_RNA_seq/Data/"

analysisPath <- "Analysis_2014-12-03"
analysisPath <- paste0("/home/Shared/data/seq/Shimizu_RNA_seq/", analysisPath)
dir.create(analysisPath, showWarnings = FALSE)
setwd(analysisPath)

load("Shimizu_workspace.Rdata")



### do not consider 990 and 8212 because they were observed too few times
elim.samps <- new.samps$sample_name[grepl(pattern = "990|8212", new.samps$sample_name)]
x <- x[,!names(x) %in% elim.samps]
new.samps <- new.samps[!new.samps$sample_name %in% elim.samps, ]

all(colnames(x)==rownames(new.samps))


dim(x)


library(edgeR)
d <- DGEList(x, group=new.samps$tree_ID)
d <- calcNormFactors(d)

# make sure a gene is expressed (CPM > 1) in more than 2 samples
cps <- cpm(d, normalized.lib.sizes=TRUE)
d <- d[ rowSums( cps > 10 ) > 10, ]
dim(d$counts)

dcpm <- cpm(d, normalized.lib.sizes=TRUE)




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
source(paste0(RPath, "heatmap2.r"))



source(paste0(RPath, "Kmeans_clustering.R"))


pam.x <- dcpm
dim(pam.x)



# all
kmeans.clustering(pam.x, new.samps, genes.full.description, colors.wp=new.samps, norm.method="norm", out.path="Clustering/Kmeans_all2/", out.name="CPM_norm_all", prior.nr.cl=30:31, mc.cores=1, ylim=c(-2, 5), clustering.samps="all", colors=c("yellow","blue"))



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




