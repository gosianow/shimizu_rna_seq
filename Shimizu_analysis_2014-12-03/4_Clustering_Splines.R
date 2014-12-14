
#####################################################################################################
### BioC 3.0
### clustering the splines fitting

#####################################################################################################


### load data

RPath <- "/home/gosia/R/R_Shimizu_RNA_seq/Shimizu_analysis_2014-12-03/"
dataPath <- "/home/Shared/data/seq/Shimizu_RNA_seq/Data/"

analysisPath <- "Analysis_2014-12-03"
analysisPath <- paste0("/home/Shared/data/seq/Shimizu_RNA_seq/", analysisPath)
dir.create(analysisPath, showWarnings = FALSE)
setwd(analysisPath)

load("Shimizu_workspace.Rdata")
load(paste0("Plots_Splines/" , "Spline_fitting_01" ,".RData"))

dir.create("Plots_Splines_Clustering/", showWarnings=F, recursive=T)



##################### clustering splines

source(paste0(RPath, "Kmeans_clustering_splines.R"))

x <- dspl
samps <- new.samps.spl
out.path="Plots_Splines_Clustering/Kmeans100/"
out.name="Kmeans"
prior.nr.cl=5:50

# x <- dspl[1:1000, ]
# samps <- new.samps.spl
# out.path="Plots_Splines_Clustering/Test/"
# out.name="Kmeans"
# prior.nr.cl=5:50


allClusters <- kmeans.clustering.splines(x, samps, genes.full.description, out.path, out.name="Kmeans", prior.nr.cl=5:50, mc.cores=15, iter.max = 100)  


load(paste0(out.path ,"Kmeans_allClusters.RData"))


plot.kmeans.clustering.splines(allClusters, x, samps, genes.full.description, out.path, out.name="Kmeans", prior.nr.cl=5:50)


##################### clustering splines normalized (that the expression for the first time point is 0)


xzero <- normalize.2zero.start(x = dspl, samps = new.samps.spl)

x <- xzero
samps <- new.samps.spl
out.path="Plots_Splines_Clustering/ZeroStart/"
out.name="Kmeans"
prior.nr.cl=5:50

allClusters <- kmeans.clustering.splines(x, samps, genes.full.description, out.path, out.name="Kmeans", prior.nr.cl=5:50, mc.cores=10, iter.max = 100)  


load(paste0(out.path ,"Kmeans_allClusters.RData"))


plot.kmeans.clustering.splines(allClusters, x, samps, genes.full.description, out.path, out.name="Kmeans", prior.nr.cl=5:50, ylim=c(-0.3, 0.3))




















































