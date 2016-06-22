
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

load(paste0("Plots_Splines/Spline_fitting_01/" , "Spline_fitting_01.RData"))

dir.create("Plots_Splines_Clustering/", showWarnings=F, recursive=T)

source(paste0(RPath, "Kmeans_clustering_splines.R"))

########################################################################
##################### clustering splines
########################################################################
# x <- dspl
# samps <- new.samps.spl
# out.path="Plots_Splines_Clustering/Kmeans100/"
# out.name="Kmeans"
# prior.nr.cl=5:50
# ylim=c(0,1)


x <- dspl
samps <- new.samps.spl
out.path="Plots_Splines_Clustering/Splines/"
out.name="Kmeans"


allClusters <- kmeans.clustering.splines(x, samps, genes.full.description, out.path, out.name="Kmeans", prior.nr.cl=5:200, mc.cores=15, iter.max = 1000)  


load(paste0(out.path ,"Kmeans_allClusters.RData"))


plot.kmeans.clustering.splines(allClusters, x, samps, genes.full.description, out.path, out.name="Kmeans", prior.nr.cl=5:200, ylim=c(0,1))


plot.kmeans.clustering.splines.Exact_Expr(allClusters, x, samps, genes.full.description, out.path, out.name="Kmeans", prior.nr.cl=5:50, ylim=c(0,1), cl2plot = "CL35")




########################################################################
##################### clustering splines normalized (that the expression for the first time point is 0)
########################################################################

xzero <- normalize.2zero.start(x = dspl, samps = new.samps.spl)

x <- xzero
samps <- new.samps.spl
out.path="Plots_Splines_Clustering/ZeroStart/"
out.name="Kmeans"
prior.nr.cl=5:50

allClusters <- kmeans.clustering.splines(x, samps, genes.full.description, out.path, out.name="Kmeans", prior.nr.cl=5:50, mc.cores=10, iter.max = 100)  


load(paste0(out.path ,"Kmeans_allClusters.RData"))


plot.kmeans.clustering.splines(allClusters, x, samps, genes.full.description, out.path, out.name="Kmeans", prior.nr.cl=5:50, ylim=c(-0.5, 0.5))



########################################################################
##################### clustering splines normalized (that the minimal expression is 0)
########################################################################

xzero <- normalize.2zero.min(x = dspl, samps = new.samps.spl)

x <- xzero
samps <- new.samps.spl
out.path="Plots_Splines_Clustering/ZeroMin/"
out.name="Kmeans"


allClusters <- kmeans.clustering.splines(x, samps, genes.full.description, out.path, out.name="Kmeans", prior.nr.cl=5:100, mc.cores=10, iter.max = 1000)  


load(paste0(out.path ,"Kmeans_allClusters.RData"))


plot.kmeans.clustering.splines(allClusters, x, samps, genes.full.description, out.path, out.name="Kmeans", prior.nr.cl=5:100, ylim=c(0, 0.5))



######### check genes with small sd

sd.xzero <- apply(xzero, 1, sd)
hist(sd.xzero)
sd.xzeroSort <- sort(sd.xzero, decreasing = FALSE)

pdf(paste0(out.path, out.name,"_Avg_Expr.pdf"), w=10, h=5) 

gene <- names(sd.xzeroSort)[20]
expr <- dspl[gene, ]
plot.expr(new.samps = new.samps.spl, expr, trees.order, main=gene, ylim=c(0, 1), month.days)
expr <- xzero[gene, ]
plot.expr(new.samps = new.samps.spl, expr, trees.order, main=gene, ylim=c(0, 1), month.days)

dev.off()




########################################################################
##################### clustering splines normalized (that the minimal expression is 0)
##################### remove genes with small variance TODO
########################################################################

xzero <- normalize.2zero.min(x = dspl, samps = new.samps.spl)
out.path="Plots_Splines_Clustering/ZeroMin_BigSD/"
out.name="Kmeans"

dir.create(out.path, showWarnings=F, recursive=T)


######### check genes with small sd

sd.xzero <- apply(xzero, 1, sd)

pdf(paste0(out.path, out.name,"_HistSD.pdf"), w=10, h=5) 
hist(sd.xzero)
dev.off()
sd.xzeroSort <- sort(sd.xzero, decreasing = FALSE)

sd.keep <- which( sd.xzeroSort > 0.1 )
length(sd.keep)


######### check genes with small FC (calculate on dcpm)

fc.xzero <- apply(xzero, 1, function(g){
  
  max(g)/min(g)
  
})



pdf(paste0(out.path, out.name,"_Avg_Expr.pdf"), w=10, h=5) 

gene <- names(sd.keep)[1]
expr <- dspl[gene, ]
plot.expr(new.samps = new.samps.spl, expr, trees.order, main=gene, ylim=c(0, 1), month.days)
expr <- xzero[gene, ]
plot.expr(new.samps = new.samps.spl, expr, trees.order, main=gene, ylim=c(0, 1), month.days)

dev.off()



######### cluster

x <- xzero
samps <- new.samps.spl




allClusters <- kmeans.clustering.splines(x, samps, genes.full.description, out.path, out.name="Kmeans", prior.nr.cl=5:100, mc.cores=10, iter.max = 1000)  


load(paste0(out.path ,"Kmeans_allClusters.RData"))


plot.kmeans.clustering.splines(allClusters, x, samps, genes.full.description, out.path, out.name="Kmeans", prior.nr.cl=5:100, ylim=c(0, 0.5))




########################################################################
##################### clustering splines normalized (that the minimal expression is 0)
##### Per tree 
########################################################################

xzero <- normalize.2zero.min(x = dspl, samps = new.samps.spl)
 dim(xzero)
cl.tree <- "8266"

samps <- new.samps.spl
samps <- samps[samps$tree_ID == cl.tree, ]

x <- xzero[ rownames(xzero) %in% rownames(FC.tree[FC.tree[, cl.tree] > 3, ]) , samps$sample_name]
dim(x)


out.path=paste0("Plots_Splines_Clustering/ZeroMin_PerTree/", cl.tree, "/")
out.name="Kmeans"
prior.nr.cl=5:50
mc.cores=20
ylim=c(0, 0.7)

allClusters <- kmeans.clustering.splines(x, samps, genes.full.description, out.path, out.name, prior.nr.cl, mc.cores, iter.max = 1000)  

load(paste0(out.path ,"Kmeans_allClusters.RData"))


plot.kmeans.clustering.splines(allClusters, x, samps, genes.full.description, out.path, out.name, prior.nr.cl, ylim)











































