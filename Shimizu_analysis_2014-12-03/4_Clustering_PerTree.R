

#####################################################################################################
### BioC 3.0
### clustering per tree

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

load(paste0("Plots_MeanForReplicates/" , "data.reduced" ,".RData"))

out.dir <- "Plots_PerTree_Clustering/"
dir.create(out.dir, showWarnings=F, recursive=T)

source(paste0(RPath, "Kmeans_clustering_splines.R"))


library(edgeR)


#####################################
gs.mDr <- read.table(paste0(dataPath, "GeneControlSets/gene_sets_for_Gosia/mDr_Day10_drought_down_regulated_genes_Harb_etal.csv"), header=T, sep=";")
dim(gs.mDr)

trees <- as.character(unique(new.samps.red$tree_ID))


for(i in 1:length(trees)){
  # i = 3
  
  x <- x.red[, new.samps.red$tree_ID == trees[i]]
  samps <- new.samps.red[new.samps.red$tree_ID == trees[i], ]
  
  ### Filtering
  
  x <- x[ MinCPM.tree[, trees[i]] & FC.tree[, trees[i]] & MinABS.tree[, trees[i]], ]
  dim(x)
  
  d <- DGEList(x)
  d <- calcNormFactors(d)
  dcpm <- cpm(d, normalized.lib.sizes=TRUE)
  
  x <- normalize.counts(counts = dcpm, norm.method = "01")
  
  cat("mDr genes: ",sum(genes.full.description$AT_ID %in% gs.mDr[,1] & genes.full.description$ID %in% rownames(x)), fill = TRUE)
  
  out.path=paste0(out.dir, "/Tree", trees[i],"/")
  out.name="Kmeans"
  prior.nr.cl=5:25
  mc.cores=4
  ylim=c(0,1)
  
  allClusters <- kmeans.clustering.splines(x, samps, genes.full.description, out.path, out.name, prior.nr.cl, mc.cores, iter.max = 1000)  
  
# load(paste0(out.path ,"Kmeans_allClusters.RData"))
  
  plot.kmeans.clustering.splines(allClusters, x, samps, genes.full.description, out.path, out.name, prior.nr.cl, ylim)
  
   
}

















