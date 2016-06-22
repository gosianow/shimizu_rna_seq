#####################################################################################################
# BioC 3.0

### Self organizing maps

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


source(paste0(RPath, "Kmeans_clustering.R"))


dcpm.norm <- normalize.counts(counts = dcpm, norm.method = "01")





library(kohonen)


# Create the SOM Grid 
som_grid <- somgrid(xdim = 50, ydim = 50, topo="hexagonal")


# Train the SOM
som_model <- som(dcpm.norm, grid=som_grid, rlen=200, keep.data = TRUE)



dir.create("Plots_SOM/", showWarnings=F, recursive=T)

library(RColorBrewer)
colgrey <- function(n) colorRampPalette(c("white", "black"))(n)
colcluster <- function(n) colorRampPalette(brewer.pal(12,"Paired"))(n)
coolBlueHotRed <- function(n, alpha = 1) rainbow(n, end=4/6, alpha=alpha)[n:1]



pdf(paste0("Plots_SOM/" , "SOM" ,".pdf"))

plot(som_model, type="changes")

plot(som_model, type="counts")

plot(som_model, type="dist.neighbours", palette.name = colgrey)

plot(som_model, type="codes")


mydata <- som_model$codes 
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var)) 
for (i in 2:40) {
  wss[i] <- sum(kmeans(mydata, centers=i, iter.max = 100)$withinss)
}
plot(wss)

nrCl <- 30
km <- kmeans(mydata, centers=nrCl, iter.max = 100)
som_cluster <- km$cluster

plot(som_model, type="mapping", bgcol = colcluster(nrCl)[som_cluster] , main = "Clusters")



### use hierarchical clustering to cluster the codebook vectors
nrCl <- 10
som_cluster <- cutree(hclust(dist(som_model$codes)), nrCl)
### plot these results:
plot(som_model, type="mapping", bgcol = colcluster(nrCl)[som_cluster] , main = "Clusters")
# add.cluster.boundaries(som_model, som_cluster)


dev.off()





pdf(paste0("Plots_SOM/" , "SOM_property" ,".pdf"))

for(var in 1:ncol(som_model$data))
plot(som_model, type = "property", property = som_model$codes[,var], main=colnames(som_model$data)[var], palette.name=coolBlueHotRed)


dev.off()




flowering <- as.numeric(genes.full.description[rownames(dcpm.norm),"Flowering"])
flowering[is.na(flowering)] <- 0


fl <- split(flowering, som_model$unit.classif)

flfull <- list()

for(i in 1:2500){
  # i = 2497
#   print(i)
  if(is.null(fl[[paste0(i)]]))
    fl[[paste0(i)]] <- 0
  flfull[[i]] <- fl[[paste0(i)]]
}



flproperty <- unlist(lapply(fl, mean))


pdf(paste0("Plots_SOM/" , "SOM_property_flowering" ,".pdf"))

  plot(som_model, type = "property", property = flproperty, main="Flowering", palette.name=coolBlueHotRed)

dev.off()





prop <- as.character(genes.full.description[rownames(dcpm.norm),"Drought_regulation"])
prop[is.na(prop)] <- 0
prop[prop == "Down"] <- 1
prop[prop == "Up"] <- 0

prop <- as.numeric(prop)

pr <- split(prop, som_model$unit.classif)

prfull <- list()

for(i in 1:2500){
  # i = 2497
  #   print(i)
  if(is.null(pr[[paste0(i)]]))
    pr[[paste0(i)]] <- 0
  prfull[[i]] <- pr[[paste0(i)]]
}



property <- unlist(lapply(prfull, mean))


pdf(paste0("Plots_SOM/" , "SOM_property_Drought_Down" ,".pdf"))

plot(som_model, type = "property", property = property, main="", palette.name=coolBlueHotRed)

dev.off()











