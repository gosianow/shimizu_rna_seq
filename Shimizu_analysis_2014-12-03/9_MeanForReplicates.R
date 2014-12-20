#####################################################################################################
# BioC 3.0

### Fit splines 

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

out.dir <- "Plots_MeanForReplicates/"
dir.create(out.dir, showWarnings=F, recursive=T)


source(paste0(RPath, "Kmeans_clustering_splines.R"))

load(paste0("Plots_VarianceMeanTrend/", "ABSMeanScale.RData"))



###############################################

### do not consider 990 and 8212 because they were observed too few times
elim.samps <- c(new.samps$sample_name[grepl(pattern = "990|8212", new.samps$sample_name)], new.samps$sample_name[grepl(pattern = "flower_bud", new.samps$developmental_stage)])

x <- x[,!names(x) %in% elim.samps]
new.samps <- new.samps[!new.samps$sample_name %in% elim.samps, ]


all(colnames(x)==rownames(new.samps))


trees.order <- trees.order[!grepl(pattern = "990|8212", trees.order$tree_ID), ]
trees.order <- trees.order[trees.order$legend != "8266-drought-flower_bud", ]

dim(x)


library(edgeR)
d <- DGEList(x, group=new.samps$tree_ID)
d <- calcNormFactors(d)


###############################################
# Filters 
###############################################

dcpm <- cpm(d, normalized.lib.sizes=TRUE)

trees <- unique(new.samps$tree_ID)

### Min FC
FC.tree <- matrix(0, nrow(dcpm), length(trees))
colnames(FC.tree) <- trees
rownames(FC.tree) <- rownames(dcpm)

### Min CPM
MinCPM.tree <- matrix(0, nrow(dcpm), length(trees))
colnames(MinCPM.tree) <- trees
rownames(MinCPM.tree) <- rownames(dcpm)


### Min ABS
MinABS.tree <- matrix(0, nrow(dcpm), length(trees))
colnames(MinABS.tree) <- trees
rownames(MinABS.tree) <- rownames(dcpm)



for(i in 1:length(trees)){
  # i = 1  
  
  dcpm.t <- dcpm[, new.samps$tree_ID == trees[i] ]
  
  ### FC
  fc.all <- apply(dcpm.t, 1, function(g){
    
    if(max(g)  == 0)
      fc <- 0
    else
      fc <- max(g) / min(g[ g != 0 ])
    
    return(fc)
  })
  
  FC.tree[,i] <- fc.all > 2
  
  
  ### CPM  
  MinCPM.tree[,i] <- rowSums( dcpm.t > 1 ) > (dim(dcpm.t)[2]-1) 

  
  ### ABS
  abs.all <- apply(dcpm.t, 1, function(g){
    # g = dcpm.t[100, ]
    
    if(!sum( g > 1 ) > (dim(dcpm.t)[2]-1))
      return(FALSE)
    else
      ps <- predict(s, min(g))$y
    
    abs <- (max(g) - min(g)) > 3*ps
    
    return(abs)
  })
  
  sum(abs.all)
  
  MinABS.tree[,i] <- abs.all

  
}



###############################################
# Mean Expr for replicates
###############################################


pairs <- list(c("H6_8266_20081208", "H7_8266_20081208"), c("H8_8266_20090105", "H9_8266_20090105"), c("E4_8266_20090202", "F4_8266_20090202"), c("N1_8266_20090224", "N5_8266_20090224"), c("I9_970_20081216", "J1_970_20081216"), c("J2_970_20090207", "K9_970_20090207"))



x.red <- matrix(0, nrow(x), length(pairs))
colnames(x.red) <- unlist(pairs)[rep(c(TRUE, FALSE), length(pairs))]


for(i in 1:length(pairs)){
  
  x.red[, i] <- rowMeans(x[, pairs[[i]]])
  
}

x.red <- cbind(x.red, x[, !colnames(x) %in% unlist(pairs)])

new.samps.red <- new.samps[colnames(x.red), ]
new.samps.red$sample_name <- new.samps.red$short.name
rownames(new.samps.red) <- new.samps.red$short.name

colnames(x.red) <- new.samps.red$short.name

### order samples 
samps.order <- row.names(new.samps.red[order( factor(new.samps.red$developmental_stage, levels=as.character(sort(unique(new.samps.red$developmental_stage), decreasing = T))), new.samps.red$drough.control, new.samps.red$tree_ID, new.samps.red$time_nr), ])
x.red <- x.red[, samps.order]
new.samps.red <- new.samps.red[samps.order,]


write.table(x.red, paste0(out.dir , "x.red" ,".txt"), quote = FALSE, sep="\t")
write.table(new.samps.red, paste0(out.dir , "new.samps.red" ,".txt"), quote = FALSE, sep="\t")

save(x.red, new.samps.red, FC.tree, MinCPM.tree, MinABS.tree, file = paste0(out.dir , "data.reduced" ,".RData"))























