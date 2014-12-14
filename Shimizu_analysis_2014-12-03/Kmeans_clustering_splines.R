library(parallel)
library(clusterSim)

# x <- dspl[1:1000,]
# samps <- new.samps.spl
# out.path="Plots_Splines_Clustering/Test/"
# out.name=""
# prior.nr.cl=2:5 
# mc.cores=3
# iter.max <- 100


kmeans.clustering.splines <- function(x, samps, genes.full.description, out.path, out.name="Kmeans", prior.nr.cl=2:5, mc.cores=4, iter.max = 100)  
{
  
  dir.create(out.path, showWarnings=FALSE, recursive=TRUE)
  
  ### order samples 
  samps.order <- row.names(samps[order( factor(samps$developmental_stage, levels=as.character(sort(unique(samps$developmental_stage), decreasing = T))), samps$drough.control, samps$tree_ID, samps$time_nr), ])
  x <- x[, samps.order]
  samps <- samps[samps.order,]
  
  
  allClusters <- mclapply(1:length(prior.nr.cl), function(i){
    # i=1
    cat("Clustering into ", prior.nr.cl[i], "clusters \n")
    cl <- kmeans(x=x, centers=prior.nr.cl[i], iter.max = iter.max , nstart = 2) 
    return(cl)
    
  }, mc.cores=mc.cores)
  
  names(allClusters) <- paste0("CL", prior.nr.cl)
  
  df.cl <- matrix(0, nrow(x), length(prior.nr.cl))
  df.cent <- list()
  
  for(i in 1:length(allClusters)){
    cl <- allClusters[[i]]
    df.cl[, i] <- cl$cluster
    rownames(cl$centers) <- paste0("CL",prior.nr.cl[i],":",1:prior.nr.cl[i])
    df.cent[[i]] <- cl$centers
  }
  
  df.cent <- do.call(rbind, df.cent)
  rownames(df.cl) <- rownames(x)
  colnames(df.cl) <- paste0("CL", prior.nr.cl)
  
  write.table(df.cl, paste0(out.path,"/",out.name,"_all_clustering.txt"), quote = FALSE, sep = "\t" )
  write.table(df.cent, paste0(out.path,"/",out.name,"_all_clustering_centers.txt"), quote = FALSE, sep = "\t")
  
  
  cat("Calculate GAP indexes... \n")
  
  GapIndx <- mclapply(1:(length(allClusters)-1), function(p){ 
    # p=2   
    cat("... gap index for", prior.nr.cl[p] , "...\n")
    
    clall <- cbind(allClusters[[p]]$cluster, allClusters[[p+1]]$cluster)
    
    pam.indexGap.pc  <- index.Gap(x=x, clall=clall , method="k-means", centrotypes="centroids", reference.distribution="pc", B=100)$diffu
    pam.indexGap.unif <- index.Gap(x=x, clall=clall, method="k-means", centrotypes="centroids", reference.distribution="unif", B=100)$diffu
    
    return(c(pam.indexGap.pc, pam.indexGap.unif))
    
  }, mc.cores=mc.cores)
  
  
  GapIndx <- do.call(rbind, GapIndx)
  colnames(GapIndx) <- paste0("GapIndx", c("PC", "UNIF"))
  rownames(GapIndx) <- paste0("CL",prior.nr.cl[-length(prior.nr.cl)])
  
  nc <- prior.nr.cl[-length(prior.nr.cl)]
  
  pdf(paste0(out.path,"/",out.name,"_idexes.pdf"))

  plot(nc, GapIndx[,1], main= paste0("Gap pc statistic \n final nr of clusters = min{nc: diffu(nc, nc+1) >= 0} \n "),  xlab="nc", ylab="diffu", pch=19, col=ifelse(GapIndx[,1] >= 0, 2,1), xaxt = "n")
  axis(side=1, at=nc, labels=nc)
  abline(h=0)
  abline(v=nc, lty=2, lwd=0.5, col=ifelse(GapIndx[,1] >= 0, 2,"grey"))
  
  plot(nc, GapIndx[,2], main=paste0("Gap unif statistic \n final nr of clusters = min{nc: diffu(nc, nc+1) >= 0} \n "), xlab="nc", ylab="diffu", pch=19, col=ifelse(GapIndx[,2] >= 0, 2,1), xaxt = "n")
  axis(side=1, at=nc, labels=nc)
  abline(h=0)
  abline(v=nc, lty=2, lwd=0.5, col=ifelse(GapIndx[,2] >= 0, 2,"grey"))
  
  dev.off()
  
  save(allClusters, file = paste0(out.path,"/",out.name,"_allClusters.RData"))
  
  return(invisible(allClusters))
  
}



library(reshape2)
library(ggplot2)
library(scales)
library(ggdendro)
library(lattice)
library(parallel)
library(gtools)
library(gplots) 
library(RColorBrewer)

source(paste0(RPath, "heatmap2centers.R"))


plot.kmeans.clustering.splines <- function(allClusters, x, samps, genes.full.description, out.path, out.name="Kmeans", prior.nr.cl=2:5, ylim=c(0, 1) ){
  
  cl.names <- names(allClusters)
  
  
  for(i in 1:length(allClusters)){
    # i = 1
    
    out.cl <- paste0(out.path,"/", cl.names[i])
#     out.cl <- paste0(out.path,"/", cl.names[i],"/") ### create folder for each custering
#     dir.create(out.cl, showWarnings = FALSE, recursive = TRUE)
    
    
    ##### heatmaps
    
    cent <- allClusters[[i]]$centers
    rownames(cent) <- paste0("CL", 1:nrow(cent))
    
    cent.melt <- melt(cent)
    colnames(cent.melt) <- c("X1", "X2", "value")
    cent.melt$X1 <- factor(cent.melt$X1, levels = unique(cent.melt$X1))
    
    
# myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
  myPalette <- colorRampPalette(c("green", "black", "red"))


# pdf(paste0(out.cl, out.name,"_Heatmap.pdf"), w=10, h=5) 
# p <- ggplot(cent.melt, aes(x = X2, y = X1, fill = value)) + geom_tile() + scale_fill_gradientn(colours = myPalette(100))
# p + labs(x = "Samples", y = "Clusters") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))  + theme_bw() + theme(axis.text.x = element_text(angle = 330, hjust = 0, colour = new.samps.spl$tree_col, size = 8))
# dev.off()


######################## cluster the groups (to see how similar they are)

hc <- hclust(dist(cent))

cent.hc <- cent[hc$order, ]  
cl.ord <- rownames(cent.hc)
cent.melt <- melt(cent.hc)
colnames(cent.melt) <- c("X1", "X2", "value")
cent.melt$X1 <- factor(cent.melt$X1, levels = cl.ord)



ggp <- ggplot(cent.melt, aes(x = X2, y = X1, fill = value)) + geom_tile() + scale_fill_gradientn(colours = myPalette(100)) + labs(x = "Samples", y = "Clusters") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))  + theme_bw() + theme(axis.text.x = element_text(angle = 330, hjust = 0, colour = new.samps.spl$tree_col, size = 9))
ggsave(paste0(out.cl, out.name,"_Heatmap_hcOrder.pdf"), width = 13, height = 10, units = "in")

######### plot dendogram


ggd <- ggdendrogram(hc, size=4, theme_dendro = TRUE, rotate=TRUE)
ggsave(paste0(out.cl, out.name,"_Heatmap_hcOrder_dendo.pdf"), width = 10, height = 10, units = "in")


######### plot heatmap & dendogram

pdf(paste0(out.cl, out.name,"_Heatmap&dendo.pdf"), w=15, h=11) 

heatmap.2centers(cent, Rowv=as.dendrogram(hc), Colv=NA, dendrogram = "row", scale = "none", trace="none", col=colorRampPalette(c("green", "black", "red"))(250), density.info="none", xlab = "", ylab = "", key = TRUE, keysize = 0.5, margins = c(12, 5), cexCol = 1.3, cexRow = 1, col.axis.x = new.samps.spl$tree_col, sepcolor = NA, sepwidth = c(0, 0))

# heatmap(cent, Rowv=as.dendrogram(hc), Colv=NA, scale = "none", col=colorRampPalette(c("green", "black", "red"))(250), ColSideColors = new.samps.spl$tree_col)

dev.off()




###### plot gene expression in each cluster separetely 

pdf(paste0(out.cl, out.name,"_Avg_Expr.pdf"), w=10, h=5) 

### to keep the order from top to bottom on heatmap
cl.ordR <- rev(cl.ord)
cent.hcR <- cent.hc[cl.ordR, ]

all(rownames(new.samps.spl) == colnames(cent.hcR))

for(c in 1:length(cl.ordR)){
  # c=1
  
  plot(0, type="n", main=paste0("Cluster ", cl.ordR[c]) , xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=ylim, xlab="Time", ylab="Normalized Gene Expression", xaxt = "n")
  axis(side=1, at=month.days[,2], labels=month.days[,1])
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[246])
   
  
  trees.order.plot <- trees.order[trees.order$legend %in% unique(new.samps.spl$tree_legend), ]
  
  for(t in trees.order.plot$legend){
    
    lines(new.samps.spl$time_nr[new.samps.spl$tree_legend == t], cent.hcR[ c, new.samps.spl$tree_legend == t] , col=trees.order.plot$color[trees.order.plot$legend==t], type="l", pch=trees.order.plot$pch[trees.order.plot$legend==t], cex=trees.order.plot$cex[trees.order.plot$legend==t], lwd=5)
    
  }
  
  
  legend("topright", legend = trees.order.plot$legend, col=trees.order.plot$color, cex=0.8, text.col=trees.order.plot$color) 
  
}

dev.off()


    
    
  }
  

  
}



normalize.2zero.start <- function(x, samps){
  
  trees <- unique(samps$tree_ID)
  
  x.t.f <- NULL
  
  for(i in trees){
    # i = trees[1]
    
    x.t <- x[, samps$tree_ID == i]
    
    for(j in 1:nrow(x.t))
      x.t[j,] <- x.t[j,] - x.t[j,1]
    
    x.t.f <- cbind(x.t.f, x.t)
    
  }
  
  return(x.t.f)
  
}

























