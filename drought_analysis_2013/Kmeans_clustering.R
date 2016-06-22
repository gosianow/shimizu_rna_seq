############################################
### help functions:
############################################



plot.mean.sd <- function(d){ 
  d.mean <- apply(d, 1, mean)
  d.sd <- apply(d, 1, sd)
  plot(d.mean, d.sd)
}

normalize.counts <- function(counts, norm.method=c("norm", "01", "range")[1], c=0, d=1){
  
  t(apply(counts, 1, function(g){
    
    if(norm.method=="norm"){
      m <- mean(g, na.rm=T)
      sd <- sd(g, na.rm=T) 
      return((g-m)/sd)
    } else if(norm.method=="01"){
      return((g-min(g,na.rm=T))/(max(g,na.rm=T)-min(g,na.rm=T)))
    }
    else if(norm.method=="range"){
      return(c*(1-(g-min(g,na.rm=T))/(max(g,na.rm=T)-min(g,na.rm=T)))+d*((g-min(g,na.rm=T))/(max(g,na.rm=T)-min(g,na.rm=T))))
    }
  }))
}


plot.error.bars <- function(x, y, delta, col=1, barwidth=100){
  segments(x,y-delta,x,y+delta, col=col)
  segments(x-barwidth,y+delta,x+barwidth,y+delta, col=col)
  segments(x-barwidth,y-delta,x+barwidth,y-delta, col=col)
}


############################################
### kmeans.clustering
############################################


# pam.x <- count.data$d.cpm.org
# pam.x <- pam.x[1:200,]
# 
# norm.method="norm"
# out.path="Clustering/Kmeans_test/"
# out.name=""
# prior.nr.cl=2:5
# mc.cores=1
# ylim=c(-2,6)
# clustering.samps=drought.samp
# colors=c("yellow","green4")


kmeans.clustering <- function(pam.x, new.samps, genes.full.description, colors.wp, norm.method="norm", out.path="Clustering/Kmeans/", out.name="", prior.nr.cl=2:5, mc.cores=1, ylim=c(-2,6), clustering.samps="all", colors=c("yellow","green4"))
{
  
  dir.create(out.path, showWarnings=FALSE, recursive=TRUE)
  
  samps.order <- row.names(new.samps[order( factor(new.samps$developmental_stage, levels=as.character(sort(unique(new.samps$developmental_stage), decreasing = T))), new.samps$drough.control, new.samps$tree_ID, new.samps$time_nr), ])
  pam.x <- pam.x[, samps.order]
  new.samps <- new.samps[samps.order,]
 
  if(!any(clustering.samps=="all"))
    new.samps <- new.samps[new.samps$sample_name %in% clustering.samps,]
  
  if(norm.method!="none")
    pam.x <- normalize.counts(counts=pam.x, norm.method=norm.method)
  
  pam.x <- na.omit(pam.x) 
  pam.df <- as.data.frame(pam.x)
  
  if(!any(clustering.samps=="all"))
    pam.df <- as.data.frame(pam.x[, new.samps$sample_name])
  
  # cat("Clac dist matrix \n")  
  # pam.dist <- dist(x=pam.x , method = "euclidean")

  
  pam.obj.all <- mclapply(prior.nr.cl, function(i){
    # i=3
    cat("Clustering into ", i, "clusters \n")
    pam.obj <- kmeans(x=pam.x, centers=i, iter.max = 100, nstart = 2) 
    pam.cl <- pam.obj$cluster
    
    # Calinski-Harabasz index
    # pam.obj$index.G1 <- index.G1(x=pam.x, cl=pam.cl, centrotypes="centroids")
    # Silhouette index
    # pam.obj$index.S <- index.S(d=pam.dist, cl=pam.cl)   
    # tot.withinss
    # pam.obj$index.S <- log(pam.obj$tot.withinss)
       
    pam.df[, paste0("clustering", i)] <- pam.cl
    
    pam.df.sort <- pam.df[order(pam.df[, paste0("clustering", i)]), ]
    heat.matrix <- as.matrix(pam.df.sort[,-ncol(pam.df.sort)])   
    colnames(heat.matrix) <- new.samps$short.name
    
    # ColSideColors = tree colors
    #     png(paste0(out.path,"/",out.name,"_Heatmap_k" ,i, ".png"), w=1000, h=1000)
    #     heatmap.2axis(x=heat.matrix, Rowv=NA, Colv=NA, dendrogram="none", labRow=NA, labCol=NULL, scale="none", trace="none", col=colorRampPalette(colors)(250), key=T, keysize =c(0.5, 0.5), density.info="none", margins = c(5,5), xlab = "Samples", ylab = "Clusters of DE Genes", RowSideColors=colorRampPalette(brewer.pal(12,"Set3"))(i)[pam.df.sort[, ncol(pam.df.sort)]], ColSideColors=new.samps$tree_col) 
    #     dev.off()
      
    
    lab.cl <- pam.df.sort[, ncol(pam.df.sort)] 
    labRow <- matrix(NA, 1, length(lab.cl) )
    
    for( c in unique(lab.cl)){
      # c=1
      c.ind <- which(lab.cl==c)
      c.ind <- c.ind[round(length(c.ind)/2)]
      labRow[c.ind] <- c          
    }

    
    # ColSideColors = heat map of WP
    png(paste0(out.path,"/",out.name,"_Heatmap_k" ,i, ".png"), w=1000, h=1000)
    heatmap.2axis(x=heat.matrix, Rowv=NA, Colv=NA, dendrogram="none", labRow=labRow, labCol=NULL, scale="none", trace="none", col=colorRampPalette(colors)(250), key=T, keysize =c(0.5, 0.5), density.info="none", margins = c(5,5), xlab = "Samples", ylab = "Clusters of DE Genes", RowSideColors=colorRampPalette(brewer.pal(12,"Set3"))(i)[pam.df.sort[, ncol(pam.df.sort)]], ColSideColors=colors.wp[new.samps$sample_name,"colors.wp"], cexRow=2) 
    dev.off()
    

    
    if(length(ylim)==2){
      min.pam.df <- ylim[1]
      max.pam.df <- ylim[2]   
    } else {
      min.pam.df <- min(pam.df[, -ncol(pam.df)])
      max.pam.df <- max(pam.df[, -ncol(pam.df)])
    }
    

pdf(paste0(out.path,"/",out.name,"_Avg_Expression_in_cluster" ,i, ".pdf"), w=10, h=5) 
for(c in 1:i){
  # c=1
  pam.df.tmp <- pam.df[pam.cl==c, -ncol(pam.df)]
  pam.df.tmp <-  pam.df.tmp[, rownames(new.samps)]
  
  pam.df.tmp.avg <- apply(pam.df.tmp, 2, mean)
  se <- function(v){sd(v)/sqrt(length(v))}
  pam.df.tmp.sd <- apply(pam.df.tmp, 2, sd)
  
  gg.df <- data.frame(Time=new.samps$time_nr, Tree=new.samps$tree_ID, Mean=pam.df.tmp.avg, SD=pam.df.tmp.sd)
  
  plot(0, type="n", main=paste0("Clustering into ", i, " groups \n Cluster ", c) ,xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min.pam.df, max.pam.df), xlab="Time", ylab="Normalized Gene Expression", col.main=colorRampPalette(brewer.pal(12,"Set3"))(i)[c], xaxt = "n")
  axis(side=1, at=month.days[,2], labels=month.days[,1])
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[246])
  
  for(t in trees.order$legend[trees.order$legend %in% unique(new.samps$tree_legend)]){
    lines(new.samps$time_nr[new.samps$tree_legend == t], pam.df.tmp.avg[ new.samps$tree_legend == t] , col=trees.order$color[trees.order$legend==t], type="b", pch=trees.order$pch[trees.order$legend==t], cex=trees.order$cex[trees.order$legend==t]) 
    plot.error.bars(x=new.samps$time_nr[new.samps$tree_legend == t], y=pam.df.tmp.avg[new.samps$tree_legend == t], delta=pam.df.tmp.sd[new.samps$tree_legend == t], col=trees.order$color[trees.order$legend==t], barwidth=100000)
    
  }
  legend("topleft", legend = trees.order$legend[trees.order$legend %in% unique(new.samps$tree_legend)], col=trees.order$color[trees.order$legend %in% unique(new.samps$tree_legend)], cex=0.5, text.col=trees.order$color[trees.order$legend %in% unique(new.samps$tree_legend)])  
}
dev.off()
    


pdf(paste0(out.path,"/",out.name,"_Avg_Expression_plusWP_in_cluster" ,i, ".pdf"), w=10, h=5) 
for(c in 1:i){
  # c=1
  pam.df.tmp <- pam.df[pam.cl==c, -ncol(pam.df)]
  pam.df.tmp <-  pam.df.tmp[, rownames(new.samps)]
  
  pam.df.tmp.avg <- apply(pam.df.tmp, 2, mean)
  se <- function(v){sd(v)/sqrt(length(v))}
  pam.df.tmp.sd <- apply(pam.df.tmp, 2, sd)
  
  gg.df <- data.frame(Time=new.samps$time_nr, Tree=new.samps$tree_ID, Mean=pam.df.tmp.avg, SD=pam.df.tmp.sd)

  WP.org <- new.samps$Water.Potential
  WP.scaled <- normalize.counts(counts=t(as.matrix(WP.org)), norm.method="range", c=min.pam.df, d=max.pam.df)
  WP.at <- seq(0, 1, by=0.1)
  WP.at <- WP.at[WP.at <= max(WP.org,na.rm=T ) & WP.at >= min(WP.org, na.rm=T)]
  WP.at.scaled <- normalize.counts(counts=t(as.matrix(c(WP.org, WP.at))), norm.method="range", c=min.pam.df, d=max.pam.df)
  WP.at.scaled <- WP.at.scaled[-(1:length(WP.org))]
  
  #par()$mar
  par(mar=c(5.1,4.1,4.1,5.1))
      
  plot(0, type="n", main=paste0("Clustering into ", i, " groups \n Cluster ", c) ,xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min.pam.df, max.pam.df), xlab="Time", ylab="Normalized Gene Expression", col.main=colorRampPalette(brewer.pal(12,"Set3"))(i)[c], xaxt = "n", las=1)
  axis(side=1, at=month.days[,2], labels=month.days[,1])
  axis(side=4, at=WP.at.scaled, labels=WP.at, las=1)
  mtext("Water potential", side=4, line=3, las=0)
#   mtext("Water potential", side=4, adj=0, outer = TRUE,  line=3)
  
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[246])
  
  for(t in trees.order$legend[trees.order$legend %in% unique(new.samps$tree_legend)]){
    lines(new.samps$time_nr[new.samps$tree_legend == t], pam.df.tmp.avg[ new.samps$tree_legend == t] , col=trees.order$color[trees.order$legend==t], type="b", pch=trees.order$pch[trees.order$legend==t], cex=1.5*trees.order$cex[trees.order$legend==t]) 

    # WP
    lines(new.samps$time_nr[new.samps$tree_legend == t],  WP.scaled[ new.samps$tree_legend == t]  , col=trees.order$color[trees.order$legend==t], lwd=0.5, lty=2) 
    
  }
  legend("topleft", legend = trees.order$legend[trees.order$legend %in% unique(new.samps$tree_legend)], col=trees.order$color[trees.order$legend %in% unique(new.samps$tree_legend)], cex=0.5, text.col=trees.order$color[trees.order$legend %in% unique(new.samps$tree_legend)])  
}
dev.off()


       
    
    # genes.flowering <- genes.full.description[!is.na(genes.full.description$Flowering), ]
    
    
    #     pdf(paste0(out.path,"/",out.name,"_Expression_Flowering_genes_k" ,i, ".pdf"), w=10, h=5)
    #     for(c in 1:i){
    #       # c=5
    #       pam.df.tmp <- pam.df[pam.cl==c, -ncol(pam.df)]
    #       pam.df.tmp <-  pam.df.tmp[, rownames(new.samps)]
    #       
    #       pam.df.tmp.fl <- merge(cbind(rownames(pam.df.tmp), pam.df.tmp), genes.flowering, by=1, all=FALSE)
    #       
    #       if(nrow(pam.df.tmp.fl) >= 1){
    #         for(j in 1:nrow(pam.df.tmp.fl)){
    #           # j=1
    #           expr <- pam.df.tmp.fl[j, rownames(new.samps)]
    #           
    #           plot(0, type="n", main=paste("Clustering into ", i, " groups: Cluster ", c, "\n", pam.df.tmp.fl[j,1], pam.df.tmp.fl[j,"AT_ID"], "\n", pam.df.tmp.fl[j,"Description2"]) ,xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min.pam.df, max.pam.df), xlab="Time", ylab="Normalized Gene Expression", col.main=colorRampPalette(brewer.pal(12,"Set3"))(i)[c])
    #           rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[246])
    #           for(t in levels(new.samps$tree_ID)){
    #             lines(new.samps$time_nr[new.samps$tree_ID==t], expr[new.samps$tree_ID==t] , col=which(levels(new.samps$tree_ID)==t), type="b", pch=ifelse(new.samps$drough.control[new.samps$tree_ID==t]=="drought", 18, 16))    
    #           }
    #           legend("topleft", legend = levels(new.samps$tree_ID), col=1:length(levels(new.samps$tree_ID)), cex=1, text.col= 1:length(levels(new.samps$tree_ID)), pch=c(18, 16, 16, 16, 18, 18))
    #           
    #         }
    #       }
    #     } 
    #     dev.off()
    
    return(pam.obj)
    
  }, mc.cores=1)
  
  
  save(pam.obj.all, file=paste0(out.path,"/",out.name,"_pam_obj.Rdata"))
  
  cat("Calculate GAP indexes... \n")
  
  Gap.inds <- mclapply(1:(length(pam.obj.all)-1), function(p){ 
    # p=2
    cat("... gap index for", p+1, "...\n")
    clall=cbind(pam.obj.all[[p]]$cluster, pam.obj.all[[p+1]]$cluster)
    pam.indexGap.pc  <- index.Gap(x=pam.x, clall=clall , method="k-means", centrotypes="centroids", reference.distribution="pc", B=100)$diffu
    pam.indexGap.unif <- index.Gap(x=pam.x, clall=clall, method="k-means", centrotypes="centroids", reference.distribution="unif", B=100)$diffu
    
    return(c(pam.indexGap.pc, pam.indexGap.unif))
    
  }, mc.cores=mc.cores)
  
  
#   pam.indexG1 <- unlist(lapply(pam.obj.all, function(p) p$index.G1))
#   pam.indexS <- unlist(lapply(pam.obj.all, function(p) p$index.S))
  pam.indexGap.pc <- unlist(lapply(Gap.inds, function(p) p[1]))
  pam.indexGap.unif <- unlist(lapply(Gap.inds, function(p) p[2]))

names(pam.indexGap.pc) <- prior.nr.cl[-length(prior.nr.cl)]
names(pam.indexGap.unif) <- prior.nr.cl[-length(prior.nr.cl)]
# nr.cl.Gap.pc <- names(which.min(which(pam.indexGap.pc >= 0)))
# nr.cl.Gap.unif <- names(which.min(which(pam.indexGap.unif >= 0)))


#   save(pam.obj.all, pam.indexG1, pam.indexS, pam.indexGap.pc, pam.indexGap.unif, file=paste0(out.path,"/",out.name,"_pam_obj.Rdata"))

save(pam.obj.all, pam.indexGap.pc, pam.indexGap.unif, file=paste0(out.path,"/",out.name,"_pam_obj.Rdata"))
  
  pdf(paste0(out.path,"/",out.name,"_kmeans_idexes.pdf"))
#   plot(prior.nr.cl, pam.indexG1, main="Calinski-Harabasz index - max", xlab="nr of clusters")
#   plot(prior.nr.cl, pam.indexS, main="Silhouette index - max", xlab="nr of clusters")
plot(prior.nr.cl[-length(prior.nr.cl)], pam.indexGap.pc, main= paste0("Gap pc statistic \n final nr of clusters = min{nc: diffu(nc, nc+1) >= 0} \n "),  xlab="nc", ylab="diffu", pch=19, col=ifelse(pam.indexGap.pc >= 0, 2,1))
abline(h=0)
abline(v=prior.nr.cl[-length(prior.nr.cl)], lty=2, lwd=0.5, col=ifelse(pam.indexGap.pc >= 0, 2,"grey"))
plot(prior.nr.cl[-length(prior.nr.cl)], pam.indexGap.unif, main=paste0("Gap unif statistic \n final nr of clusters = min{nc: diffu(nc, nc+1) >= 0} \n "), xlab="nc", ylab="diffu", pch=19, col=ifelse(pam.indexGap.unif >= 0, 2,1))
abline(h=0)
abline(v=prior.nr.cl[-length(prior.nr.cl)], lty=2, lwd=0.5, col=ifelse(pam.indexGap.unif >= 0, 2,"grey"))
  dev.off()
  
  
  pam.cls <- do.call("cbind", lapply(pam.obj.all, function(p) p$cluster))
  colnames(pam.cls) <- paste0("clustering", prior.nr.cl)
  pam.df <- cbind(ID=rownames(pam.df), pam.cls,  pam.df)
  
  pam.df <- merge(genes.full.description, pam.df, by="ID", all.y=T)
  
  write.table(pam.df, paste0(out.path,"/",out.name,"_kmeans_table.xls" ), sep="\t", row.names=F, quote=F)
write.table(pam.df, paste0(out.path,"/",out.name,"_kmeans_table.txt" ), sep="\t", row.names=F, quote=F)

  
}

