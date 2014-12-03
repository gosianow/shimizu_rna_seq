############################################
### help functions:
############################################



plot.mean.sd <- function(d){ 
  d.mean <- apply(d, 1, mean)
  d.sd <- apply(d, 1, sd)
  plot(d.mean, d.sd)
}

normalize.counts <- function(counts, norm.method=c("norm", "01")[1]){
  
  t(apply(counts, 1, function(g){
    
    if(norm.method=="norm"){
      m <- mean(g)
      sd <- sd(g) 
      return((g-m)/sd)
    } else if(norm.method=="01"){
      return((g-min(g))/(max(g)-min(g)))
    }
  }))
}

plot.error.bars <- function(x, y, delta, col=1, barwidth=100){
  segments(x,y-delta,x,y+delta, col=col)
  segments(x-barwidth,y+delta,x+barwidth,y+delta, col=col)
  segments(x-barwidth,y-delta,x+barwidth,y-delta, col=col)
}


############################################
### PAM.clustering
############################################


# pam.x <- count.data$d.cpm.org
# pam.x <- pam.x[1:200,]
# 
# dist.method= "euclidean"
# norm.method="norm"
# out.path="Clustering/PAM/"
# out.name=""
# prior.nr.cl=2:5
# mc.cores=1
# ylim=c(-2,6)
# clustering.samps=NA
# colors=c("yellow","blue")



PAM.clustering <- function(pam.x, new.samps, genes.full.description,  dist.method= "euclidean", norm.method="norm", out.path="Clustering/PAM/", out.name="", prior.nr.cl=2:5, mc.cores=1, ylim=c(-2,6), clustering.samps=NA, colors=c("yellow","blue") ){
  
  dir.create(out.path, showWarnings=FALSE, recursive=TRUE)
  
  samps.order <- row.names(new.samps[order( factor(new.samps$developmental_stage, levels=as.character(sort(unique(new.samps$developmental_stage), decreasing = T))), new.samps$drough.control, new.samps$tree_ID, new.samps$time_nr), ])
  pam.x <- pam.x[, samps.order]
  new.samps <- new.samps[samps.order,]
  
  # !!! correct clustering.samps thingy

  if(!is.na(clustering.samps))
    new.samps <- new.samps[clustering.samps,]
  
  if(norm.method!="none")
    pam.x <- normalize.counts(counts=pam.x, norm.method=norm.method)
  
  pam.x <- na.omit(pam.x)
  
  pam.df <- as.data.frame(pam.x)
  
  if(!is.na(clustering.samps))
    pam.df <- as.data.frame(pam.x[, clustering.samps])
  
  cat("Calculating distance matrix...\n")
  
  if(dist.method %in% c("spearman" ,"kendall", "pearson") ){
    
    sim.obj.g <- cor(t(pam.x), method=dist.method)
    pam.dist <- as.dist(1-sim.obj.g)
    
  } else {
    
    pam.dist <- dist(x=pam.x , method = dist.method)
  }
  
  
  pam.obj.all <- mclapply(prior.nr.cl, function(i){
    # i=2
    cat("Clustering into ", i, "clusters \n")
    pam.obj <- pam(x=pam.dist, k=i) 
    pam.cl <- pam.obj$clustering
    
    #     pam.obj$index.S <- index.S(d=pam.dist, cl=pam.cl)
    #     pam.obj$index.G1 <- index.G1(x=pam.x, d=pam.dist, cl=pam.cl, centrotypes="medoids")
    #     
    pam.df[, paste0("clustering", i)] <- pam.cl
    
    
    pdf(paste0(out.path, "/",out.name,"_Sillhouette_plot_k" ,i, ".pdf"), w=15, h=10)
    pam.sill <- silhouette(x=pam.cl, dist=pam.dist)
    
    summ.pam.sill <- summary(pam.sill)
    pam.obj$index.S <- summ.pam.sill$avg.width
    
    plot(pam.sill, col=colorRampPalette(brewer.pal(12,"Set3"))(i))
    dev.off()
    
    pam.df.sort <- pam.df[order(pam.df[, paste0("clustering", i)]), ]
    
    png(paste0(out.path,"/",out.name,"_Heatmap_plot_k" ,i, ".png"), w=1000, h=1000)
    heatmap.2axis(x=as.matrix(pam.df.sort[,-ncol(pam.df.sort)]), Rowv=NA, Colv=NA, dendrogram="none", labRow=NA, labCol=NULL, scale="none", trace="none", col=colorRampPalette(colors)(250), key=T, keysize =c(0.5, 0.5), density.info="none", margins = c(5,5), xlab = "Samples", ylab = "Clusters of DE Genes", RowSideColors=colorRampPalette(brewer.pal(12,"Set3"))(i)[pam.df.sort[, ncol(pam.df.sort)]], ColSideColors=new.samps$tree_col) 
    dev.off()
    
    
    if(length(ylim)==2){
      min.pam.df <- ylim[1]
      max.pam.df <- ylim[2]   
    } else {
      min.pam.df <- min(pam.df[, -ncol(pam.df)])
      max.pam.df <- max(pam.df[, -ncol(pam.df)])
    }
    
    
    pdf(paste0(out.path,"/",out.name,"_Expression_k" ,i, ".pdf"), w=10, h=5) 
    for(c in 1:i){
      # c=1
      pam.df.tmp <- pam.df[pam.cl==c, -ncol(pam.df)]
      pam.df.tmp <-  pam.df.tmp[, rownames(new.samps)]
      
      pam.df.tmp.avg <- apply(pam.df.tmp, 2, mean)
      se <- function(v){sd(v)/sqrt(length(v))}
      pam.df.tmp.sd <- apply(pam.df.tmp, 2, sd)
      
      
      gg.df <- data.frame(Time=new.samps$time_nr, Tree=new.samps$tree_ID, Mean=pam.df.tmp.avg, SD=pam.df.tmp.sd)
      
      #         pd <- position_dodge(.1) # move them .05 to the left and right
      #       print(ggplot(gg.df, aes(x=Time, y=Mean, colour=Tree, group=Tree)) + geom_errorbar(aes(ymax=Mean+SD, ymin=Mean-SD), width=.1, position=pd) + geom_line(position=pd) + geom_point(position=pd, size=3) + xlab("Time") +  ylab("Normalized gene expression") + ggtitle(paste0("Clustering into ", i, " groups \n Cluster ", c)) + theme(legend.justification=c(1,1), legend.position="right"))
      
      plot(0, type="n", main=paste0("Clustering into ", i, " groups \n Cluster ", c) ,xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min.pam.df, max.pam.df), xlab="Time", ylab="Normalized Gene Expression", col.main=colorRampPalette(brewer.pal(12,"Set3"))(i)[c])
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[246])
      for(t in levels(new.samps$tree_ID)){
        lines(new.samps$time_nr[new.samps$tree_ID==t], pam.df.tmp.avg[new.samps$tree_ID==t] , col=which(levels(new.samps$tree_ID)==t), type="b", pch=ifelse(new.samps$drough.control[new.samps$tree_ID==t]=="drought", 18, 16)) 
        plot.error.bars(x=new.samps$time_nr[new.samps$tree_ID==t], y=pam.df.tmp.avg[new.samps$tree_ID==t], delta=pam.df.tmp.sd[new.samps$tree_ID==t], col=which(levels(new.samps$tree_ID)==t), barwidth=100000)
        
      }
      legend("topleft", legend = levels(new.samps$tree_ID), col=1:length(levels(new.samps$tree_ID)), cex=1, text.col= 1:length(levels(new.samps$tree_ID)), pch=c(18, 16, 16, 16, 18, 18))  
    }
    dev.off()
    
    
    genes.flowering <- genes.full.description[!is.na(genes.full.description$Flowering), ]
    
    
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
    #       
    #       
    #       
    #     } 
    #     dev.off()
    
    
    return(pam.obj)
    
  }, mc.cores=1)
  
  
  #   save(pam.obj.all, file=paste0(out.path,"/",out.name,"_pam_obj.Rdata"))
  #   
  #   cat("Calculate indexes... \n")
  #   
  #   Gap.inds <- mclapply(1:(length(pam.obj.all)-1), function(p){ 
  #     # p=2
  #     clall=cbind(pam.obj.all[[p]]$clustering, pam.obj.all[[p+1]]$clustering)
  #     pam.indexGap.pc  <- index.Gap(x=pam.x, d=pam.dist, clall=clall , method="pam", centrotypes="medoids", reference.distribution="pc")$diffu
  #     pam.indexGap.unif <- index.Gap(x=pam.x, d=pam.dist, clall=clall, method="pam", centrotypes="medoids", reference.distribution="unif")$diffu
  #     
  #     return(c(pam.indexGap.pc, pam.indexGap.unif))
  #     
  #   }, mc.cores=mc.cores)
  
  
  pam.indexS <- unlist(lapply(pam.obj.all, function(p) p$index.S))
  #   pam.indexG1 <- unlist(lapply(pam.obj.all, function(p) p$index.G1))
  #   pam.indexGap.pc <- unlist(lapply(Gap.inds, function(p) p[1]))
  #   pam.indexGap.unif <- unlist(lapply(Gap.inds, function(p) p[2]))
  #   
  #   save(pam.obj.all, pam.indexS, pam.indexG1, pam.indexGap.pc, pam.indexGap.unif, file=paste0(out.path,"/",out.name,"_pam_obj.Rdata"))
  
  pdf(paste0(out.path,"/",out.name,"_pam_idexes.pdf"))
  plot(prior.nr.cl, pam.indexS, main="Silhouette index - max", xlab="nr of clusters")
  #   plot(prior.nr.cl, pam.indexG1, main="Calinski-Harabasz index - max", xlab="nr of clusters")
  #   plot(prior.nr.cl[-length(prior.nr.cl)], pam.indexGap.pc, main="Gap pc statistic - max", xlab="nr of clusters")
  #   plot(prior.nr.cl[-length(prior.nr.cl)], pam.indexGap.unif, main="Gap unif statistic - max", xlab="nr of clusters")
  dev.off()
  
  pam.cls <- do.call("cbind", lapply(pam.obj.all, function(p) p$clustering))
  colnames(pam.cls) <- paste0("clustering", prior.nr.cl)
  pam.df <- cbind(ID=rownames(pam.df), pam.cls,  pam.df)
  pam.df <- merge(genes.full.description, pam.df, by="ID", all.y=T)
  
  write.table(pam.df, paste0(out.path,"/",out.name,"_pam_table.xls" ), sep="\t", row.names=F, quote=F)
  
}

