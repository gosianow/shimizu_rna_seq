
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



# setwd("/home/Shared/data/seq/Shimizu_RNA_seq/")
# load("Shimizu_workspace.Rdata")
# 
# plot.data <- d.cpm
# genesID <- rownames(plot.data)[1:10]
# samps.to.plot <- new.samps$sample_name
# out.path <- "Plots_of_selected_genes"
# out.name <- "test"
# plot.WP=TRUE
# ylab="Gene expression in CPM"


plot.genes <- function(plot.data, genesID=rownames(plot.data)[1:10], new.samps, samps.to.plot=new.samps$sample_name, genes.full.description, trees.order, month.days, out.path="Plots", out.name="Expression", plot.WP=TRUE, ylab="Gene expression in CPM"){
  
  dir.create(out.path, showWarnings=F, recursive=T)
  
  
  pdf(paste0(out.path,"/",out.name,".pdf"), w=10, h=5) 
  for(c in 1:length(genesID)){
    # c=1

    genes.description.tmp <- genes.full.description[genes.full.description$ID==genesID[c], , drop=FALSE]
        
    if(nrow(genes.description.tmp)==0){      
      plot.main <- paste0(genesID[c])
    }else{      
      plot.main <- paste0(genesID[c], "  ",ifelse(is.na(genes.description.tmp$AT_ID), "", genes.description.tmp$AT_ID), "\n", ifelse(is.na(genes.description.tmp$AT_symbol), "", genes.description.tmp$AT_symbol), "\n", ifelse(is.na(genes.description.tmp$Description), "", genes.description.tmp$Description))
    }
   
    ylim.min <- min(plot.data[genesID[c],],na.rm=T)
    ylim.max <- max(plot.data[genesID[c],],na.rm=T)
    
    if(all(c(ylim.min, ylim.max)==0)){
      ylim.min <- -1
      ylim.max <- 1
    }
    
    if(plot.WP){
      WP.org <- new.samps$Water.Potential
      WP.scaled <- normalize.counts(counts=t(as.matrix(WP.org)), norm.method="range", c=ylim.min, d=ylim.max)
      WP.at <- seq(0, 1, by=0.1)
      WP.at <- WP.at[WP.at <= max(WP.org,na.rm=T ) & WP.at >= min(WP.org, na.rm=T)]
      WP.at.scaled <- normalize.counts(counts=t(as.matrix(c(WP.org, WP.at))), norm.method="range", c=ylim.min, d=ylim.max)
      WP.at.scaled <- WP.at.scaled[-(1:length(WP.org))]     
    }

    #par()$mar
    par(mar=c(5.1,4.1,4.1,5.1))
    
    plot(0, type="n", main=plot.main ,xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(ylim.min, ylim.max), xlab="Time", ylab=ylab, xaxt = "n", las=1)
    axis(side=1, at=month.days[,2], labels=month.days[,1])
    if(plot.WP){
      axis(side=4, at=WP.at.scaled, labels=WP.at, las=1)
      mtext("Water potential", side=4, line=3, las=0)     
    }

    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[246])
    
    for(t in trees.order$legend[trees.order$legend %in% unique(new.samps$tree_legend)]){
      lines(new.samps$time_nr[new.samps$tree_legend == t], plot.data[ genesID[c], new.samps$tree_legend == t] , col=trees.order$color[trees.order$legend==t], type="b", pch=trees.order$pch[trees.order$legend==t], cex=1.5*trees.order$cex[trees.order$legend==t]) 
      
      if(plot.WP){
        lines(new.samps$time_nr[new.samps$tree_legend == t],  WP.scaled[ new.samps$tree_legend == t]  , col=trees.order$color[trees.order$legend==t], lwd=0.5, lty=2)         
      }
      
    }
    legend("topleft", legend = trees.order$legend[trees.order$legend %in% unique(new.samps$tree_legend)], col=trees.order$color[trees.order$legend %in% unique(new.samps$tree_legend)], cex=0.5, text.col=trees.order$color[trees.order$legend %in% unique(new.samps$tree_legend)])  
  }
  dev.off()
  
  
}



