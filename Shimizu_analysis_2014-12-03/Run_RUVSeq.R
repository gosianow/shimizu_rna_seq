
library(ggplot2)
library(limma)

run.RUVSeq <- function(x, new.samps, neg.contrl.ID, elim.samps, plots.path="Plots_RUN_RUVSeq_normalization", plot.name=""){ 
  
  
  dir.create(plots.path, recursive=T, showWarnings = FALSE)
  

  x.tmp <- x[, new.samps$sample_name]  
  x.tmp <- x.tmp[,!names(x.tmp) %in% elim.samps]
  new.samps.tmp <- new.samps[!new.samps$sample_name %in% elim.samps, ]
  
  
  d <- DGEList(x.tmp, group=new.samps.tmp$tree_ID)
  d <- calcNormFactors(d)
  
  
  # make sure a gene is expressed (CPM > 1) in more samples
  cps <- cpm(d, normalized.lib.sizes=TRUE)
  d <- d[ rowSums(cps>1) >  10 , ] # ncol(cps)/2
  
  
  # normalization with betweenLaneNormalization
  set.org <- newSeqExpressionSet(as.matrix(d$counts), phenoData = data.frame(tree=as.factor(new.samps.tmp$tree_ID), row.names = colnames(x.tmp))) 
  
  set <- betweenLaneNormalization(set.org, which = "upper")
  
  # RUVg
cIdx <- intersect(neg.contrl.ID, rownames(counts(set)))
  set1 <- RUVg(set, cIdx=cIdx, k = 1, round=TRUE) 


pdf( paste(plots.path, "/",  plot.name, "_W", ".pdf", sep=""))
barplot(pData(set1)$W_1, beside=TRUE, col=new.samps.tmp$tree_col)
dev.off()



pdf( paste(plots.path, "/", plot.name, "_RLE_NegativeControl", ".pdf", sep=""), width = 15, height = 7)

plotRLE(counts(set.org)[cIdx,], outline = FALSE, col = new.samps.tmp$tree_col, main="No normalisation")  
plotRLE(normCounts(set)[cIdx,], outline = FALSE, col = new.samps.tmp$tree_col, main="Upperquantile normalisation")
plotRLE(normCounts(set1)[cIdx,], outline = FALSE, col = new.samps.tmp$tree_col, main="RUVg normalisation")

dev.off()




pdf( paste(plots.path, "/", plot.name, "_MDS_pc2", ".pdf", sep=""), width = 10, height = 10)


mds <- plotMDS(DGEList(counts(set.org), lib.size=rep(1, ncol(counts(set.org)) )) , col=new.samps.tmp$tree_col, top=500, labels=new.samps.tmp$short.name, prior.count=2, main="No normalisation")
mds <- plotMDS(DGEList( normCounts(set) , lib.size=rep(1, ncol( normCounts(set) ))) , col=new.samps.tmp$tree_col, top=500, labels=new.samps.tmp$short.name, prior.count=2, main="Upperquantile normalisation")
mds <- plotMDS(DGEList( normCounts(set1) , lib.size=rep(1, ncol( normCounts(set1) ))) , col=new.samps.tmp$tree_col, top=500, labels=new.samps.tmp$short.name, prior.count=2, main="RUVg normalisation")
mds <- plotMDS(d, col=new.samps.tmp$tree_col, top=500, labels=new.samps.tmp$short.name, prior.count=2, main="d")

dev.off()




pdf( paste(plots.path, "/",plot.name, "_MDS_pc10", ".pdf", sep=""), width = 10, height = 10)


plotMDS(DGEList(counts(set.org), lib.size=rep(1, ncol(counts(set.org)) )) , col=new.samps.tmp$tree_col, top=500, labels=new.samps.tmp$short.name, prior.count=10, main="No normalisation")
plotMDS(DGEList( normCounts(set) , lib.size=rep(1, ncol( normCounts(set) ))) , col=new.samps.tmp$tree_col, top=500, labels=new.samps.tmp$short.name, prior.count=10, main="Upperquantile normalisation")
plotMDS(DGEList( normCounts(set1) , lib.size=rep(1, ncol( normCounts(set1) ))) , col=new.samps.tmp$tree_col, top=500, labels=new.samps.tmp$short.name, prior.count=10, main="RUVg normalisation")
plotMDS(d, col=new.samps.tmp$tree_col, top=500, labels=new.samps.tmp$short.name, prior.count=10, main="d")

dev.off()




pdf( paste(plots.path, "/",  plot.name, "_MDS_bcv", ".pdf", sep=""), width = 10, height = 10)


plotMDS(DGEList(counts(set.org), lib.size=rep(1, ncol(counts(set.org)) )) , col=new.samps.tmp$tree_col, top=500, labels=new.samps.tmp$short.name, method="bcv", main="No normalisation")
plotMDS(DGEList( normCounts(set) , lib.size=rep(1, ncol( normCounts(set) ))) , col=new.samps.tmp$tree_col, top=500, labels=new.samps.tmp$short.name, method="bcv", main="Upperquantile normalisation")
plotMDS(DGEList( normCounts(set1) , lib.size=rep(1, ncol( normCounts(set1) ))) , col=new.samps.tmp$tree_col, top=500, labels=new.samps.tmp$short.name, method="bcv", main="RUVg normalisation")
plotMDS(d, col=new.samps.tmp$tree_col, top=500, labels=new.samps.tmp$short.name, method="bcv", main="d")

dev.off()


  
#   pdf( paste(plots.path, "/", plot.name, "_MV_cpm_filtering", ".pdf", sep=""))  
#   meanVarPlot(set.org, log=T)
#   meanVarPlot(set, log=T)
#   meanVarPlot(set1, log=T) 
#   dev.off()

  

  pdf( paste(plots.path, "/",  plot.name, "_RLE", ".pdf", sep=""), width = 10, height = 5)
  
  plotRLE(set.org, outline = FALSE, col = new.samps.tmp$tree_col, main="No normalisation") 
  plotRLE(set, outline = FALSE, col = new.samps.tmp$tree_col, main="Upperquantile normalisation")
  plotRLE(set1, outline = FALSE, col = new.samps.tmp$tree_col, main="RUVg normalisation")
  
  dev.off()
  
  

  pdf( paste(plots.path, "/",plot.name, "_PCA", ".pdf", sep=""), width = 10, height = 10)
  
  plotPCA(set.org, col = new.samps.tmp$tree_col, main="No normalisation")
  plotPCA(set, col = new.samps.tmp$tree_col, main="Upperquantile normalisation")
  plotPCA(set1, col = new.samps.tmp$tree_col, main="RUVg normalisation")
  
  plotPCA(set.org, col = new.samps.tmp$colors.wp, main="No normalisation")
  plotPCA(set, col = new.samps.tmp$colors.wp, main="Upperquantile normalisation")
  plotPCA(set1, col = new.samps.tmp$colors.wp, main="RUVg normalisation")

  dev.off()
    
return(set1)  

}





run.edgeR.RUVSeq <- function(set1, new.samps, model.formula, varialbs, elim.samps, FDR=0.1, plots.path="Plots_RUN_edgeR_RUVSeq", plot.name="", varialbs.plot, fit.clr = 1, LRTcoef=2, trees.order=trees.order, genes.full.description=genes.full.description, AT.id=AT.id){ 
  
  
  dir.create(plots.path, recursive=T, showWarnings = FALSE)
  
  library(stringr)
  model.char <- str_replace_all(as.character(model.formula), " ", "") # can be gsub()
  model.char <- str_replace_all(model.char[2], "\\.", "_")
  
  new.samps.tmp <- new.samps[row.names(pData(set1)), ]
  new.samps.tmp <- new.samps.tmp[!new.samps.tmp$sample_name %in% elim.samps, ]
  
  
  # DE
  
  d <- DGEList(counts = counts(set1)[, new.samps.tmp$sample_name])
  d <- calcNormFactors(d, method = "upperquartile")
  
  
  # design model matrix
  design <- model.matrix(model.formula, data=merge(new.samps.tmp, pData(set1), by=0, sort=FALSE) )
  print(head(design))

  
  
  d <- estimateGLMCommonDisp(d, design)
  d <- estimateGLMTagwiseDisp(d, design) 
  d <- estimateGLMTrendedDisp(d,design)
  fit <- glmFit(d, design)
  lrt <- glmLRT(fit, coef = LRTcoef) 
  
  top.tags <- topTags(lrt, n=nrow(lrt$table))
  top.tags <- top.tags$table
  top.genes <- rownames(top.tags[top.tags$FDR < FDR,])
  length(top.genes)
  
  
  pdf( paste(plots.path, "/", model.char, plot.name, "_HISTpvs", ".pdf", sep=""), width = 7, height =7)  
  hist(top.tags$PValue, breaks=50, col="deeppink")
  # plotSmear(d)
  dev.off()
  
  
  if(length(top.genes)==0){
    cat("* NO genes found significant \n")
    invisible(top.tags)
  }
  
  cat(paste0("* ",length(top.genes)," genes found significant \n"))
  
  
  
  write.table( merge(fit$coefficients[top.genes,], top.tags[top.genes,], by=0, sort=FALSE) ,paste(plots.path, "/", model.char,  plot.name, "_top_fitting",".xls", sep=""), sep="\t", row.names=F, quote=F)
  
  
  Top.fitting.list <- list()
  UP.coeffs <- matrix(NA, length(top.genes), ncol(design))
  rownames(UP.coeffs) <- top.genes  
  DOWN.coeffs <- matrix(NA, length(top.genes), ncol(design))
  rownames(DOWN.coeffs) <- top.genes
  
rownames(fit$offset) <- rownames(fit$counts)  


  pdf( paste(plots.path, "/", model.char, plot.name, "_top_fitting", ".pdf", sep=""), width = 10, height =7.5)
  par(mfrow=c(2,3))
  
  for(i in top.genes){
    # i=top.genes[1];i 
    
    
    # log.cpm.raw <- as.numeric(log(d.cpm[i, , drop=FALSE])) - fit$offset[i,]
    log.cpm.raw <- as.numeric(log(d$counts[i, , drop=FALSE])) - fit$offset[i,]
    log.cpm.fit <- as.numeric(fit$coefficients[i, , drop=FALSE] %*% t(fit$design)) 
    
    
    plot.data <- data.frame( new.samps.tmp[, varialbs.plot], log.cpm.raw, log.cpm.fit )
    colnames(plot.data) <- c(varialbs.plot, "Log.Counts.in.CPM.RAW", "Log.Counts.in.CPM.FIT")
    rownames(plot.data) <- row.names(fit$samples)
    
    Top.fitting.list[[i]] <- plot.data
    
    coeffs <- round(fit$coefficients[i,], 2)
    
    AT <- AT.id[AT.id[,1] == i, 2]
    if.UP <- if.DOWN <- FALSE
    if(length(AT)!=0){
      if.UP <- AT %in% UP.genes[,1]
      if(if.UP)
        UP.coeffs[i,] <- coeffs
      if.DOWN <- AT %in% DOWN.genes[,1]
      if(if.DOWN)
        DOWN.coeffs[i,] <- coeffs
    }
    
    new.samps.tmp$tree_ID_nr <- new.samps.tmp$tree_ID
    levels(new.samps.tmp$tree_ID_nr) <- 1:length(levels(new.samps.tmp$tree_ID)) +1
    
    plot(Top.fitting.list[[i]][,1], Top.fitting.list[[i]][,2], main=paste(i, "/", AT, "\n", "Coefffs:",paste(coeffs, collapse=", "), "\n FDR:", top.tags[i, "FDR"], "\n UP:", if.UP, "DOWN:", if.DOWN), xlab = varialbs.plot, ylab="Log(Counts in CPM)", cex.main=0.7, pch=ifelse(new.samps.tmp$drough.control=="drought", 5, 1), col=new.samps.tmp$tree_col)
    points(Top.fitting.list[[i]][,1], Top.fitting.list[[i]][,3], col=1, pch="*", cex=2)
    
  }
  
  dev.off()
  
  
#   pdf(paste(plots.path, "/", model.char, plot.name, "_top_genes", ".pdf", sep=""), width = 10, height =5)
#   
#   for(j in 1:length(top.genes)){
#     # j=1
#     
#     plot(0, type="n", main=paste0(top.genes[j]) ,xlim=c(min(new.samps.tmp$time_nr), max(new.samps.tmp$time_nr)), ylim=c(min(na.omit(d.cpm[top.genes[j], ])), max(na.omit(d.cpm[top.genes[j], ]))), xlab="Time", ylab="Gene Expression in cpm", xaxt = "n")
#     axis(side=1, at=month.days[,2], labels=month.days[,1])
#     rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[246])
#     for(t in trees.order$legend){
#       # t=trees.order$legend[1]
#       lines(new.samps.tmp$time_nr[new.samps.tmp$tree_legend == t], d.cpm[top.genes[j], new.samps.tmp$tree_legend == t] , col=trees.order$color[trees.order$legend==t], type="b", pch=trees.order$pch[trees.order$legend==t], cex=trees.order$cex[trees.order$legend==t]) 
#     }
#     legend("topright", legend = trees.order$legend, col=trees.order$color, cex=0.5, text.col=trees.order$color)
#     
#   }
#   
#   dev.off()
  
invisible(top.tags)

}






