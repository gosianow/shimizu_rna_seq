# model.formula = ~ Water.Potential
# varialbs=c("Water.Potential")
# varialbs.plot = "Water.Potential"

# plots.path="Plots_RUN_edgeR"
# fit.clr = 1
# LRTcoef=2
# FDR=0.1
# 
# elim.samps = c(flowered.samps, control.samps)
# plot.name="_drought_samps"



run.edgeR <- function(x, new.samps, model.formula, varialbs, elim.samps, FDR=0.1, plots.path="Plots_RUN_edgeR", plot.name="", varialbs.plot, fit.clr = 1, LRTcoef=2, trees.order=trees.order, genes.full.description=genes.full.description, AT.id=AT.id){ 
  
  dir.create(plots.path, recursive=T, showWarnings = FALSE)
  
  library(stringr)
  model.char <- str_replace_all(as.character(model.formula), " ", "") # can be gsub()
  model.char <- str_replace_all(model.char[2], "\\.", "_")
  
  
  x.tmp <- x[, new.samps$sample_name]  
  x.tmp <- x.tmp[,!names(x.tmp) %in% elim.samps]
  new.samps.tmp <- new.samps[!new.samps$sample_name %in% elim.samps, ]
  
  #   for(i in varialbs){
  #     x.tmp <- x.tmp[,!is.na(new.samps.tmp[, i])]
  #     new.samps.tmp <- new.samps.tmp[!is.na(new.samps.tmp[, i]), ]
  #   } 
  
  
  d <- DGEList(x.tmp, group=new.samps.tmp$tree_ID)
  d <- calcNormFactors(d)
  
  ### make sure a gene is expressed (CPM > 1) in more than 2 samples
  d.cpm <- cpm(d, normalized.lib.sizes=TRUE)
  
#   eps <- min(d.cpm[d.cpm!=0])
  
#   new.samps.tmp$FT_gene <- log(d.cpm["GID031739_3527284",]+eps)
#   new.samps.tmp$SVP_gene <- log(d.cpm["GID037469_3529277",]+eps)
#   new.samps.tmp$FUL_gene <- log(d.cpm["GID056430_3537430",]+eps)
  
  d <- d[ rowSums(d.cpm>1) > ncol(d.cpm)/2, ]
  # cat(paste("*", nrow(d), "genes expressed in more than two samples \n"))
  
  
  # design model matrix
  design <- model.matrix(model.formula, data=new.samps.tmp)
  print(head(design))
  
  # estimate dispersion  
  d <- estimateGLMCommonDisp(d,design)
  d <- estimateGLMTrendedDisp(d,design)
  d <- estimateGLMTagwiseDisp(d,design)
  
  pdf( paste(plots.path, "/", model.char, plot.name, "_MV", ".pdf", sep=""), width = 7, height =7)  
  plotMeanVar(d)
  plotBCV(d)
  dev.off()
  
  # glmFit
  fit <- glmFit(d,design)
  # glmLRT
  lrt <- glmLRT(fit, coef=LRTcoef)
  
  ## plot TOP genes // expression vs meteo var
  top.tags <- topTags(lrt, n=nrow(lrt$table))
  top.tags <- top.tags$table
  top.genes <- rownames(top.tags[top.tags$FDR < FDR,])
  length(top.genes)
  #fit$coefficients[top.genes,]
  
  pdf( paste(plots.path, "/", model.char, plot.name, "_HISTpvs", ".pdf", sep=""), width = 7, height =7)  
  hist(top.tags$PValue, breaks=50, col="orange")
  # plotSmear(d)
  dev.off()
  
  
  if(length(top.genes)==0){
    cat("* NO genes found significant \n")
    invisible(top.tags)
  }
  
  cat(paste0("* ",length(top.genes)," genes found significant \n"))
  
#   pdf( paste(plots.path, "/", model.char, plot.name, "_MA", ".pdf", sep=""), width = 7, height =7)  
#   plotSmear(d, de.tags=top.genes)
#   dev.off()
  
  
  write.table( merge(top.tags[top.genes,] , fit$coefficients[top.genes,], by=0, sort=FALSE), paste(plots.path, "/", model.char,  plot.name, "_top_fitting",".xls", sep=""), sep="\t", row.names=F, quote=F)
  
  
#   d.fit <- DGEList(fit$fitted.values, group=new.samps.tmp$tree_ID[match(names(x.tmp), new.samps.tmp$sample_name)])
#   d.fit <- calcNormFactors(d.fit)
#   d.cpm <- cpm(d, normalized.lib.sizes=TRUE)
#   d.fit.cpm <- cpm(d.fit, normalized.lib.sizes=TRUE)

  
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
    

#     Top.fitting.list[[i]] <- data.frame(new.samps.tmp[, varialbs.plot], log(d.cpm[i,]), log(d.fit.cpm[i,]))
#     colnames(Top.fitting.list[[i]]) <- c(varialbs.plot, "Log.Counts.in.CPM.RAW", "Log.Counts.in.CPM.FIT")
#     rownames(Top.fitting.list[[i]]) <- row.names(fit$samples)
        
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
    
    plot(Top.fitting.list[[i]][,1], Top.fitting.list[[i]][,2], main=paste(i, "/", AT, "\n", "Coefffs:",paste(coeffs, collapse=", "), "\n FDR:", top.tags[i, "FDR"], "\n UP:", if.UP, "DOWN:", if.DOWN), xlab = varialbs.plot, ylab="Log(Counts in CPM)", cex.main=0.7, pch=ifelse(new.samps.tmp$drough.control=="drought", 5, 1), col=new.samps.tmp$tree_col, cex=1.5)
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
#   
  invisible(top.tags)
  
}

