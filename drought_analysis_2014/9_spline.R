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

out.dir <- "Plots_Splines/"
dir.create(out.dir, showWarnings=F, recursive=T)


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

########### make sure a gene is expressed (CPM > 1) in more than 2 samples
cps <- cpm(d, normalized.lib.sizes=TRUE)
d <- d[ rowSums( cps > 10 ) > 10, ]
# sum(rowSums( cps > 1 ) > 46)
dim(d$counts)

########### make sure FC is higher than 2 (per tree)

dcpm <- cpm(d, normalized.lib.sizes=TRUE)

trees <- unique(new.samps$tree_ID)
FC.tree <- matrix(0, nrow(dcpm), length(trees))
colnames(FC.tree) <- trees
rownames(FC.tree) <- rownames(dcpm)


for(i in 1:length(trees)){
  # i = 1  
  
  dcpm.t <- dcpm[, new.samps$tree_ID == trees[i] ]
  
  fc.all <- apply(dcpm.t, 1, function(g){
    
    if(max(g)  == 0)
      fc <- 0
    else
      fc <- max(g) / min(g[ g != 0 ])
    
    return(fc)
  })
  
  FC.tree[,i] <- fc.all
  
}


d <- d[ rowSums( FC.tree > 3 ) > 0 , ]
# sum( rowSums( FC.tree > 4 ) > 0 )
dim(d$counts)



dcpm <- cpm(d, normalized.lib.sizes=TRUE)


########### log transformation 

d$counts <- d$counts
dcpm <- cpm(d, normalized.lib.sizes=TRUE)
dcpm <- log2(dcpm + 1)


########### normalize to 01

dcpm.norm <- normalize.counts(counts = dcpm, norm.method = "01")

# ### weights if weighted splines
# d1 <- d
# d1$counts <- d1$counts + 1
# weights <- normalize.counts(counts = log(cpm(d1, normalized.lib.sizes=TRUE)), norm.method = "01")

d.cpm <- dcpm.norm

genes <- rownames(d.cpm)

dspl <- vector("list", 4)
dspl <- lapply(c(12, 6, 12, 11), function(n) matrix(0, length(genes), n) )
names(dspl) <-  trees.order$legend 

### standard deviation
sd.spl <- matrix(0, length(genes), 4)
colnames(sd.spl) <- trees.order$legend
### mean expression
m.spl <- matrix(0, length(genes), 4)
colnames(m.spl) <- trees.order$legend


spar=0.1


out.path <- "Plots_Splines/Spline_fitting_01_sp0.1_CPM_predict/"
dir.create(out.path, showWarnings=F, recursive=T)


for(j in 1:length(genes)){
  # j=3
  
  pdf(paste0(out.path, genes[j],".pdf"), h=5, w=10)
  
  plot(0, type="n", main=paste0(genes[j]) ,xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(d.cpm[genes[j], ])), max(na.omit(d.cpm[genes[j], ]))), xlab="Time", ylab="Gene Expression", xaxt = "n")
  axis(side=1, at=month.days[,2], labels=month.days[,1])
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[246])
  
  for(t in trees.order$legend){
    # t=trees.order$legend[2]
    
    time <- new.samps$time_nr[new.samps$tree_legend == t]
    expr <- d.cpm[genes[j], new.samps$tree_legend == t]
    #     w <- weights[genes[j], new.samps$tree_legend == t]
    
    ### plot raw expression
    lines(time, expr , col=trees.order$color[trees.order$legend==t], type="b", pch=trees.order$pch[trees.order$legend==t], cex=trees.order$cex[trees.order$legend==t]/2, lwd=1, lty = 3) 
    
    ### smooth with splines
    
    #     sm.spl <- smooth.spline(time, expr, w = w , spar=spar)
    sm.spl <- smooth.spline(time, expr , spar=spar)
    
    sm.p <- predict(sm.spl, x = sort(all.days$days.nr[all.days$days.nr <= max(time) & all.days$days.nr >= min(time)], decreasing = FALSE))
    
    lines(sm.p, col=trees.order$color[trees.order$legend==t], type="l", pch=trees.order$pch[trees.order$legend==t], cex=trees.order$cex[trees.order$legend==t], lwd=4) 
    
    dspl[[t]][j,] <- sm.spl$y
    sd.spl[j, t] <- sd(sm.spl$y)
    m.spl[j, t] <- mean(sm.spl$y)
    
    if(j == 1)
      colnames(dspl[[t]]) <- paste0(unique(new.samps$short.name[new.samps$tree_legend == t]))  
    
  }
  legend("topright", legend = trees.order$legend, col=trees.order$color, cex=0.5, text.col=trees.order$color)
  
  dev.off()
  
}




dspl <- do.call(cbind, dspl)

rownames(dspl) <- genes
rownames(sd.spl) <- genes
rownames(m.spl) <- genes


new.samps.spl <- unique(new.samps[-which(colnames(new.samps) %in% c("sample_num", "sample_name", "sample_ID"))])
rownames(new.samps.spl) <- new.samps.spl$short.name

new.samps.spl <- new.samps.spl[colnames(dspl),]

new.samps.spl$sample_name <- new.samps.spl$short.name


save(dspl, sd.spl, m.spl, new.samps.spl, FC.tree, file = paste0(out.path , "Spline_fitting_01" ,".RData"))



















