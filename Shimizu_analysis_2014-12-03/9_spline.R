#####################################################################################################
# BioC 3.0

### Fit splines 

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


trees.order <- trees.order[!grepl(pattern = "990|8212", trees.order$tree_ID), ]
trees.order <- trees.order[trees.order$legend != "8266-drought-flower_bud", ]

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



d.cpm <- dcpm.norm


genes <- rownames(d.cpm)[1:20]


dir.create("Plots_Splines/", showWarnings=F, recursive=T)



pdf(paste0("Plots_Splines/" , "Spline_fitting" ,".pdf"), h=5, w=10)


for(j in 1:length(genes)){
  # j=1
  
  plot(0, type="n", main=paste0(genes[j]) ,xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(d.cpm[genes[j], ])), max(na.omit(d.cpm[genes[j], ]))), xlab="Time", ylab="Gene Expression", xaxt = "n")
  axis(side=1, at=month.days[,2], labels=month.days[,1])
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[246])
  for(t in trees.order$legend){
    # t=trees.order$legend[1]
    
    time <- new.samps$time_nr[new.samps$tree_legend == t]
    expr <- d.cpm[genes[j], new.samps$tree_legend == t]
    
    lines(time, expr , col=trees.order$color[trees.order$legend==t], type="b", pch=trees.order$pch[trees.order$legend==t], cex=trees.order$cex[trees.order$legend==t]/2, lwd=1, lty = 3) 
    
    
    sm.spl <- smooth.spline(time, expr, spar=0.4)
    lines(sm.spl, col=trees.order$color[trees.order$legend==t], type="l", pch=trees.order$pch[trees.order$legend==t], cex=trees.order$cex[trees.order$legend==t], lwd=4) 
    
    
  }
  legend("topright", legend = trees.order$legend, col=trees.order$color, cex=0.5, text.col=trees.order$color)
  
}


dev.off()

































