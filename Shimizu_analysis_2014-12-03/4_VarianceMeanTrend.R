#####################################################################################################
# BioC 3.0

### Check the variance mean relationship for the replicate samples (tree 8266 and 970)

#####################################################################################################

### load data

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

out.dir <- "Plots_VarianceMeanTrend/"
dir.create(out.dir, showWarnings=F, recursive=T)

####################################


all(colnames(x)==rownames(new.samps))

dim(x)

# new.samps[,c("sample_name","tree_ID", "time_ch")]

pairs <- list(c("H6_8266_20081208", "H7_8266_20081208"), c("H8_8266_20090105", "H9_8266_20090105"), c("E4_8266_20090202", "F4_8266_20090202"), c("N1_8266_20090224", "N5_8266_20090224"), c("I9_970_20081216", "J1_970_20081216"), c("J2_970_20090207", "K9_970_20090207"))

pairs.samps <- unlist(pairs)


library(edgeR)
d <- DGEList(x[, pairs.samps])
d <- calcNormFactors(d)

# make sure a gene is expressed (CPM > 1) in more samples
cps <- cpm(d, normalized.lib.sizes=TRUE)

# sum(rowSums( cps > 0 ) > 11)
d <- d[ rowSums( cps > 0 ) > 11, ]
dim(d$counts)


dcpm <- cpm(d, normalized.lib.sizes=TRUE)





for(i in 1:length(pairs)){
  # i = 1
  png(paste0(out.dir, "SDMeanScaleLog_",paste0(pairs[[i]], collapse = "_"),".png"), 600, 600)
  
#   v <- abs(dcpm[, pairs[[i]][1]] - dcpm[,pairs[[i]][2]]) 
v <- apply(dcpm[, pairs[[i]]], 1, sd)
  m <- rowMeans(dcpm[, pairs[[i]]])
  
  smoothScatter( log10(m), log10(v) , nrpoints = Inf, nbin = 500, xlab="mean" , ylab="sd", main = paste0(pairs[[i]], collapse = " , "))
  abline(a = 0, b = 1, col = "red")
  
  dev.off()
}


for(i in 1:length(pairs)){
  # i = 1
png(paste0(out.dir, "SDMeanScale50",paste0(pairs[[i]], collapse = "_"),".png"), 600, 600)

#   v <- abs(dcpm[, pairs[[i]][1]] - dcpm[,pairs[[i]][2]]) 
v <- apply(dcpm[, pairs[[i]]], 1, sd)
m <- rowMeans(dcpm[, pairs[[i]]])

l <- lowess(m, v, f = 2/3)

smoothScatter( m, v , nrpoints = Inf, nbin = 500, xlab="mean" , ylab="sd", main = paste0(pairs[[i]], collapse = " , ") , xlim=c(0, 50), ylim=c(0, 50))
abline(a = 0, b = 1, col = "red")
lines(l, col = "blue")

dev.off()

}

#################################

v <- NULL
m <- NULL

for(i in 1:length(pairs)){
  # i = 1
  
  v <- c(v, apply(dcpm[, pairs[[i]]], 1, sd))
  m <- c(m, rowMeans(dcpm[, pairs[[i]]]))
  
}


SDMean.lowess <- l <- lowess(m[m >=1 ],v[m >=1 ])
# s <- smooth.spline(m[m >=1 ], v[m >=1 ], spar=2/3)
s <- smooth.spline(l$x, l$y, spar=2/3)


png(paste0(out.dir, "SDMeanScale50.png"), 600, 600)

smoothScatter( m, v , nrpoints = Inf, nbin = 500, xlab="mean" , ylab="sd", main = paste0(pairs[[i]], collapse = " , ") , xlim=c(0, 50), ylim=c(0, 50))
abline(a = 0, b = 1, col = "red")
lines(l, col = "blue")
lines(s, col = "green")

dev.off()



d <- unique(data.frame(mean = l$x, sd = l$y))

write.table(d, paste0(out.dir, "SDMeanScale.txt"), quote = FALSE, row.names = FALSE)
save(s, file = paste0(out.dir, "SDMeanScale.RData"))



#################################

v <- NULL
m <- NULL

for(i in 1:length(pairs)){
  # i = 1
  
  v <- c(v, apply(dcpm[, pairs[[i]]], 1, function(x) abs(diff(x))))
  m <- c(m, rowMeans(dcpm[, pairs[[i]]]))
  
}


SDMean.lowess <- l <- lowess(m[m >=1 ],v[m >=1 ])
# s <- smooth.spline(m[m >=1 ], v[m >=1 ], spar=2/3)
s <- smooth.spline(l$x, l$y, spar=2/3)


png(paste0(out.dir, "ABSMeanScale50.png"), 600, 600)

smoothScatter( m, v , nrpoints = Inf, nbin = 500, xlab="mean" , ylab="sd", main = paste0(pairs[[i]], collapse = " , ") , xlim=c(0, 50), ylim=c(0, 50))
abline(a = 0, b = 1, col = "red")
lines(l, col = "blue")
lines(s, col = "green")

dev.off()



d <- unique(data.frame(mean = l$x, sd = l$y))

write.table(d, paste0(out.dir, "ABSMeanScale.txt"), quote = FALSE, row.names = FALSE)
save(s, file = paste0(out.dir, "ABSMeanScale.RData"))


#################################


for(i in 1:length(pairs)){
  # i = 1
  png(paste0(out.dir, "FCMeanScale",paste0(pairs[[i]], collapse = "_"),".png"), 600, 600)
  
  fc <- apply(dcpm[, pairs[[i]]], 1, function(x) max(x)/min(x))
  m <- rowMeans(dcpm[, pairs[[i]]])
  
  plot(log10(m), log10(fc), pch=".")
  
#   smoothScatter( m, fc , nrpoints = Inf, nbin = 1000, xlab="mean" , ylab="sd", main = paste0(pairs[[i]], collapse = " , "))
  abline(a = 0, b = 1, col = "red")
  
  
  dev.off()
  
}



########################### log

dcpm <- log2(cpm(d, normalized.lib.sizes=TRUE))


for(i in 1:length(pairs)){
  # i = 1
  png(paste0(out.dir, "SDMeanLog_",paste0(pairs[[i]], collapse = "_"),".png"), 600, 600)
  
  #   v <- abs(dcpm[, pairs[[i]][1]] - dcpm[,pairs[[i]][2]]) 
  v <- apply(dcpm[, pairs[[i]]], 1, sd)
  m <- rowMeans(dcpm[, pairs[[i]]])
  
  smoothScatter(m, v, nrpoints = Inf, nbin = 500, xlab="mean" , ylab="sd", main = paste0(pairs[[i]], collapse = " , "))
  
  dev.off()
}












