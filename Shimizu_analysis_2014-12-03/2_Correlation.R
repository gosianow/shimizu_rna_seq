#####################################################################################################
# BioC 3.0

### Check the correlation between samples from the same day for 8266 and 970

#####################################################################################################

### load data

RPath <- "/home/gosia/R/R_Shimizu_RNA_seq/Shimizu_analysis_2014-12-03/"
dataPath <- "/home/Shared/data/seq/Shimizu_RNA_seq/Data/"

analysisPath <- "Analysis_2014-12-03"
analysisPath <- paste0("/home/Shared/data/seq/Shimizu_RNA_seq/", analysisPath)
dir.create(analysisPath, showWarnings = FALSE)
setwd(analysisPath)

load("Shimizu_workspace.Rdata")

dir.create("Plots_Correlation/", showWarnings=F, recursive=T)

####################################


all(colnames(x)==rownames(new.samps))

dim(x)


library(edgeR)
d <- DGEList(x, group=new.samps$tree_ID)
d <- calcNormFactors(d)

# make sure a gene is expressed (CPM > 1) in more samples
cps <- cpm(d, normalized.lib.sizes=TRUE)

d <- d[ rowSums( cps > 10 ) > 10, ]
dim(d$counts)

d$counts  <- d$counts + 1
dcpm <- log10(cpm(d, normalized.lib.sizes=TRUE))



new.samps[,c("sample_name","tree_ID", "time_ch")]

pairs <- list(c("H6_8266_20081208", "H7_8266_20081208"), c("H8_8266_20090105", "H9_8266_20090105"), c("E4_8266_20090202", "F4_8266_20090202"), c("N1_8266_20090224", "N5_8266_20090224"), c("I9_970_20081216", "J1_970_20081216"), c("J2_970_20090207", "K9_970_20090207"))


pdf(paste0("Plots_Correlation/" , "Scatter" ,".pdf"))

for(i in 1:length(pairs)){
  # i = 1 
  smoothScatter( dcpm[, pairs[[i]][1] ] , dcpm[, pairs[[i]][2]], nrpoints = Inf, nbin = 500, xlab=pairs[[i]][1], ylab=pairs[[i]][2])
  abline(a = 0, b = 1, col = "red")
  
}

dev.off()



# library(gclus)
# 
# dta <- dcpm
# dta.r <- abs(cor(dta)) # get correlations
# dta.col <- dmat.color(dta.r) # get colors
# # reorder variables so those with highest correlation
# # are closest to the diagonal
# dta.o <- order.single(dta.r) 
# 
# #### takes too much time to plot 
# png(paste0("Plots_Correlation/" , "Pairs" ,".png"), w = 1000, h = 1000)
# cpairs(dta, dta.o, panel.colors=dta.col, gap=.5, main="Variables Ordered and Colored by Correlation" )
# dev.off()



library(corrgram)

pdf(paste0("Plots_Correlation/" , "Corrgram" ,".pdf"), 10, 10)

corrgram(dcpm[,new.samps$tree_ID %in% c("970")], order=TRUE, lower.panel=panel.shade, upper.panel=panel.pie, text.panel=panel.txt, main="", labels = new.samps[new.samps$tree_ID %in% c("970"), "time_ch"])

corrgram(dcpm[,new.samps$tree_ID %in% c("970")], order=FALSE, lower.panel=panel.shade, upper.panel=panel.pie, text.panel=panel.txt, main="", labels = new.samps[new.samps$tree_ID %in% c("970"), "time_ch"])

corrgram(dcpm[,new.samps$tree_ID %in% c("8266")], order=TRUE, lower.panel=panel.shade, upper.panel=panel.pie, text.panel=panel.txt, main="", labels = new.samps[new.samps$tree_ID %in% c("8266"), "time_ch"])

corrgram(dcpm[,new.samps$tree_ID %in% c("8266")], order=FALSE, lower.panel=panel.shade, upper.panel=panel.pie, text.panel=panel.txt, main="", labels = new.samps[new.samps$tree_ID %in% c("8266"), "time_ch"], col.regions=colorRampPalette(c("red","salmon","white","royalblue","navy")))

dev.off()




















