#####################################################################################################
# NEW data (from Oct 18th)

# plots of expression
# clustering PAM per tree

#####################################################################################################

#####################################################################################################

### Load data 

#####################################################################################################

# setwd("/home/gosia/Shimizu_RNA_seq/")

### Robust packages
# install.packages(pkgs="R/R_packages_Robust_versions/edgeR_4.0.33.tar.gz", lib="R/Library/")
# install.packages(pkgs="R/R_packages_Robust_versions/limma_4.17.12.tar.gz", lib="R/Library/")

# library(package=limma, lib.loc="R/Library/")
# library(package=edgeR, lib.loc="R/Library/")

setwd("/Users/gosia/Analysis/Shimizu_RNA_seq/")

### samples & factors
# load(file="Samples_out/new_samps_interpolation.RData")
new.samps <- read.table("Samples_out/new_samps_interpolation_November.csv", header=TRUE, sep=";", stringsAsFactors=FALSE)

rownames(new.samps) <- new.samps$sample_name
new.samps$tree_ID <- as.factor(new.samps$tree_ID)
new.samps$tree_col <- new.samps$tree_ID
levels(new.samps$tree_col) <- 1:length(levels(new.samps$tree_col))
levels(new.samps$tree_col) <- apply(col2rgb(levels(new.samps$tree_col)), 2, function(c2r){ rgb(c2r[1], c2r[2], c2r[3], maxColorValue=255)})
new.samps$tree_col <- as.character(new.samps$tree_col)

### counts
x <- read.table("Data/new_data_18Oct/raw_data-clean.csv", sep=",", header=T, row.names=1)
x <- x[, new.samps$sample_name]


flowered.samps <- c("E7_8266_20090416") # flowered sample
control.samps <- rownames(new.samps[new.samps$drough.control=="control",])
November.samps <- c("N9_970_20090114","M1_1099_20090114","I5_8212_20090119","F8_8266_20090119")

###########################
### files with control genes
###########################
# UP.genes <- read.table("Data/genes_descr_control/mDr_Day10_drought_up_regulated_genes_Harb_etal.csv", sep=";", header=T, stringsAsFactors=FALSE)
# DOWN.genes <- read.table("Data/genes_descr_control/mDr_Day10_drought_down_regulated_genes_Harb_etal.csv", sep=";", header=T, stringsAsFactors=FALSE)
# Drought.genes <- read.table("Data/genes_descr_control/mDr_Day10_drought_regulated_genes_Harb_etal.csv", sep=";", header=T, stringsAsFactors=FALSE)

### DE unigenes for S.beccariana overlap with BLAST A.thaliana, all clusters
# MolEcol.DE.genes <- read.table("Data/genes_descr_control/gene_list_in_all_clusters_DEgenes_Mol_Ecol.csv", sep=";", header=T, stringsAsFactors=FALSE)

# AT.id <- read.table("Data/new_data_18Oct/trinity_tair10_obh.csv", sep=",", stringsAsFactors=FALSE)
# ATs.tmp <- apply(AT.id, MARGIN=1, function(x){  
#   # x=AT.id[1,]
#   strsplit(as.character(x[2]), ".", fixed=T)[[1]][1]
# })
# AT.id.clean <- cbind(AT.id[,1], ATs.tmp)
# write.table(AT.id.clean, "Data/new_data_18Oct/trinity_tair10_obh_clean.csv", quote=F, col.names=F, row.names=F,sep=",")

AT.id <- read.table("Data/new_data_18Oct/trinity_tair10_obh_clean.csv", sep=",", stringsAsFactors=FALSE)

# library(org.At.tair.db)
# syms <- mget(AT.id[,2], org.At.tairSYMBOL)
# syms <- sapply(syms, function(u){ if(is.na(u[1])) "" else paste(u,collapse=";") })
# tfd <- read.table("Data/genes_descr_control/TAIR10_functional_descriptions",sep="\t", header=TRUE, quote="", comment.char="", stringsAsFactors=FALSE)
# t10id <- gsub("\\.[1-9]","",tfd$Model_name)
# m <- match(AT.id[,2],t10id)
# descr <- as.character(tfd$Short_description[m])
# genes.description.tmp <- data.frame(Cid=AT.id[,1], ATid=AT.id[,2], ATsymbol=syms, Atdescription=descr)
# Cid <- rownames(x)
# genes.description <- merge(Cid, genes.description.tmp, by=1, all.x=T)
# write.table(genes.description, "Data/new_data_18Oct/genes_description.xls", quote=F, sep="\t", row.names=F)
# # Note: delete ' from excel file

genes.description <- read.table("Data/new_data_18Oct/genes_description.xls", header=T, stringsAsFactors=F, sep="\t")


### Flowering genes - Table S4
Athaliana.flowering.genes <- read.table("Data/genes_descr_control/gene_list_flowering_related_MolEcol_Athaliana_S4.txt", sep="\t", header=F, stringsAsFactors=FALSE)

# Athaliana.flowering.genes[,1] <- gsub(" ", "", Athaliana.flowering.genes[,1])
# # "AT1G65480" - FT gene, "AT2G22540" - SVP gene
# Athaliana.flowering.genes.FT.SVP <- Athaliana.flowering.genes[Athaliana.flowering.genes[,1] %in% c("AT1G65480", "AT2G22540"), ]


# genes.x <- rownames(x)
# genes.full.description <- data.frame(ID=genes.x)
# genes.full.description <- merge(genes.full.description, genes.description, by=1, all.x=TRUE)
# genes.full.description <- merge(genes.full.description, cbind(Flowering="FL", Athaliana.flowering.genes), by.x=2, by.y=2, all.x=TRUE)
# colnames(genes.full.description) <- c("AT_ID", "ID", "AT_symbol", "Description", "Flowering", "Description2")
# genes.full.description <- genes.full.description[,c("ID", "AT_ID","AT_symbol", "Description", "Flowering", "Description2")]
# write.table(genes.full.description, "Data/new_data_18Oct/genes_description_full.txt", sep="\t", quote=F, row.names=F)

genes.full.description <- read.table("Data/new_data_18Oct/genes_description_full.txt", header=T, sep="\t", stringsAsFactors=F)



#####################################################################################################

### plots of expression for flowering genes

#####################################################################################################

#AT.genes <- Athaliana.flowering.genes.FT.SVP
AT.genes <- Athaliana.flowering.genes

#elim.samps=c(flowered.samps)
elim.samps=NULL

new.samps <- new.samps[order(new.samps$time_nr), ]

x <- x[, new.samps$sample_name] 
x <- x[,!names(x) %in% elim.samps]
new.samps <- new.samps[!new.samps$sample_name %in% elim.samps, ]

new.samps$tree_ID <- as.factor(new.samps$tree_ID)

library(edgeR)
d.org <- DGEList(x, group=new.samps$tree_ID)
d.org <- calcNormFactors(d.org)


d.cpm <- cpm(d.org, normalized.lib.sizes=TRUE)
d.cpm.l <- log(d.cpm + min(d.cpm[d.cpm != 0]))


####################################
### all plots in one file
####################################

dir.create("Plots_of_flowering_genes/", showWarnings=F, recursive=T)

pdf(paste0("Plots_of_flowering_genes/" , "Flowering_genes_from_S4_table_Nov_fl" ,".pdf"), h=5, w=10)

info.table <- NULL

for(g in 1:nrow(AT.genes)){
  # g=8
  cat(paste(g, ", "))
  genes <- AT.id[AT.id[,2] == AT.genes[g, 1], 1]
  
  if(length(genes!=0)){
    for(j in 1:length(genes)){
      # j=1
      
      info.table <- rbind(info.table, data.frame(AT.ID=AT.genes[g, 1], AT.description=AT.genes[g, 2] , ID=genes[j]))
      
      plot(0, type="n", main=paste0(AT.genes[g, 1] ," - ", AT.genes[g, 2] , "\n",  genes[j]) ,xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(d.cpm[genes[j], ])), max(na.omit(d.cpm[genes[j], ]))), xlab="Time", ylab="Gene Expression in cpm")
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[246])
      for(t in levels(new.samps$tree_ID)){
        lines(new.samps$time_nr[new.samps$tree_ID==t], d.cpm[genes[j], new.samps$tree_ID==t] , col=which(levels(new.samps$tree_ID)==t), type="b", pch=ifelse(new.samps$drough.control[new.samps$tree_ID==t]=="drought", 18, 16)) 
      }
      legend("topleft", legend = levels(new.samps$tree_ID), col=1:length(levels(new.samps$tree_ID)), cex=1, text.col= 1:length(levels(new.samps$tree_ID)))
      
    }
    
  }
}

plot.vars <- names(new.samps)[20:34]

for(v in plot.vars){
  
  plot(0, type="n", xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(new.samps[,v])), max(na.omit(new.samps[,v]))), xlab="Time", ylab=v)
  for(t in levels(new.samps$tree_ID)){
    lines(new.samps$time_nr[new.samps$tree_ID==t], new.samps[new.samps$tree_ID==t, v], col=which(levels(new.samps$tree_ID)==t), type="b", pch=ifelse(new.samps$drough.control[new.samps$tree_ID==t]=="drought", 18, 16)) 
  }
  legend("topleft", legend = levels(new.samps$tree_ID), col=1:length(levels(new.samps$tree_ID)), cex=0.5, text.col= 1:length(levels(new.samps$tree_ID)))
  
}

# write.table(info.table,"Plots_of_flowering_genes/info_table.xls", quote=FALSE, sep="\t", row.names=FALSE)

dev.off()


### log cpm

pdf(paste0("Plots_of_flowering_genes/" , "Flowering_genes_from_S4_table_Nov_fl_log" ,".pdf"), h=5, w=10)

info.table <- NULL

for(g in 1:nrow(AT.genes)){
  # g=8
  cat(paste(g, ", "))
  genes <- AT.id[AT.id[,2] == AT.genes[g, 1], 1]
  
  if(length(genes!=0)){
    for(j in 1:length(genes)){
      # j=1
      
      info.table <- rbind(info.table, data.frame(AT.ID=AT.genes[g, 1], AT.description=AT.genes[g, 2] , ID=genes[j]))
      
      plot(0, type="n", main=paste0(AT.genes[g, 1] ," - ", AT.genes[g, 2] , "\n",  genes[j]) ,xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(d.cpm.l[genes[j], ])), max(na.omit(d.cpm.l[genes[j], ]))), xlab="Time", ylab="Gene Expression in log(cpm)")
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[246])
      for(t in levels(new.samps$tree_ID)){
        lines(new.samps$time_nr[new.samps$tree_ID==t], d.cpm.l[genes[j], new.samps$tree_ID==t] , col=which(levels(new.samps$tree_ID)==t), type="b", pch=ifelse(new.samps$drough.control[new.samps$tree_ID==t]=="drought", 18, 16)) 
      }
      legend("topleft", legend = levels(new.samps$tree_ID), col=1:length(levels(new.samps$tree_ID)), cex=1, text.col= 1:length(levels(new.samps$tree_ID)))
      
    }
    
  }
}

plot.vars <- names(new.samps)[20:34]

for(v in plot.vars){
  
  plot(0, type="n", xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(new.samps[,v])), max(na.omit(new.samps[,v]))), xlab="Time", ylab=v)
  for(t in levels(new.samps$tree_ID)){
    lines(new.samps$time_nr[new.samps$tree_ID==t], new.samps[new.samps$tree_ID==t, v], col=which(levels(new.samps$tree_ID)==t), type="b", pch=ifelse(new.samps$drough.control[new.samps$tree_ID==t]=="drought", 18, 16)) 
  }
  legend("topleft", legend = levels(new.samps$tree_ID), col=1:length(levels(new.samps$tree_ID)), cex=0.5, text.col= 1:length(levels(new.samps$tree_ID)))
  
}

#write.table(info.table,"Plots_of_flowering_genes/info_table.xls", quote=FALSE, sep="\t", row.names=FALSE)

dev.off()




#####################################################################################################

### clustering

#####################################################################################################

############################################
### select genes for clustering
############################################

# elim.samps=c(flowered.samps, "K4_990_20081205", "K5_990_20090511")
elim.samps=NULL
x <- x[,!names(x) %in% elim.samps]
new.samps <- new.samps[!new.samps$sample_name %in% elim.samps, ]


library(limma)
library(edgeR)

d.org <- DGEList(x, group=new.samps$tree_ID)
d.org <- calcNormFactors(d.org)

### make sure a gene is expressed (CPM > 1) in more than 2 samples
d.cpm.org <- cpm(d.org, normalized.lib.sizes=TRUE)
dim(d.cpm.org)

# sum(d.cpm.org["GID033785_3528015",] > 0.1)
# sum( rowSums(d.cpm.org > 1) > 2 )

# cpm cutoff hihg:
#d.org <- d.org[ rowSums(d.cpm.org>1) > 2, ]
# cpm cutoff relaxed:
d.org <- d.org[ rowSums(d.cpm.org > 1) > 2, ]
d.cpm.org <- cpm(d.org, normalized.lib.sizes=TRUE)
d.cpm.org.l <- log(d.cpm.org  + min(d.cpm.org [d.cpm.org  != 0]))

count.data <- list()
count.data$d.cpm.org <- d.cpm.org
count.data$d.cpm.org.l <- d.cpm.org.l
selected.genes <- list()

# save.image("Clustering_workspace_all_genes.Rdata")

### select DE genes

# estimate dispersion  
design <- model.matrix(~-1+tree_ID, data=new.samps)
d.org <- estimateGLMCommonDisp(d.org,design)
d.org <- estimateGLMTrendedDisp(d.org,design)
d.org <- estimateGLMTagwiseDisp(d.org,design)

plot(d.org$AveLogCPM, d.org$tagwise.dispersion)

# # glmFit
# fit <- glmFit(d.org,design)
# 
# # make contrast to test any difference from control
# # mc <- makeContrasts(contrasts=c("(tree_ID970+tree_ID8212)/2-(tree_ID1099+tree_ID1377)/2", "tree_ID8266-(tree_ID1099+tree_ID1377)/2"),levels=colnames(design))
# mc <- makeContrasts(contrasts=c("(tree_ID970+tree_ID8212+tree_ID8266)/3-(tree_ID1099+tree_ID1377+tree_ID990)/3", "tree_ID8266-(tree_ID1099+tree_ID1377+tree_ID990)/3"),levels=colnames(design))
# colnames(mc) <- c("sheet_vs_control","flower_vs_control")
# 
# # glmLRT conducts likelihood ratio tests for one or more coefficients in the linear model
# lrt <- glmLRT(fit, contrast = mc[,"sheet_vs_control"])
# tt <- topTags(lrt, n=nrow(d.org))
# tr=0.1
# 
# de.genes <- rownames(tt$table[tt$table[,"FDR"] <= tr,])


### genes with high variance
dir.create("Clustering/Variance/", showWarnings=F, recursive=T)

pdf("Clustering/Variance/CPM_Mean-variance_trend.pdf")
mean <- apply(d.cpm.org.l, 1, mean)
sd <- apply(d.cpm.org.l, 1, sd)
plot(mean,sd, pch=16, main="Log(cpm)")

mean <- apply(d.cpm.org, 1, mean)
sd <- apply(d.cpm.org, 1, sd)
plot(mean,sd, pch=16, main="cpm", log="xy")

plotMeanVar(d.org, show.raw.vars=T)

dev.off()


### assings weights based on mean-variance trend - voom

pdf("Clustering/Variance/Voom_estimation_Mean-variance_trend.pdf")
d.voom <- voom(counts=d.org, plot = TRUE, span=0.6)
dev.off()

d.voom.weight <- (d.voom$E * d.voom$weights) # values in log2

# row.weights <- rowSums(d.voom$weights)

# pdf("Clustering/Voom_weighted_Mean-variance_w_trend.pdf")
# d.voom.mean <- unlist(lapply(1:nrow(d.voom.weight), function(v){ 
#   sum(d.voom.weight[v, ])/row.weights[v] 
# }))
# d.voom.sd <- unlist(lapply(1:nrow(d.voom.weight), function(v){ 
#   sqrt(sum(d.voom$weights[v,]*(d.voom$E[v, ]-d.voom.mean[v])^2)/row.weights[v])
# }))
# plot(d.voom.mean, d.voom.sd)
# dev.off()

pdf("Clustering/Variance/Voom_Mean-variance_trend.pdf")
d.voom.mean <- apply(d.voom.weight, 1, mean)
d.voom.sd <- apply(d.voom.weight, 1, sd)
#l <- lowess(d.voom.mean, d.voom.sd)
plot(d.voom.mean, d.voom.sd, pch=16, main="Voom weighted values")
#lines(l, col=2)
dev.off()


names(d.voom.sd) <- rownames(d.voom$E)
selected.genes$var.genes.voom <- names(sort(d.voom.sd, decreasing=TRUE)[1:5000])
count.data$d.voom.weight <- d.voom.weight


pdf("Clustering/Variance/Voom_sel_CPM_Mean-variance_trend.pdf")
mean <- apply(d.cpm.org.l, 1, mean)
sd <- apply(d.cpm.org.l, 1, sd)
plot(mean,sd, pch=16, main="Log(cpm)")
points(mean[selected.genes$var.genes.voom], sd[selected.genes$var.genes.voom], pch=16, col=2)
points(mean[genes.full.description[!is.na(genes.full.description$Flowering), "ID"]], sd[genes.full.description[!is.na(genes.full.description$Flowering), "ID"]], pch=16, col=3)

mean <- apply(d.cpm.org, 1, mean)
sd <- apply(d.cpm.org, 1, sd)
plot(mean,sd, pch=16, main="cpm", log="xy")
points(mean[selected.genes$var.genes.voom], sd[selected.genes$var.genes.voom], pch=16, col=2)
points(mean[genes.full.description[!is.na(genes.full.description$Flowering), "ID"]], sd[genes.full.description[!is.na(genes.full.description$Flowering), "ID"]], pch=16, col=3)

dev.off()


### variance stabilization - DESeq 

library(DESeq)
cds <- newCountDataSet(countData=d.org$counts, conditions=new.samps$sample_ID)
cds <- estimateSizeFactors( cds )
cds <- estimateDispersions( cds, method="blind" )
vsd <- getVarianceStabilizedData( cds )

pdf("Clustering/Variance/DESeq_Mean-variance_trend.pdf")
plot(rowMeans(vsd), sqrt(genefilter::rowVars(vsd)), main="DESeq VSD")
dev.off()


selected.genes$var.genes.vs <- names(sort(genefilter::rowVars(vsd), decreasing=TRUE)[1:5000])
count.data$vsd <- vsd

# vsd["GID063259_3539580",]

pdf("Clustering/Variance/DESeq_sel_CPM_Mean-variance_trend.pdf")

mean <- apply(d.cpm.org.l, 1, mean)
sd <- apply(d.cpm.org.l, 1, sd)
plot(mean,sd, pch=16, main="Log(cpm)")
points(mean[selected.genes$var.genes.vs], sd[selected.genes$var.genes.vs], pch=16, col=2)
points(mean[genes.full.description[!is.na(genes.full.description$Flowering), "ID"]], sd[genes.full.description[!is.na(genes.full.description$Flowering), "ID"]], pch=16, col=3)

mean <- apply(d.cpm.org, 1, mean)
sd <- apply(d.cpm.org, 1, sd)
plot(mean,sd, pch=16, main="cpm", log="xy")
points(mean[selected.genes$var.genes.vs], sd[selected.genes$var.genes.vs], pch=16, col=2)
points(mean[genes.full.description[!is.na(genes.full.description$Flowering), "ID"]], sd[genes.full.description[!is.na(genes.full.description$Flowering), "ID"]], pch=16, col=3)

dev.off()


"GID053198_3536068" %in% selected.genes$var.genes.vs



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

# save.image("Clustering_workspace.Rdata")

############################################
# load clustering workspace on server
############################################

setwd("/home/gosia/Shimizu_RNA_seq/")
dir.create("Clustering", showWarnings=F, recursive=T)
load("Clustering_workspace.Rdata")

# install.packages("/home/gosia/R/packages/", "/home/gosia/R/library/3.0.1/")
library(cluster)
library(clusterSim, lib.loc="/home/gosia/R/library/3.0.1/")
library(lattice)
library(parallel)
library(gtools)
library(gplots) 
library(RColorBrewer)
source("R/heatmap2.r")

#####################

setwd("/Users/gosia/Analysis/Shimizu_RNA_seq/")
load("Clustering_workspace.Rdata")

library(cluster)
library(clusterSim)
library(lattice)
library(parallel)
library(gplots) 
library(RColorBrewer)
source("R/heatmap2.r")
#library(ggplot2)



############################################
# PAM functions
############################################


PAM.clustering <- function(pam.x, new.samps, genes.full.description,  dist.method= "spearman", norm.method="norm", out.path="Clustering/PAM/", out.name="", prior.nr.cl=2:5, mc.cores=1, ylim=c(-3,4), clustering.samps=NA){
  
  dir.create(out.path, showWarnings=FALSE, recursive=TRUE)
  
  samps.order <- row.names(new.samps[order(new.samps$drough.control, new.samps$tree_ID, new.samps$time_nr), ])
  pam.x <- pam.x[, samps.order]
  new.samps <- new.samps[samps.order,]
  if(!is.na(clustering.samps))
  new.samps <- new.samps[clustering.samps,]
  
  if(norm.method!="none")
    pam.x <- normalize.counts(counts=pam.x, norm.method=norm.method)
  
  pam.x <- na.omit(pam.x)

  pam.df <- as.data.frame(pam.x[, clustering.samps])
  
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
    # i=8
    cat("Clustering into ", i, "clusters \n")
    pam.obj <- pam(x=pam.dist, k=i) 
    pam.cl <- pam.obj$clustering
    
    pam.obj$index.S <- index.S(d=pam.dist, cl=pam.cl)
    pam.obj$index.G1 <- index.G1(x=pam.x, d=pam.dist, cl=pam.cl, centrotypes="medoids")
    
    pam.df[, paste0("clustering", i)] <- pam.cl
    
    
    #     pdf(paste0(out.path, "/",out.name,"_Sillhouette_plot_k" ,i, ".pdf"), w=15, h=10)
    #     pam.sill <- silhouette(x=pam.cl, dist=pam.dist)
    #     plot(pam.sill, col=colorRampPalette(brewer.pal(12,"Set3"))(i))
    #     dev.off()
    
    #     png(paste0(out.path,"/",out.name,"_Parallel1_plot_k" ,i, ".png"), w=600, h=1000)
    #     plot(parallel(~pam.df[, -ncol(pam.df)], pam.df, groups=clustering))
    #     dev.off()
    #     
    #     png(paste0(out.path,"/",out.name,"_Parallel2_plot_k" ,i, ".png"), w=1000, h=1000)
    #     plot(parallel(~pam.df[, -ncol(pam.df)] | clustering, pam.df))
    #     dev.off()
    
    #     pdf(paste0(out.path, "/",out.name,"_Parcoord_plot_k" ,i, ".pdf"), w=10, h=5)
    #     for(j in 1:i){
    #       parcoord(x=pam.x[pam.cl==j, ], col = colorRampPalette(brewer.pal(12,"Set3"))(i)[j], lty = 1, var.label = FALSE)
    #     }
    #     dev.off()
    
    #     pdf(paste0(out.path, "/",out.name,"_Parcoord_plot_k" ,i, ".pdf"), w=10, h=5)
    #     for(j in 1:i){
    #      plot(1:ncol(pam.x), rep(1, ncol(pam.x)), type="n", ylim=c(ylimL,ylimU))
    #       apply(pam.x[pam.cl==j, ], 1, function(g){
    #         lines(1:ncol(pam.x), g, col=colorRampPalette(brewer.pal(12,"Set3"))(i)[j])
    #       }) 
    #     }
    #     dev.off()
    
    pam.df.sort <- pam.df[order(pam.df[, paste0("clustering", i)]), ]
    
    png(paste0(out.path,"/",out.name,"_Heatmap_plot_k" ,i, ".png"), w=1000, h=1000)
    heatmap.2axis(x=as.matrix(pam.df.sort[,-ncol(pam.df.sort)]), Rowv=NA, Colv=NA, dendrogram="none", labRow=NA, labCol=NULL, scale="none", trace="none", col=colorRampPalette(brewer.pal(9,"Blues"))(100), key=T, keysize =c(0.5, 0.5), density.info="none", margins = c(5,5), xlab = "Samples", ylab = "Clusters of DE Genes", RowSideColors=colorRampPalette(brewer.pal(12,"Set3"))(i)[pam.df.sort[, ncol(pam.df.sort)]], ColSideColors=new.samps$tree_col) 
    dev.off()
    
    
    
    # ColSideColors=colorRampPalette(brewer.pal(8,"Pastel2"))(length(levels(new.samps$tree_ID)))[new.samps[samps.order, "tree_ID"]]
    
    #pam.medoids <- pam.obj$medoids
    
    #     png(paste0(out.path,"/",out.name,"_Heatmap_plot_k" ,i, "b.png"), w=1000, h=1000)
    #     heatmap.2axis(x=as.matrix(pam.df.sort[pam.medoids,-ncol(pam.df.sort)]), Rowv=NA, Colv=NA, dendrogram="none", labRow=NA, labCol=NULL, scale="none", trace="none", col=colorRampPalette(brewer.pal(9,"Blues"))(100), key=T, keysize =c(0.5, 0.5), density.info="none", margins = c(5,5), xlab = "Samples", ylab = "Medoids", RowSideColors=colorRampPalette(brewer.pal(12,"Set3"))(i)[pam.df.sort[pam.medoids, ncol(pam.df.sort)]], ColSideColors=colorRampPalette(brewer.pal(8,"Pastel2"))(length(levels(new.samps$tree_ID)))[new.samps[samps.order, "tree_ID"]]) 
    #     dev.off()
    
    
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
    
    
    pdf(paste0(out.path,"/",out.name,"_Expression_Flowering_genes_k" ,i, ".pdf"), w=10, h=5)
    for(c in 1:i){
      # c=5
      pam.df.tmp <- pam.df[pam.cl==c, -ncol(pam.df)]
      pam.df.tmp <-  pam.df.tmp[, rownames(new.samps)]
      
      pam.df.tmp.fl <- merge(cbind(rownames(pam.df.tmp), pam.df.tmp), genes.flowering, by=1, all=FALSE)
      
      if(nrow(pam.df.tmp.fl) >= 1){
        for(j in 1:nrow(pam.df.tmp.fl)){
          # j=1
          expr <- pam.df.tmp.fl[j, rownames(new.samps)]
          
          plot(0, type="n", main=paste("Clustering into ", i, " groups: Cluster ", c, "\n", pam.df.tmp.fl[j,1], pam.df.tmp.fl[j,"AT_ID"], "\n", pam.df.tmp.fl[j,"Description2"]) ,xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min.pam.df, max.pam.df), xlab="Time", ylab="Normalized Gene Expression", col.main=colorRampPalette(brewer.pal(12,"Set3"))(i)[c])
          rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[246])
          for(t in levels(new.samps$tree_ID)){
            lines(new.samps$time_nr[new.samps$tree_ID==t], expr[new.samps$tree_ID==t] , col=which(levels(new.samps$tree_ID)==t), type="b", pch=ifelse(new.samps$drough.control[new.samps$tree_ID==t]=="drought", 18, 16))    
          }
          legend("topleft", legend = levels(new.samps$tree_ID), col=1:length(levels(new.samps$tree_ID)), cex=1, text.col= 1:length(levels(new.samps$tree_ID)), pch=c(18, 16, 16, 16, 18, 18))
          
        }
      }
      
      
      
    } 
    dev.off()
    
    
    return(pam.obj)
    
  }, mc.cores=mc.cores)
  
  
  save(pam.obj.all, file=paste0(out.path,"/",out.name,"_pam_obj.Rdata"))
  
  cat("Calculate indexes... \n")
  
  Gap.inds <- mclapply(1:(length(pam.obj.all)-1), function(p){ 
    # p=2
    clall=cbind(pam.obj.all[[p]]$clustering, pam.obj.all[[p+1]]$clustering)
    pam.indexGap.pc  <- index.Gap(x=pam.x, d=pam.dist, clall=clall , method="pam", centrotypes="medoids", reference.distribution="pc")$diffu
    pam.indexGap.unif <- index.Gap(x=pam.x, d=pam.dist, clall=clall, method="pam", centrotypes="medoids", reference.distribution="unif")$diffu
    
    return(c(pam.indexGap.pc, pam.indexGap.unif))
    
  }, mc.cores=mc.cores)
  
  
  pam.indexS <- unlist(lapply(pam.obj.all, function(p) p$index.S))
  pam.indexG1 <- unlist(lapply(pam.obj.all, function(p) p$index.G1))
  pam.indexGap.pc <- unlist(lapply(Gap.inds, function(p) p[1]))
  pam.indexGap.unif <- unlist(lapply(Gap.inds, function(p) p[2]))
  
  save(pam.obj.all, pam.indexS, pam.indexG1, pam.indexGap.pc, pam.indexGap.unif, file=paste0(out.path,"/",out.name,"_pam_obj.Rdata"))
  
  pdf(paste0(out.path,"/",out.name,"_pam_idexes.pdf"))
  plot(prior.nr.cl, pam.indexS, main="Silhouette index - max", xlab="nr of clusters")
  plot(prior.nr.cl, pam.indexG1, main="Calinski-Harabasz index - max", xlab="nr of clusters")
  plot(prior.nr.cl[-length(prior.nr.cl)], pam.indexGap.pc, main="Gap pc statistic - max", xlab="nr of clusters")
  plot(prior.nr.cl[-length(prior.nr.cl)], pam.indexGap.unif, main="Gap unif statistic - max", xlab="nr of clusters")
  dev.off()
  
  pam.cls <- do.call("cbind", lapply(pam.obj.all, function(p) p$clustering))
  colnames(pam.cls) <- paste0("clustering", prior.nr.cl)
  pam.df <- cbind(ID=rownames(pam.df), pam.cls,  pam.df)
  pam.df <- merge(genes.full.description, pam.df, by=1, all.y=T)
  
  write.table(pam.df, paste0(out.path,"/",out.name,"_pam_table.csv" ), sep=";", row.names=F, quote=F)
  
  
}


############################################
# PAM run
############################################

# sel genes: DESeq vsd, data: cpm

pam.x <- count.data$d.cpm.org[selected.genes$var.genes.vs, ]


PAM.clustering(pam.x=pam.x, new.samps=new.samps, genes.full.description, dist.method= "euclidian", norm.method="norm", out.path="Clustering/PAM_DESeq_CPM/euclidian_norm/", out.name="8266", prior.nr.cl=2:20, mc.cores=1, ylim=c(-3,6), clustering.samps=new.samps$sample_name[new.samps$tree_ID=="8266"])


PAM.clustering(pam.x=pam.x, new.samps=new.samps, genes.full.description, dist.method= "euclidian", norm.method="norm", out.path="Clustering/PAM_DESeq_CPM/euclidian_norm/", out.name="970", prior.nr.cl=2:20, mc.cores=1, ylim=c(-3,6), clustering.samps=new.samps$sample_name[new.samps$tree_ID=="970"])


# flowering clustering

pam.x <- count.data$d.cpm.org[genes.full.description[!is.na(genes.full.description$Flowering), "ID"], ]

PAM.clustering(pam.x, new.samps, genes.full.description, dist.method= "euclidian", norm.method="norm", out.path="Clustering/PAM_flowering/", out.name="FL_CPM_EU_norm", prior.nr.cl=2:20, mc.cores=1, ylim=c(-3,6))

