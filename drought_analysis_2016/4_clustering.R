##############################################################################
## <<4_clustering.R>>

# BioC 3.3
# Created 4 Aug 2016

##############################################################################
Sys.time()
##############################################################################

library(limma)
library(edgeR)
library(DESeq2)
library(ggplot2)
library(reshape2)
library(BiocParallel)
library(pheatmap)
library(FlowSOM) # for SOM
library(kohonen) # for som
library(clusterSim) # for index.Gap
library(ConsensusClusterPlus)
library(RColorBrewer)

##############################################################################
# Test arguments
##############################################################################

rwd='/home/Shared/data/seq/shimizu_rna_seq'
rcode='/home/gosia/R/shimizu_rna_seq/drought_analysis_2016'
clust_data_dir='Analysis_2016-06-22/4_Clustering/data/expr_rlog_deseq2_interaction_any_time_any_tree.txt'
out_dir='Analysis_2016-06-22/4_Clustering/expr_rlog_deseq2_interaction_any_time_any_tree'
normalize_counts_function <- "normalize_counts.R"

##############################################################################
# Read in the arguments
##############################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(args)

##############################################################################

setwd(rwd)

dir.create(file.path(out_dir, "clustering"), recursive = TRUE)

source(paste0(rcode, "/", normalize_counts_function))

rand_seed <- 7
maxK <- 20

# -----------------------------------------------------------------------------
### Load data
# -----------------------------------------------------------------------------


data_org <- read.table(clust_data_dir, header = TRUE, sep = "\t", as.is = TRUE)

rownames(data_org) <- data_org$gene_id

data <- as.matrix(data_org[, -1])


data_norm <- normalize.counts(data, norm.method = "norm")

thr <- 2
data_norm[data_norm < -thr] <- -thr
data_norm[data_norm > thr] <- thr


# --------------------------------------------------------------------------
# run Levine et al. 2015 pca scoring,
# --------------------------------------------------------------------------

z <- data_norm

npc <- 5

pr <- prcomp(z)  # default is to center but not scale

score <- outer( rep(1,ncol(z)), pr$sdev[1:npc]^2 ) * abs(pr$rotation[, 1:npc])

score_sum <- rowSums(score)

scores_out <- data.frame(sample_id = names(score_sum), pca_score = score_sum)

write.table(scores_out, file = file.path(out_dir, paste0("pca_scores.xls")), sep="\t", row.names=FALSE, quote=FALSE)


## plot PCA scores

scores_out$sample_id <- factor(scores_out$sample_id, levels = scores_out$sample_id)

ggp <- ggplot(scores_out, aes(x = sample_id, y = pca_score)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=12))


pdf(file.path(out_dir, paste0("pca_scores.pdf")), width = 10, height = 7)
print(ggp)
dev.off()


samps2keep <- score_sum > 2
table(samps2keep)


# -----------------------------------------------------------------------------
# SOM from FlowSOM pckg
# -----------------------------------------------------------------------------


set.seed(rand_seed)
flowsom_som <- FlowSOM::SOM(data_norm)


flowsom_codes <- flowsom_som$codes
rownames(flowsom_codes) <- 1:nrow(flowsom_codes)


# plot codes
cluster_rows <- hclust(dist(flowsom_codes), method = "average")

pheatmap(flowsom_codes, cluster_cols = FALSE, cluster_rows = cluster_rows, border_color = NA, fontsize_row = 6, fontsize_col = 10, fontsize = 7, filename = file.path(out_dir, paste0("flowsom_codes_wd.pdf")), width = 10, height = 10)



# ----------------------------------------
# SOM + ConsensusClusterPlus
# ----------------------------------------


### Run ConsensusClusterPlus on SOM codes

# ## try out writing my own function
# my_hc<- function(this_dist,k){
#   tmp <- hclust(this_dist, method = "complete")
#   assignment <- cutree(tmp, k)
#   return(assignment)
# }
# 
# pdf(file.path(out_dir, paste0("flowsom_ConsensusClusterPlus_hc.pdf")), width = 7, height = 7)
# 
# codes_ccp <- ConsensusClusterPlus(t(flowsom_codes), maxK = maxK, reps = 100, pItem = 0.9, pFeature = 1, plot = NULL, verbose = FALSE, clusterAlg = "my_hc", innerLinkage="average", finalLinkage="average", distance = "euclidean", seed = rand_seed)
# 
# dev.off()


## using hc + average linkage
pdf(file.path(out_dir, paste0("flowsom_ConsensusClusterPlus_hc.pdf")), width = 7, height = 7)

codes_ccp <- ConsensusClusterPlus(t(flowsom_codes), maxK = maxK, reps = 100, pItem = 0.9, pFeature = 1, plot = NULL, verbose = FALSE, clusterAlg = "hc", innerLinkage="average", finalLinkage="average", distance = "euclidean", seed = rand_seed)

dev.off()


pdf(file.path(out_dir, paste0("flowsom_ConsensusClusterPlus_hc_ICL.pdf")), width = 7, height = 7)

codes_icl <- calcICL(codes_ccp)

dev.off()



pdf(file.path(out_dir, paste0("flowsom_ConsensusClusterPlus_hc_pheatmap.pdf")), width = 10, height = 7)

for(k in 2:maxK){
  # k <- 5
  
  fm <- codes_ccp[[k]]$ml
  cluster_rows <- hclust(as.dist( 1 - fm ), method = "average")
  
  pheatmap(flowsom_codes, cluster_cols = FALSE, cluster_rows = cluster_rows, border_color = NA, fontsize_row = 6, fontsize_col = 10, fontsize = 7, filename = NA, main = paste0("Consensus cluster number k = ", k))
  
}

dev.off()




### Save the clustering results for genes

clust_out <- data.frame(gene_id = rownames(data_norm), stringsAsFactors = FALSE)

for(k in 2:maxK){
  # k = 2
  
  ## get cluster ids for the codes
  clust_codes <- codes_ccp[[k]]$consensusClass
  
  ## get cluster ids for the genes
  clust <- clust_codes[flowsom_som$mapping[,1]]
  
  clust_out[, paste0("clust_", k)] <- clust
  
}

write.table(clust_out, file.path(out_dir, "clustering", paste0("flowsom_ConsensusClusterPlus_hc.txt")), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



# ----------------------------------------
# SOM + kmeans + GAP index
# ----------------------------------------

clust_out <- data.frame(gene_id = rownames(data_norm), stringsAsFactors = FALSE)
clust_codes_out <- data.frame(code_id = 1:nrow(flowsom_codes), stringsAsFactors = FALSE)
wss_codes <- rep(0, maxK-1)


for(k in 2:maxK){
  # k = 3
  
  ## run kmeans
  codes_kmeans <- kmeans(flowsom_codes, centers = k, iter.max = 100)
  ## get cluster ids for the codes
  clust_codes <- codes_kmeans$cluster
  ## get cluster ids for the genes
  clust <- clust_codes[flowsom_som$mapping[,1]]
  
  clust_out[, paste0("clust_", k)] <- clust
  clust_codes_out[, paste0("clust_", k)] <- clust_codes
  
  wss_codes[k-1] <- sum(codes_kmeans$withinss)
  
  
}

write.table(clust_out, file.path(out_dir, "clustering", paste0("flowsom_kmeans.txt")), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


## plot the within sum of squares 
pdf(file.path(out_dir, paste0("flowsom_kmeans_wss_codes.pdf")), width = 10, height = 7)
plot(2:maxK, wss_codes)
dev.off()



## calculate the GAP index

gap_codes <- rep(0, maxK-2)

for(k in 2:(maxK-1)){
  # k = 2
  gap_codes[k-1] <- index.Gap(flowsom_codes, clall = clust_codes_out[, paste0("clust_", c(k, k+1))], method="k-means", centrotypes="centroids", reference.distribution="unif", B=100)$diffu
  
}


## plot the gap index 

pdf(file.path(out_dir, paste0("flowsom_kmeans_gap_codes.pdf")), width = 10, height = 7)

nc <- 2:(maxK-1)
plot(nc, gap_codes, main=paste0("Gap unif statistic \n final nr of clusters = min{nc: diffu(nc, nc+1) >= 0} \n "), xlab="nc", ylab="diffu", pch=19, col = ifelse(gap_codes >= 0, 2, 1), xaxt = "n")
axis(side=1, at=nc, labels=nc)
abline(h=0)
abline(v=nc, lty=2, lwd=0.5, col = ifelse(gap_codes >= 0, 2,"grey"))

dev.off()



## plot pheatmaps with clustering

pdf(file.path(out_dir, paste0("flowsom_kmeans_pheatmap.pdf")), width = 10, height = 10)

for(k in 2:maxK){
  # k = 2
  
  clust_codes <- clust_codes_out[, paste0("clust_", k)]
  clust_order <- order(clust_codes, decreasing = FALSE)
  
  annotation_row <- data.frame(clustering = factor(clust_codes), row.names = 1:nrow(clust_codes_out))
  
  pheatmap(flowsom_codes[clust_order, ], cluster_cols = FALSE, cluster_rows = FALSE, annotation_row = annotation_row, gaps_row = cumsum(table(annotation_row$clustering)), annotation_legend = FALSE, border_color = NA, fontsize_row = 6, fontsize_col = 10, fontsize = 6, filename = NA, main = paste0("Clustering into ", k, " groups"))
  
}

dev.off()



# ------------------------------------------------------------
# SOM + PAM + silhouette index
# ------------------------------------------------------------


clust_out <- data.frame(gene_id = rownames(data_norm), stringsAsFactors = FALSE)
silhouette_codes <- rep(0, maxK-1)


pdf(file.path(out_dir, paste0("flowsom_pam_silhouette.pdf")), width = 10, height = 7)

for(k in 2:maxK){
  # k = 3
  
  ## run PAM
  pam_dist <- dist(x = flowsom_codes, method = "euclidean")
  codes_pam <- pam(pam_dist, k=k)
  ## get cluster ids for the codes
  clust_codes <- codes_pam$clustering
  ## get cluster ids for the genes
  clust <- clust_codes[flowsom_som$mapping[,1]]
  
  clust_out[, paste0("clust_", k)] <- clust
  
  ## get the silhouette index
  pam_sill <- silhouette(clust_codes, dist = pam_dist)
  silhouette_codes[k-1] <- summary(pam_sill)$avg.width
  
  plot(pam_sill)
  
}

dev.off()

write.table(clust_out, file.path(out_dir, "clustering", paste0("flowsom_pam.txt")), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


## plot the within sum of squares 
pdf(file.path(out_dir, paste0("flowsom_pam_silhouette_codes.pdf")), width = 10, height = 7)
plot(2:maxK, silhouette_codes)
dev.off()





# -----------------------------------------------------------------------------
# som from kohonen pckg - but it does not create codes that cluster nicely...
# -----------------------------------------------------------------------------

## create the SOM grid 
koho_grid <- somgrid(xdim = 10, ydim = 10, topo = "rectangular")

## train the SOM

# koho_som <- som(data_norm, grid = koho_grid, rlen = 200, radius = flowsom_som$radius[[1]])

koho_som <- som(data_norm, grid = koho_grid, rlen = 200)

koho_codes <- koho_som$codes


## som default plots

color_palette <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))

pdf(file.path(out_dir, paste0("kohonen_som.pdf")))

plot(koho_som, type="changes", main = "changes")

plot(koho_som, type="counts", main = "counts", palette.name = color_palette)

plot(koho_som, type="dist.neighbours", main = "dist.neighbours", palette.name = color_palette)

plot(koho_som, type="codes", main = "codes", palette.name = color_palette)

dev.off()



# plot codes

cluster_rows <- hclust(dist(koho_codes), method = "average")

pheatmap(koho_codes, cluster_cols = FALSE, cluster_rows = cluster_rows, border_color = NA, fontsize_row = 6, fontsize_col = 10, fontsize = 7, filename = file.path(out_dir, paste0("kohonen_codes.pdf")), width = 10, height = 10)





























