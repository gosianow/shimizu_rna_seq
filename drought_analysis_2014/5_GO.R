
#####################################################################################################

### GO analysis with topGO

#####################################################################################################


setwd("/home/Shared/data/seq/Shimizu_RNA_seq/")
load("Shimizu_workspace.Rdata")

dir.create("GO", showWarnings=F, recursive=T)

library(topGO)
library(org.At.tair.db)


go.name <- "Kmeans_all"
whole.clustering <- read.table("Clustering/Kmeans_all/CPM_norm_all_kmeans_table.xls", sep="\t", header=T)
nr.clusters <- 32


# go.name <- "Kmeans_drought_trees"
# whole.clustering <- read.table("Clustering/Kmeans_drought_trees/CPM_norm_drought_kmeans_table.xls", sep="\t", header=T)
# nr.clusters <- 22



head(whole.clustering)
whole.clustering <- whole.clustering[!is.na(whole.clustering$AT_ID), ]

fun.gene.sel <- function(gene.vector) {
  return(gene.vector <- ifelse(gene.vector==0, FALSE, TRUE))
}

# genes must be in AT_IDs 
assayed.genes <- whole.clustering$AT_ID

allRes.Fisher.merged <- NULL
allRes.Fisher.elim.merged <- NULL



for(cluster in 1:nr.clusters){
  
  # cluster <- 1
  sel.cluster.genes <- whole.clustering$AT_ID[whole.clustering[, paste0("clustering", nr.clusters)]==cluster]
  
  gene.vector=as.numeric(assayed.genes %in% sel.cluster.genes)
  names(gene.vector) <- assayed.genes
  names(assayed.genes) <- gene.vector
  
  
  for(go in c("BP","MF","CC")){
    # go = "BP"
    
    sampleGOdata <- new("topGOdata", description = "Simple session", ontology = go, allGenes = gene.vector, geneSel = fun.gene.sel , nodeSize = 10, annot = annFUN.org, mapping = "org.At.tair.db")
    # Fisher's exact test which is based on gene counts, and a Kolmogorov-Smirnov like test which computes enrichment based on gene scores
    # For the method classic each GO category is tested independently
    resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
    resultFisher.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher")
    
    
    pValues.Fisher <- score(resultFisher)
    #topNodes.Fisher <- sum(pValues.Fisher < 0.05)
    topNodes.Fisher <- length(pValues.Fisher)
    pValues.Fisher.elim <- score(resultFisher.elim)
    #topNodes.Fisher.elim <- sum(pValues.Fisher.elim < 0.0001)
    topNodes.Fisher.elim <- length(pValues.Fisher.elim)
    
    # pdf("GO/pValues.pdf")
    # colMap <- function(x) {
    #   .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
    #   return(.col[match(1:length(x), order(x))])
    # }
    # 
    # gstat <- termStat(sampleGOdata, names(pValues.Fisher))
    # gSize <- gstat$Annotated / max(gstat$Annotated) * 4
    # gCol <- colMap(gstat$Significant)
    # 
    # plot(pValues.Fisher, pValues.Fisher.elim[names(pValues.Fisher)], xlab = "p-value classic", ylab = "p-value elim",pch = 19, cex = gSize, col = gCol)
    # plot(pValues.Fisher, pValues.Fisher.elim[names(pValues.Fisher)], xlab = "p-value classic", ylab = "p-value elim",pch = 19)
    # dev.off()
    
   
      allRes.Fisher <- GenTable(sampleGOdata, classicFisher = resultFisher, elimFisher = resultFisher.elim, orderBy = "classicFisher", ranksOf = "elimFisher", topNodes = topNodes.Fisher)      
      allRes.Fisher$GO <- go
      allRes.Fisher$cluster <- cluster      
      allRes.Fisher.merged <- rbind(allRes.Fisher.merged, allRes.Fisher)
      

      allRes.Fisher.elim <- GenTable(sampleGOdata, classicFisher = resultFisher, elimFisher = resultFisher.elim, orderBy = "elimFisher", ranksOf = "classicFisher", topNodes = topNodes.Fisher.elim)      
      allRes.Fisher.elim$GO <- go
      allRes.Fisher.elim$cluster <- cluster   
      allRes.Fisher.elim.merged <- rbind(allRes.Fisher.elim.merged, allRes.Fisher.elim)

    
    
    #     pdf(paste("GO/GO_", go, ".pdf", sep=""), width = 10, height = 10)
    #     showSigOfNodes(sampleGOdata, score(resultKS), firstSigNodes = 5, useInfo = 'all')
    #     dev.off()
    
    
  }
  
  
}



allRes.Fisher.merged$classicFisher.adj <- p.adjust(p=allRes.Fisher.merged$classicFisher, method = "BH")


write.table(allRes.Fisher.merged, paste("GO/",go.name,"_GO_Fisher" , ".xls", sep=""), sep="\t", row.names=F)

#try(write.table(allRes.Fisher.merged[allRes.Fisher.merged$classicFisher.adj <= 0.05, ], paste("GO/",go.name,"_GO_Fisher_FDR0.05" , ".xls", sep=""), sep="\t", row.names=F))

allRes.Fisher.elim.merged$elimFisher.adj <- p.adjust(p=allRes.Fisher.elim.merged$elimFisher, method = "BH")

write.table(allRes.Fisher.elim.merged, paste("GO/",go.name,"_GO_Fisher_elim" , ".xls", sep=""), sep="\t", row.names=F)

#try(write.table(allRes.Fisher.elim.merged[allRes.Fisher.elim.merged$elimFisher.adj <= 0.05, ], paste("GO/",go.name,"_GO_Fisher_elim_FDR0.05" , ".xls", sep=""), sep="\t", row.names=F))

