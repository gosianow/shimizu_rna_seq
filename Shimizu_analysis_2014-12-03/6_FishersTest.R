
#####################################################################################################

### Fisher's exact tests
### BioC 3.0
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

load(paste0("Plots_MeanForReplicates/" , "data.reduced" ,".RData"))


##############################
### Athaliana_flowering_genes
##############################


out.dir <- "Plots_FishersTest/Athaliana_flowering_genes/"
dir.create(out.dir, showWarnings=F, recursive=T)


gene.set <- genes.full.description[!is.na(genes.full.description$Flowering), "AT_ID"] 
name.out <- "Athaliana_flowering_genes"


trees <- as.character(unique(new.samps.red$tree_ID))
nr.clusters <- c(17, 9, 16, 15)
names(nr.clusters) <- trees


for(i in 1:length(trees)){
  # i = 4
  
  go.name <- paste0("Kmeans_tree", trees[i])
  
  whole.clustering.org <- read.table(paste0("Plots_PerTree_Clustering/", "/Tree", trees[i],"/", "Kmeans_all_clustering.txt"), sep="\t", header=T)
  
  whole.clustering <- data.frame(ID = rownames(whole.clustering.org), AT_ID = genes.full.description[rownames(whole.clustering.org), "AT_ID"] , clustering = whole.clustering.org[, paste0("CL", nr.clusters[i])])
  
  whole.clustering <- whole.clustering[!is.na(whole.clustering$AT_ID),]

  
  fisher.pvalues <- matrix(0, 1 ,nr.clusters[i])
  colnames(fisher.pvalues) <- paste0( "CL", 1:nr.clusters[i])
  
  
  fisher.table <- matrix(0,nr.clusters[i],6)
  colnames(fisher.table) <- c("X1", "X2", "X3", "X4", "pvalues", "cluster.size")
  
  for(cluster in 1:nr.clusters[i]){
    
    #cluster <- 1
    
    clustering.tmp <- whole.clustering
    
    clustering.tmp$gene.set <- ifelse(whole.clustering$AT_ID %in%  gene.set, 1, 2)
    clustering.tmp$clustering <- ifelse(whole.clustering$clustering==cluster, 1, 2)
    
    
    fisher.x <- as.matrix(with(clustering.tmp, table(gene.set, clustering)))
    
    fisher.out <- fisher.test(fisher.x, alternative="greater")
    
    fisher.table[cluster, "pvalues"] <- fisher.out$p.value
    fisher.pvalues[1, cluster] <- fisher.out$p.value
    
    if(any(fisher.x < 15)){    
      fisher.table[cluster, "pvalues"] <- NA
      fisher.pvalues[1, cluster] <- NA
    }
    
    
    fisher.table[cluster, "X1"] <- fisher.x[1,1]
    fisher.table[cluster, "X2"] <- fisher.x[1,2]
    fisher.table[cluster, "X3"] <- fisher.x[2,1]
    fisher.table[cluster, "X4"] <- fisher.x[2,2]
    fisher.table[cluster,  "cluster.size"] <- fisher.x[1,1] + fisher.x[2,1]
    
    
  }
   
  cl.ord <- read.table(paste0("Plots_PerTree_Clustering/", "/Tree", trees[i],"/Clustering/", "CL",nr.clusters[i],"Kmeans_ClustOrder.txt"), sep="\t", header=FALSE)
  
  cl.ord <- as.numeric(gsub("CL","" ,cl.ord[,2]))
  
  fisher.table <- data.frame(cluster=1:nr.clusters[i], fisher.table)[cl.ord,]
  
  write.table(fisher.table, paste0(out.dir,"/",  go.name , "_" , name.out, "_Fisher_pvalues.xls"), sep="\t", row.names=F)
  
  
  fisher.pvalues.org <- fisher.pvalues
  
  fisher.pvalues <- fisher.pvalues.org
  
  
  fisher.pvalues.adj <- p.adjust(p=fisher.pvalues, method = "BH")
  fisher.pvalues.adj <- matrix(fisher.pvalues.adj, dim(fisher.pvalues)[1], dim(fisher.pvalues)[2])
  colnames(fisher.pvalues.adj) <- paste0( "CL", 1:nr.clusters[i])
  
  fisher.pvalues.adj <- data.frame(gene.set="Athaliana.flowering.genes", fisher.pvalues.adj)
  
  fisher.pvalues <- data.frame(gene.set="Athaliana.flowering.genes", fisher.pvalues)
  
  write.table(fisher.pvalues, paste0(out.dir,"/", go.name , "_" ,"ALL_PVALUES_Fisher_pvalues.xls"), sep="\t", row.names=F)
  
  write.table(fisher.pvalues.adj, paste0(out.dir,"/", go.name , "_" ,"ALL_PVALUES_Fisher_pvalues_ADJ.xls"), sep="\t", row.names=F)
  
  
}
  


##############################
### gene sets
##############################


out.dir <- "Plots_FishersTest/Gene_sets/"
dir.create(out.dir, showWarnings=F, recursive=T)

gene.sets.dir <- paste0(dataPath, "GeneControlSets/gene_sets_for_Gosia/")
gene.sets.all <- dir(gene.sets.dir)

trees <- as.character(unique(new.samps.red$tree_ID))
nr.clusters <- c(17, 9, 16, 15)
names(nr.clusters) <- trees


for(i in 1:length(trees)){
  # i = 1
  
  go.name <- paste0("Kmeans_tree", trees[i])
  
  whole.clustering.org <- read.table(paste0("Plots_PerTree_Clustering/", "/Tree", trees[i],"/", "Kmeans_all_clustering.txt"), sep="\t", header=T)
  
  whole.clustering <- data.frame(ID = rownames(whole.clustering.org), AT_ID = genes.full.description[rownames(whole.clustering.org), "AT_ID"] , clustering = whole.clustering.org[, paste0("CL", nr.clusters[i])])
  
  whole.clustering <- whole.clustering[!is.na(whole.clustering$AT_ID),]

  fisher.pvalues <- matrix(0, length(gene.sets.all),nr.clusters[i])
  colnames(fisher.pvalues) <- paste0( "CL", 1:nr.clusters[i])
  
  for(j in 1:length(gene.sets.all)){
    # j=29
    
    cat("\n", gene.sets.all[j], "\n") 
    file.data <- read.table(paste0(gene.sets.dir, "/", gene.sets.all[j]), header=T, sep=";", stringsAsFactors = FALSE)
    
    print(head(file.data))
    
    gene.set <- file.data[,1]
    name.out <- gsub(".csv", "",gene.sets.all[j])
    
    fisher.table <- matrix(0,nr.clusters[i],6)
    colnames(fisher.table) <- c("X1", "X2", "X3", "X4", "pvalues", "cluster.size")
    
    for(cluster in 1:nr.clusters[i]){
      
      # cluster <- 1
      
      clustering.tmp <- whole.clustering
      
      clustering.tmp$gene.set <- ifelse(whole.clustering$AT_ID %in%  gene.set, 1, 2)
      clustering.tmp$clustering <- ifelse(whole.clustering$clustering==cluster, 1, 2)
      
      
      fisher.x <- as.matrix(with(clustering.tmp, table(gene.set, clustering)))
      
      ### !!! very inportant to set up alternative hipothesis
      fisher.out <- fisher.test(fisher.x, alternative="greater")
      
      fisher.table[cluster, "pvalues"] <- fisher.out$p.value
      fisher.pvalues[j, cluster] <- fisher.out$p.value
      
      if(any(fisher.x <= 15)){
        fisher.table[cluster, "pvalues"] <- NA
        fisher.pvalues[j, cluster] <- NA
      }
      
      fisher.table[cluster, "X1"] <- fisher.x[1,1]
      fisher.table[cluster, "X2"] <- fisher.x[1,2]
      fisher.table[cluster, "X3"] <- fisher.x[2,1]
      fisher.table[cluster, "X4"] <- fisher.x[2,2]
      fisher.table[cluster,  "cluster.size"] <- fisher.x[1,1] + fisher.x[2,1]
      
      
    }
    
    
    cl.ord <- read.table(paste0("Plots_PerTree_Clustering/", "/Tree", trees[i],"/Clustering/", "CL",nr.clusters[i],"Kmeans_ClustOrder.txt"), sep="\t", header=FALSE)
    
    cl.ord <- as.numeric(gsub("CL","" ,cl.ord[,2]))
    
    
    fisher.table <- data.frame(cluster=1:nr.clusters[i], fisher.table)
    
    write.table(fisher.table[cl.ord,], paste0(out.dir,"/",  go.name , "_" , name.out, "_Fisher_pvalues.xls"), sep="\t", row.names=F)
    
    
  }
  
  fisher.pvalues.org <- fisher.pvalues  
  fisher.pvalues <- fisher.pvalues.org
  
  fisher.pvalues.adj <- p.adjust(p=fisher.pvalues, method = "BH")
  fisher.pvalues.adj <- matrix(fisher.pvalues.adj, dim(fisher.pvalues)[1], dim(fisher.pvalues)[2])
  colnames(fisher.pvalues.adj) <- paste0( "CL", 1:nr.clusters[i])
  
  fisher.pvalues.adj <- data.frame(gene.set=gene.sets.all, fisher.pvalues.adj[, cl.ord])
  fisher.pvalues <- data.frame(gene.set=gene.sets.all, fisher.pvalues[, cl.ord])
  
  write.table(fisher.pvalues, paste0(out.dir,"/", go.name , "_" ,"ALL_PVALUES_Fisher_pvalues.xls"), sep="\t", row.names=F)
  
  write.table(fisher.pvalues.adj, paste0(out.dir,"/", go.name , "_" ,"ALL_PVALUES_Fisher_pvalues_ADJ.xls"), sep="\t", row.names=F)
  
  
}






















##############################
### MolEcol clusters
##############################


setwd("/home/Shared/data/seq/Shimizu_RNA_seq/")
load("Shimizu_workspace.Rdata")

dir.out <- "Fisher_test/MolEcol_clusters"

dir.create(dir.out, showWarnings=F, recursive=T)


# go.name <- "Kmeans_all"
# whole.clustering <- read.table("Clustering/Kmeans_all/CPM_norm_all_kmeans_table.xls", sep="\t", header=T)
# nr.clusters <- 32


go.name <- "Kmeans_drought_trees"
whole.clustering <- read.table("Clustering/Kmeans_drought_trees/CPM_norm_drought_kmeans_table.xls", sep="\t", header=T)
nr.clusters <- 22


whole.clustering <- whole.clustering.org[!is.na(whole.clustering.org$AT_ID), c("ID", "AT_ID", paste0("clustering", nr.clusters))]
colnames(whole.clustering) <- c("ID", "AT_ID","clustering")


MolEcol.clustered.genes <- read.table("Data/genes_descr_control/gene_list_in_all_clusters_DEgenes_Mol_Ecol.csv", sep=";", header=T)


fisher.pvalues <- matrix(0, max(MolEcol.clustered.genes[,2]), nr.clusters)
colnames(fisher.pvalues) <- paste0( "CL", 1:nr.clusters)

for(i in 1:max(MolEcol.clustered.genes[,2])){
  # i=2
  
  gene.set <- MolEcol.clustered.genes[MolEcol.clustered.genes[,2]==i,1]
  name.out <- paste0("MolEcol_cluster", i)
  
  
  fisher.table <- matrix(0,nr.clusters,6)
  colnames(fisher.table) <- c("X1", "X2", "X3", "X4", "pvalues", "cluster.size")
  
  for(cluster in 1:nr.clusters){
    
    # cluster <- 1
    
    clustering.tmp <- whole.clustering
    
    clustering.tmp$gene.set <- ifelse(whole.clustering$AT_ID %in%  gene.set, 1, 2)
    clustering.tmp$clustering <- ifelse(whole.clustering$clustering==cluster, 1, 2)
    
    
    fisher.x <- as.matrix(with(clustering.tmp, table(gene.set, clustering)))
    
    ### !!! very inportant to set up alternative hipothesis
    fisher.out <- fisher.test(fisher.x, alternative="greater")
    
    fisher.table[cluster, "pvalues"] <- fisher.out$p.value
    fisher.pvalues[i, cluster] <- fisher.out$p.value
    
    if(any(fisher.x <= 15)){
      fisher.table[cluster, "pvalues"] <- NA
      fisher.pvalues[i, cluster] <- NA
    }
    
    fisher.table[cluster, "X1"] <- fisher.x[1,1]
    fisher.table[cluster, "X2"] <- fisher.x[1,2]
    fisher.table[cluster, "X3"] <- fisher.x[2,1]
    fisher.table[cluster, "X4"] <- fisher.x[2,2]
    fisher.table[cluster,  "cluster.size"] <- fisher.x[1,1] + fisher.x[2,1]
    
    
  }
  
  fisher.table <- data.frame(cluster=1:nr.clusters, fisher.table)
  
  
  write.table(fisher.table, paste0(dir.out,"/", name.out, "_", go.name, "_Fisher_pvalues.xls"), sep="\t", row.names=F)
  
  
}

fisher.pvalues.org <- fisher.pvalues

fisher.pvalues <- fisher.pvalues.org


fisher.pvalues.adj <- p.adjust(p=fisher.pvalues, method = "BH")
fisher.pvalues.adj <- matrix(fisher.pvalues.adj, dim(fisher.pvalues)[1], dim(fisher.pvalues)[2])
colnames(fisher.pvalues.adj) <- paste0( "CL", 1:nr.clusters)

fisher.pvalues.adj <- data.frame(gene.set=paste0("MolEcol_cluster", 1:max(MolEcol.clustered.genes[,2])), fisher.pvalues.adj)

fisher.pvalues <- data.frame(gene.set=paste0("MolEcol_cluster", 1:max(MolEcol.clustered.genes[,2])), fisher.pvalues)

write.table(fisher.pvalues, paste0(dir.out,"/_ALL_PVALUES_", go.name, "_Fisher_pvalues.xls"), sep="\t", row.names=F)

write.table(fisher.pvalues.adj, paste0(dir.out,"/_ALL_PVALUES_", go.name, "_Fisher_pvalues_ADJ.xls"), sep="\t", row.names=F)



##############################
### MolEcol clusters: Sh.l <-> Sh.b
##############################


setwd("/home/Shared/data/seq/Shimizu_RNA_seq/")
load("Shimizu_workspace.Rdata")

dir.out <- "Fisher_test/MolEcol_clusters_ShL_ShB"

dir.create(dir.out, showWarnings=F, recursive=T)


# go.name <- "Kmeans_all_unigenesONLY"
# whole.clustering.org <- read.table("Clustering/Kmeans_all/CPM_norm_all_kmeans_table.xls", sep="\t", header=T)
# nr.clusters <- 32


go.name <- "Kmeans_drought_trees_unigenesONLY"
whole.clustering.org <- read.table("Clustering/Kmeans_drought_trees/CPM_norm_drought_kmeans_table.xls", sep="\t", header=T)
nr.clusters <- 22


whole.clustering <- whole.clustering.org[!is.na(whole.clustering.org$ID), c("ID", "AT_ID", paste0("clustering", nr.clusters))]
colnames(whole.clustering) <- c("ID", "AT_ID","clustering")
head(whole.clustering)


MolEcol.clustered.genes <- read.table("Data/genes_descr_control/gene_list_in_all_clusters_unigenes_DEgenes_Mol_Ecol.csv", header=T, sep=",")
head(MolEcol.clustered.genes)


unigene.match <- read.table("Data/genes_descr_control/reciprocal_best_hit_Sb_vs_Sl.txt", sep="\t")
head(unigene.match)


sum(whole.clustering[, 1] %in% unigene.match[,2])


whole.clustering2 <- merge(whole.clustering, unigene.match, by.x=1, by.y=2, all.x=TRUE)
names(whole.clustering2)[4] <- "unigene_ID"

sum(!is.na(whole.clustering2[,4]))


# by AT_ID
sum(whole.clustering2[, 2] %in% MolEcol.clustered.genes[, 2])

length(unique(whole.clustering2[whole.clustering2[, 2] %in% MolEcol.clustered.genes[, 2], 2]))

# by unigene
sum(whole.clustering2[, 4] %in% MolEcol.clustered.genes[, 1])

length(unique(whole.clustering2[whole.clustering2[, 4] %in% MolEcol.clustered.genes[, 1], 4]))

###


whole.clustering2 <- whole.clustering2[!is.na(whole.clustering2$unigene_ID),]


fisher.pvalues <- matrix(0, max(MolEcol.clustered.genes[,3]), nr.clusters)
colnames(fisher.pvalues) <- paste0( "CL", 1:nr.clusters)


for(i in 1:max(MolEcol.clustered.genes[,3])){
  # i=1
  
  gene.set <- MolEcol.clustered.genes[MolEcol.clustered.genes[,3]==i,1]
  name.out <- paste0("MolEcol_cluster", i)
  
  
  fisher.table <- matrix(0,nr.clusters,6)
  colnames(fisher.table) <- c("X1", "X2", "X3", "X4", "pvalues", "cluster.size")
  
  for(cluster in 1:nr.clusters){
    
    # cluster <- 1
    
    clustering.tmp <- whole.clustering2
    
    clustering.tmp$gene.set <- ifelse(whole.clustering2$unigene_ID %in%  gene.set, 1, 2)
    clustering.tmp$clustering <- ifelse(whole.clustering2$clustering==cluster, 1, 2)
    
    
    fisher.x <- as.matrix(with(clustering.tmp, table(gene.set, clustering)))
    
    ### !!! very inportant to set up alternative hipothesis
    fisher.out <- fisher.test(fisher.x, alternative="greater")
    
    fisher.table[cluster, "pvalues"] <- fisher.out$p.value
    fisher.pvalues[i, cluster] <- fisher.out$p.value
    
    if(any(fisher.x <= 10)){
      fisher.table[cluster, "pvalues"] <- NA
      fisher.pvalues[i, cluster] <- NA
    }
    
    fisher.table[cluster, "X1"] <- fisher.x[1,1]
    fisher.table[cluster, "X2"] <- fisher.x[1,2]
    fisher.table[cluster, "X3"] <- fisher.x[2,1]
    fisher.table[cluster, "X4"] <- fisher.x[2,2]
    fisher.table[cluster,  "cluster.size"] <- fisher.x[1,1] + fisher.x[2,1]
    
    
  }
  
  fisher.table <- data.frame(cluster=1:nr.clusters, fisher.table)
  
  
  write.table(fisher.table, paste0(dir.out,"/", name.out, "_", go.name, "_Fisher_pvalues.xls"), sep="\t", row.names=F)
  
  
}

fisher.pvalues.org <- fisher.pvalues

fisher.pvalues <- fisher.pvalues.org


fisher.pvalues.adj <- p.adjust(p=fisher.pvalues, method = "BH")
fisher.pvalues.adj <- matrix(fisher.pvalues.adj, dim(fisher.pvalues)[1], dim(fisher.pvalues)[2])
colnames(fisher.pvalues.adj) <- paste0( "CL", 1:nr.clusters)

fisher.pvalues.adj <- data.frame(gene.set=paste0("MolEcol_cluster", 1:max(MolEcol.clustered.genes[,3])), fisher.pvalues.adj)

fisher.pvalues <- data.frame(gene.set=paste0("MolEcol_cluster", 1:max(MolEcol.clustered.genes[,3])), fisher.pvalues)

write.table(fisher.pvalues, paste0(dir.out,"/_ALL_PVALUES_", go.name, "_Fisher_pvalues.xls"), sep="\t", row.names=F)

write.table(fisher.pvalues.adj, paste0(dir.out,"/_ALL_PVALUES_", go.name, "_Fisher_pvalues_ADJ.xls"), sep="\t", row.names=F)


