
#####################################################################################################

### Fisher's exact tests

#####################################################################################################


##############################
### Athaliana_flowering_genes
##############################


setwd("/home/Shared/data/seq/Shimizu_RNA_seq/")
load("Shimizu_workspace.Rdata")

dir.out <- "Fisher_test/Athaliana_flowering_genes"

dir.create(dir.out, showWarnings=F, recursive=T)


go.name <- "Kmeans_all"
whole.clustering.org <- read.table("Clustering/Kmeans_all/CPM_norm_all_kmeans_table.xls", sep="\t", header=T)
nr.clusters <- 32


# go.name <- "Kmeans_drought_trees"
# whole.clustering.org <- read.table("Clustering/Kmeans_drought_trees/CPM_norm_drought_kmeans_table.xls", sep="\t", header=T)
# nr.clusters <- 22



whole.clustering <- whole.clustering.org[!is.na(whole.clustering.org$AT_ID), c("ID", "AT_ID", paste0("clustering", nr.clusters))]
colnames(whole.clustering) <- c("ID", "AT_ID","clustering")


gene.set <- Athaliana.flowering.genes[,1]
name.out <- "Athaliana_flowering_genes"


fisher.pvalues <- matrix(0, 1 ,nr.clusters)
colnames(fisher.pvalues) <- paste0( "CL", 1:nr.clusters)


fisher.table <- matrix(0,nr.clusters,6)
colnames(fisher.table) <- c("X1", "X2", "X3", "X4", "pvalues", "cluster.size")

for(cluster in 1:nr.clusters){
  
  #cluster <- 1
  
  clustering.tmp <- whole.clustering
  
  clustering.tmp$gene.set <- ifelse(whole.clustering$AT_ID %in%  gene.set, 1, 2)
  clustering.tmp$clustering <- ifelse(whole.clustering$clustering==cluster, 1, 2)
  
  
  fisher.x <- as.matrix(with(clustering.tmp, table(gene.set, clustering)))
  
  fisher.out <- fisher.test(fisher.x, alternative="greater")
  
  fisher.table[cluster, "pvalues"] <- fisher.out$p.value
  fisher.pvalues[1, cluster] <- fisher.out$p.value
  
  if(any(fisher.x < 10)){    
    fisher.table[cluster, "pvalues"] <- NA
    fisher.pvalues[1, cluster] <- NA
  }

  
  fisher.table[cluster, "X1"] <- fisher.x[1,1]
  fisher.table[cluster, "X2"] <- fisher.x[1,2]
  fisher.table[cluster, "X3"] <- fisher.x[2,1]
  fisher.table[cluster, "X4"] <- fisher.x[2,2]
  fisher.table[cluster,  "cluster.size"] <- fisher.x[1,1] + fisher.x[2,1]
  
  
}

fisher.table <- data.frame(cluster=1:nr.clusters, fisher.table)

write.table(fisher.table, paste0(dir.out,"/", name.out, "_", go.name, "_Fisher_pvalues.xls"), sep="\t", row.names=F)


fisher.pvalues.org <- fisher.pvalues

fisher.pvalues <- fisher.pvalues.org


fisher.pvalues.adj <- p.adjust(p=fisher.pvalues, method = "BH")
fisher.pvalues.adj <- matrix(fisher.pvalues.adj, dim(fisher.pvalues)[1], dim(fisher.pvalues)[2])
colnames(fisher.pvalues.adj) <- paste0( "CL", 1:nr.clusters)

fisher.pvalues.adj <- data.frame(gene.set="Athaliana.flowering.genes", fisher.pvalues.adj)

fisher.pvalues <- data.frame(gene.set="Athaliana.flowering.genes", fisher.pvalues)

write.table(fisher.pvalues, paste0(dir.out,"/_ALL_PVALUES_", go.name, "_Fisher_pvalues.xls"), sep="\t", row.names=F)

write.table(fisher.pvalues.adj, paste0(dir.out,"/_ALL_PVALUES_", go.name, "_Fisher_pvalues_ADJ.xls"), sep="\t", row.names=F)



##############################
### gene sets
##############################

setwd("/home/Shared/data/seq/Shimizu_RNA_seq/")
load("Shimizu_workspace.Rdata")

dir.out <- "Fisher_test/Gene_sets"

dir.create(dir.out, showWarnings=F, recursive=T)


# go.name <- "Kmeans_all"
# whole.clustering <- read.table("Clustering/Kmeans_all/CPM_norm_all_kmeans_table.xls", sep="\t", header=T)
# nr.clusters <- 32


go.name <- "Kmeans_drought_trees"
whole.clustering <- read.table("Clustering/Kmeans_drought_trees/CPM_norm_drought_kmeans_table.xls", sep="\t", header=T)
nr.clusters <- 22


whole.clustering <- whole.clustering.org[!is.na(whole.clustering.org$AT_ID), c("ID", "AT_ID", paste0("clustering", nr.clusters))]
colnames(whole.clustering) <- c("ID", "AT_ID","clustering")


gene.sets.dir <- "Data/genes_descr_control/gene_sets_for_Gosia"
gene.sets.all <- dir(gene.sets.dir)



fisher.pvalues <- matrix(0, length(gene.sets.all),nr.clusters)
colnames(fisher.pvalues) <- paste0( "CL", 1:nr.clusters)

for(i in 1:length(gene.sets.all)){
  # i=2
  
  cat("\n", gene.sets.all[i], "\n") 
  file.data <- read.table(paste0(gene.sets.dir, "/", gene.sets.all[i]), header=T, sep=";")
  
  print(head(file.data))
  
  gene.set <- file.data[,1]
  name.out <- gsub(".csv", "",gene.sets.all[i])
  
  
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

fisher.pvalues.adj <- data.frame(gene.set=gene.sets.all, fisher.pvalues.adj)

fisher.pvalues <- data.frame(gene.set=gene.sets.all, fisher.pvalues)

write.table(fisher.pvalues, paste0(dir.out,"/_ALL_PVALUES_", go.name, "_Fisher_pvalues.xls"), sep="\t", row.names=F)

write.table(fisher.pvalues.adj, paste0(dir.out,"/_ALL_PVALUES_", go.name, "_Fisher_pvalues_ADJ.xls"), sep="\t", row.names=F)


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


