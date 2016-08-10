##############################################################################
## <<2_MDS.R>>

# BioC 3.3
# Created 23 June 2016
# Updated 26 June 2016 

### crude analysis - MDS plots

##############################################################################
Sys.time()
##############################################################################

library(edgeR)
library(EDASeq)
library(ggplot2)
library(ggdendro)
library(DESeq)

##############################################################################
# Test arguments
##############################################################################

rwd='/home/Shared/data/seq/shimizu_rna_seq'
rcode='/home/gosia/R/shimizu_rna_seq/drought_analysis_2016'
data_dir='Data_2016'
analysis_dir='Analysis_2016-06-22'
out_dir='Plots_MDS'


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

out_dir <- paste0(analysis_dir, "/", out_dir)
dir.create(out_dir, recursive = TRUE)


####################################################
### Load data
####################################################

### Sample info
samps_new <- read.table(paste0(data_dir, "/Samples/new_samps_interpolation.xls"), header=TRUE, sep="\t", stringsAsFactors=FALSE)

samps_new <- samps_new[order(samps_new$drough_control, samps_new$tree_id, samps_new$time_nr), ]

# Change color for the sample from flower bud
samps_new$tree_col[samps_new$developmental_stage == "flower_bud"] <- "darkmagenta"
samps_new$sample_name_short <- paste0(samps_new$sample_id, "_", samps_new$sample_name_short)


# Make color palette with respect to time

all_days_short <- read.table(paste0(data_dir, "/Unique_days_short_nr.xls"), header=TRUE, sep="\t", stringsAsFactors=FALSE)

color_palette <- as.character( colorRampPalette(c("red", 'orange', "green", 'darkgreen'))(nrow(all_days_short)))
names(color_palette) <- all_days_short$days_nr

m <- match(samps_new$time_nr, all_days_short$days_nr)

samps_new$time_col <- color_palette[m]


# Make color palette with respect to water potential 

all_wp <- as.character(seq(min(round(samps_new$water_potential, 2)), max(round(samps_new$water_potential, 2)), 0.01))

color_palette <- as.character( colorRampPalette(c("red", 'orange', "dodgerblue", 'dodgerblue4'))(length(all_wp)))

m <- match(as.character(round(samps_new$water_potential, 2)), all_wp)

samps_new$water_potential_col <- color_palette[m]


### Counts
x <- read.table(paste0(data_dir, "/Data_2016-06-22/orig/Daromatica_drought_experiment_count_data_by_HTseq.csv"), sep=",", header = TRUE, as.is = TRUE)

rownames(x) <- x$gene

x <- x[, samps_new$sample_name]

dim(x)


### Load gene description

gene_description <- read.table(paste0(data_dir, "/Data_2016-06-22/genes_description_full.xls"), sep="\t", header = TRUE, as.is = TRUE)
head(gene_description)



### Load housekeeping genes
housekeeping_genes <- read.table(paste0(data_dir, "/GeneControlSets/Housekeeping_Genes/Table4_Arabidopsis_Housekeepieng_Gene_ID.txt"), sep="\t", header = TRUE, as.is = TRUE)

colnames(housekeeping_genes) <- "at_id"
housekeeping_genes$at_id <- toupper(housekeeping_genes$at_id)



### Load moderate drought response genes



####################################################
### MDS plots
####################################################


plot_mds_wrapper <- function(d, samps_new, out_dir, out_name){
  # Function that does all the plots 
  
  
  ### FC
  pdf(paste0(out_dir, "/", out_name, "mds_logfc_500_raw.pdf"), w=10, h=10)
  mdsb <- plotMDS(d, col = samps_new$tree_col, top=500, labels=samps_new$sample_name_short)
  dev.off()
  
  
  ### BVC
  pdf(paste0(out_dir, "/", out_name, "mds_bvc_500_raw.pdf"), w=10, h=10)
  mdsb <- plotMDS(d, col = samps_new$tree_col, method="bcv", top=500, labels=samps_new$sample_name_short)
  dev.off()
  
  pdf(paste0(out_dir, "/", out_name, "mds_bvc_500_raw_col_time.pdf"), w=10, h=10)
  plotMDS(mdsb, col = samps_new$time_col, labels=samps_new$sample_name_short)
  dev.off()
  
  pdf(paste0(out_dir, "/", out_name, "mds_bvc_500_raw_col_wp.pdf"), w=10, h=10)
  plotMDS(mdsb, col = samps_new$water_potential_col, labels=samps_new$sample_name_short)
  dev.off()
  
  pdf(paste0(out_dir, "/", out_name, "mds_bvc_500_adj_axes.pdf"), w=10, h=10)
  plotMDS(mdsb, col = samps_new$tree_col, labels = samps_new$sample_name_short, cex.axis = 2, xlim=c(-1, 1), ylim=c(-1, 1), cex = 1.2)
  dev.off()
  
  
  pdf(paste0(out_dir, "/", out_name, "mds_bvc_500_adj_axes_col_time.pdf"), w=10, h=10)
  plotMDS(mdsb, col = samps_new$time_col, labels = samps_new$sample_name_short, cex.axis = 2, xlim=c(-1, 1), ylim=c(-1, 1), cex = 1.2)
  dev.off()
  
  pdf(paste0(out_dir, "/", out_name, "mds_bvc_500_adj_axes_col_wp.pdf"), w=10, h=10)
  plotMDS(mdsb, col = samps_new$water_potential_col, labels = samps_new$sample_name_short, cex.axis = 2, xlim=c(-1, 1), ylim=c(-1, 1), cex = 1.2)
  dev.off()
  
  
  ### Make dengograms
  colnames(mdsb$distance.matrix) <- rownames(mdsb$distance.matrix) <- samps_new$sample_name_short
  
  hc <- hclust(as.dist(mdsb$distance.matrix), method = "single")
  
  
  ggp <- ggdendrogram(hc, size=4, theme_dendro = FALSE, rotate=TRUE) + 
    theme(axis.text.y = element_text(colour = samps_new[hc$order, "tree_col"]))
  
  pdf(paste0(out_dir, "/", out_name, "dendrogram_bvc_500.pdf"), w=10, h=10)
  print(ggp)
  dev.off()
  
  
  ggp <- ggdendrogram(hc, size=4, theme_dendro = FALSE, rotate=TRUE) + 
    theme(axis.text.y = element_text(colour = samps_new[hc$order, "time_col"]))
  
  pdf(paste0(out_dir, "/", out_name, "dendrogram_bvc_500_col_time.pdf"), w=10, h=10)
  print(ggp)
  dev.off()
  
  # PCA with EDASeq
  
  pdf(paste0(out_dir, "/", out_name, "pca_edaseq.pdf"), w=10, h=10)
  EDASeq::plotPCA(cpm(d, normalized.lib.sizes=TRUE), col = samps_new$tree_col)
  dev.off()
  
  # PCA with DESeq
  
  # cds <- newCountDataSet(countData = d$counts, conditions = samps_new$tree_id)
  # cds <- estimateSizeFactors( cds )
  # cds <- estimateDispersions( cds, method="blind" )
  # vsd <- varianceStabilizingTransformation( cds )
  # 
  # pdf(paste0(out_dir, "/", out_name, "pca_deseq.pdf"), w=10, h=10)
  # DESeq::plotPCA(vsd)
  # dev.off()
  
  ### RLE plots
  
  pdf(paste0(out_dir, "/", out_name, "rle_raw.pdf"), w=15, h=7)
  plotRLE(d$counts, col = samps_new$tree_col, outline = FALSE, main="Raw counts")
  dev.off()
  
  
  pdf(paste0(out_dir, "/", out_name, "rle_cpm.pdf"), w=15, h=7)
  plotRLE(cpm(d, normalized.lib.sizes=TRUE), col = samps_new$tree_col, outline = FALSE, main = "CPM normalization")
  dev.off()
  
}







### MDS plots for #(cpm > 1) > 10

d <- DGEList(x, group = samps_new$tree_id)
d <- calcNormFactors(d)

# make sure a gene is expressed (CPM > 1) in more samples
cps <- cpm(d, normalized.lib.sizes = TRUE)
d <- d[ rowSums( cps > 1 ) > 10, ]
dim(d$counts)


plot_mds_wrapper(d = d, samps_new = samps_new, out_dir = out_dir, out_name = "all_genes_cpm10_")




### MDS plots for #(cpm > 1) > 3

d <- DGEList(x, group = samps_new$tree_id)
d <- calcNormFactors(d)

# make sure a gene is expressed (CPM > 1) in more samples
cps <- cpm(d, normalized.lib.sizes = TRUE)
d <- d[ rowSums( cps > 1 ) > 3, ]
dim(d$counts)

plot_mds_wrapper(d = d, samps_new = samps_new, out_dir = out_dir, out_name = "all_genes_cpm3_")




### MDS plots for #(cpm > 1) > 3 for housekeeping genes

housekeeping_sl_ids <- gene_description[gene_description$at_id %in% housekeeping_genes$at_id, "sl_id"]
all(housekeeping_sl_ids %in% rownames(x))


d <- DGEList(x, group = samps_new$tree_id)
d <- calcNormFactors(d)

# make sure a gene is expressed (CPM > 1) in more samples
cps <- cpm(d, normalized.lib.sizes = TRUE)
d <- d[rowSums( cps > 1 ) > 3 & rownames(cps) %in% housekeeping_sl_ids, ]


plot_mds_wrapper(d = d, samps_new = samps_new, out_dir = out_dir, out_name = "housekeeping_genes_cpm3_")


### MDS plots for #(cpm > 1) > 3 for drought genes

table(gene_description$drought_regulation, useNA = "always")

drought_sl_ids <- gene_description[!is.na(gene_description$drought_regulation), "sl_id"]
all(drought_sl_ids %in% rownames(x))


d <- DGEList(x, group = samps_new$tree_id)
d <- calcNormFactors(d)

# make sure a gene is expressed (CPM > 1) in more samples
cps <- cpm(d, normalized.lib.sizes = TRUE)
d <- d[rowSums( cps > 1 ) > 3 & rownames(cps) %in% drought_sl_ids, ]
dim(d)



plot_mds_wrapper(d = d, samps_new = samps_new, out_dir = out_dir, out_name = "drought_genes_cpm3_")
















































