##############################################################################
## <<2_MDS.R>>

# BioC 3.3
# Created 23 June 2016

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

# rwd='/home/Shared/data/seq/shimizu_rna_seq'
# rcode='/home/gosia/R/shimizu_rna_seq/drought_analysis_2016'
# data_dir='Data_2016'
# analysis_dir='Analysis_2016-06-22'
# out_dir='Plots_MDS'


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

### Samples info
new.samps <- read.table(paste0(data_dir, "/Samples/new_samps_interpolation.xls"), header=TRUE, sep="\t", stringsAsFactors=FALSE)

new.samps <- new.samps[order(new.samps$drough.control, new.samps$tree_ID, new.samps$time_nr), ]

# Change color for the sample from flower bud
new.samps$tree_col[new.samps$developmental_stage == "flower_bud"] <- "darkmagenta"
new.samps$short.name <- paste0(new.samps$sample_ID, "_", new.samps$short.name)


# Make color palette with respect to time

all.days.short <- read.table(paste0(data_dir, "/Unique_days_short_nr.xls"), header=TRUE, sep="\t", stringsAsFactors=FALSE)

colors.palette <- as.character( colorRampPalette(c("red", 'orange', "green", 'darkgreen'))(nrow(all.days.short)))
names(colors.palette) <- all.days.short$days.nr

m <- match(new.samps$time_nr, all.days.short$days.nr)

new.samps$time_col <- colors.palette[m]


### Counts
x <- read.table(paste0(data_dir, "/Data_2016-06-22/orig/Daromatica_drought_experiment_count_data_by_HTseq.csv"), sep=",", header = TRUE, as.is = TRUE)

x <- x[, new.samps$sample_name]

dim(x)


####################################################
### MDS plots
####################################################

d <- DGEList(x, group = new.samps$tree_ID)
d <- calcNormFactors(d)

# make sure a gene is expressed (CPM > 1) in more samples
cps <- cpm(d, normalized.lib.sizes = TRUE)
d <- d[ rowSums( cps > 1 ) > 10, ]
dim(d$counts)



pdf(paste0(out_dir, "/mds_cpm10_logfc_500_raw.pdf"), w=10, h=10)
mdsb <- plotMDS(d, col = new.samps$tree_col, top=500, labels=new.samps$short.name)
dev.off()

pdf(paste0(out_dir, "/mds_cpm10_bvc_500_raw.pdf"), w=10, h=10)
mdsb <- plotMDS(d, col = new.samps$tree_col, method="bcv", top=500, labels=new.samps$short.name)
dev.off()

pdf(paste0(out_dir, "/mds_cpm10_bvc_500_adj_axes.pdf"), w=10, h=10)
plotMDS(mdsb, col = new.samps$tree_col, labels = new.samps$short.name, cex.axis = 2, xlim=c(-1, 1), ylim=c(-1, 1), cex = 1.2)
dev.off()


pdf(paste0(out_dir, "/mds_cpm10_bvc_500_adj_axes_col_time.pdf"), w=10, h=10)
plotMDS(mdsb, col = new.samps$time_col, labels = new.samps$short.name, cex.axis = 2, xlim=c(-1, 1), ylim=c(-1, 1), cex = 1.2)
dev.off()



colnames(mdsb$distance.matrix) <- rownames(mdsb$distance.matrix) <- new.samps$short.name

hc <- hclust(as.dist(mdsb$distance.matrix), method = "single")

pdf(paste0(out_dir, "/dendrogram_cpm10_bvc_500.pdf"), w=10, h=10)
ggdendrogram(hc, size=4, theme_dendro = FALSE, rotate=TRUE) + 
  theme(axis.text.y = element_text(colour = new.samps[hc$order, "tree_col"]))
dev.off()


pdf(paste0(out_dir, "/dendrogram_cpm10_bvc_500_col_time.pdf"), w=10, h=10)
ggdendrogram(hc, size=4, theme_dendro = FALSE, rotate=TRUE) + 
  theme(axis.text.y = element_text(colour = new.samps[hc$order, "time_col"]))
dev.off()


pdf(paste0(out_dir, "/pca_edaseq_cpm10.pdf"), w=10, h=10)
EDASeq::plotPCA(cpm(d, normalized.lib.sizes=TRUE), col = new.samps$tree_col)
dev.off()

# PCA with DESeq

cds <- newCountDataSet(countData = d$counts, conditions = new.samps$tree_ID)
cds <- estimateSizeFactors( cds )
cds <- estimateDispersions( cds, method="blind" )
vsd <- varianceStabilizingTransformation( cds )


pdf(paste0(out_dir, "/pca_deseq_cpm10.pdf"), w=10, h=10)
DESeq::plotPCA(vsd)
dev.off()


pdf(paste0(out_dir, "/rle_cpm10.pdf"), w=15, h=7)
plotRLE(d$counts, col = new.samps$tree_col, outline = FALSE, main="Raw counts")
plotRLE(cpm(d, normalized.lib.sizes=TRUE), col = new.samps$tree_col, outline = FALSE, main = "CPM normalization")
dev.off()


####################################################
### MDS plots for #(cpm > 1) > 3
####################################################

d <- DGEList(x, group = new.samps$tree_ID)
d <- calcNormFactors(d)

# make sure a gene is expressed (CPM > 1) in more samples
cps <- cpm(d, normalized.lib.sizes = TRUE)
d <- d[ rowSums( cps > 1 ) > 3, ]
dim(d$counts)



pdf(paste0(out_dir, "/mds_cpm3_logfc_500_raw.pdf"), w=10, h=10)
mdsb <- plotMDS(d, col = new.samps$tree_col, top=500, labels=new.samps$short.name)
dev.off()

pdf(paste0(out_dir, "/mds_cpm3_bvc_500_raw.pdf"), w=10, h=10)
mdsb <- plotMDS(d, col = new.samps$tree_col, method="bcv", top=500, labels=new.samps$short.name)
dev.off()

pdf(paste0(out_dir, "/mds_cpm3_bvc_500_adj_axes.pdf"), w=10, h=10)
plotMDS(mdsb, col = new.samps$tree_col, labels = new.samps$short.name, cex.axis = 2, xlim=c(-1, 1.2), ylim=c(-1, 1.2), cex = 1.2)
dev.off()

pdf(paste0(out_dir, "/mds_cpm3_bvc_500_adj_axes_col_time.pdf"), w=10, h=10)
plotMDS(mdsb, col = new.samps$time_col, labels = new.samps$short.name, cex.axis = 2, xlim=c(-1, 1.2), ylim=c(-1, 1.2), cex = 1.2)
dev.off()


colnames(mdsb$distance.matrix) <- rownames(mdsb$distance.matrix) <- new.samps$short.name

hc <- hclust(as.dist(mdsb$distance.matrix), method = "single")

pdf(paste0(out_dir, "/dendrogram_cpm3_bvc_500.pdf"), w=10, h=10)
ggdendrogram(hc, size=4, theme_dendro = FALSE, rotate=TRUE) + 
  theme(axis.text.y = element_text(colour = new.samps[hc$order, "tree_col"]))
dev.off()














