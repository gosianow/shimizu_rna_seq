##############################################################################
## <<3_variance_vs_mean.R>>

# BioC 3.3
# Created 28 June 2016
# Updated 11 July 2016 


##############################################################################
Sys.time()
##############################################################################

library(limma)
library(edgeR)
library(vsn)
library(ggplot2)
library(DESeq2)

##############################################################################
# Test arguments
##############################################################################

rwd='/home/Shared/data/seq/shimizu_rna_seq'
rcode='/home/gosia/R/shimizu_rna_seq/drought_analysis_2016'
data_dir='Data_2016'
analysis_dir='Analysis_2016-06-22'
out_dir='Plots_variance'
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

out_dir <- paste0(analysis_dir, "/", out_dir)
dir.create(out_dir, recursive = TRUE)

source(paste0(rcode, "/", normalize_counts_function))

####################################################
### Load data
####################################################

### Sample info
samps_new <- read.table(paste0(data_dir, "/Samples/new_samps_interpolation.xls"), header=TRUE, sep="\t", stringsAsFactors=FALSE)

samps_new <- samps_new[order(samps_new$drough_control, samps_new$tree_id, samps_new$time_nr), ]

# Change color for the sample from flower bud
samps_new$tree_col[samps_new$developmental_stage == "flower_bud"] <- "darkmagenta"
samps_new$sample_name_short <- paste0(samps_new$sample_id, "_", samps_new$sample_name_short)

### Convert into factor bcs it is needed for design matrix
samps_new$tree_id <- factor(samps_new$tree_id)
samps_new$tree_id <- relevel(samps_new$tree_id, ref="8266")

### use time as factor because otherwise it will be treated as linear but the patterns of expression do not have to be linear 
samps_new$time_ch <- factor(samps_new$time_ch)


### Counts
x <- read.table(paste0(data_dir, "/Data_2016-06-22/orig/Daromatica_drought_experiment_count_data_by_HTseq.csv"), sep=",", header = TRUE, as.is = TRUE)

rownames(x) <- x$gene

x <- x[, samps_new$sample_name]

### There are counts in this table that do not correspond to genes! Remove them!
genes2keep <- !grepl("_", rownames(x))
x <- x[genes2keep, ]


dim(x)
length(grep("g[0-9]+", rownames(x)))



####################################################
### Filtering samples
####################################################

### Do not consider 990 and 8212 and flower bud samples because they have too few replicates 

samps2keep <- grepl("970|8266|1099|1377", samps_new$tree_id) & grepl("leaf_bud", samps_new$developmental_stage)

samps_new <- samps_new[samps2keep, ]

samps_new$tree_id <- factor(samps_new$tree_id)
samps_new$time_ch <- factor(samps_new$time_ch)

x <- x[, samps2keep]



####################################################
### Variance vs mean plots - edgeR
####################################################


d <- DGEList(x, group = samps_new$tree_id)
d <- calcNormFactors(d)

table(samps_new$tree_id)


### make sure a gene is expressed (CPM > 1) in more samples
### I use 4 because there are 4 different trees that we keep in the analysis and the minimal number of replicates is 6 -> min(4, 6) = 4 to have more relaxed filtering
cps <- cpm(d, normalized.lib.sizes = TRUE)
genes2keep <- rowSums( cps > 1 ) >= 4 
d <- d[ genes2keep, ]
dim(d$counts)

cps <- cps[genes2keep, ]



### estimate dispersion for desing with trees only
design <- model.matrix(~ tree_id, data = samps_new)

disp <- estimateDisp(d, design, robust=FALSE)


out_name <- "edgeR_bcv"
pdf(paste0(out_dir, "/", out_name, ".pdf"), w=10, h=10)
plotBCV(disp)
dev.off()



fit <- glmQLFit(d, design, robust=FALSE)

out_name <- "edgeR_qldisp"
pdf(paste0(out_dir, "/", out_name, ".pdf"), w=10, h=10)
plotQLDisp(fit)
dev.off()



### estimate dispersion for desing with trees and time
design <- model.matrix(~ tree_id + time_ch, data = samps_new)

disp <- estimateDisp(d, design, robust=FALSE)


out_name <- "edgeR_trees_time_bcv"
pdf(paste0(out_dir, "/", out_name, ".pdf"), w=10, h=10)
plotBCV(disp)
dev.off()



fit <- glmQLFit(d, design, robust=FALSE)

out_name <- "edgeR_trees_time_qldisp"
pdf(paste0(out_dir, "/", out_name, ".pdf"), w=10, h=10)
plotQLDisp(fit)
dev.off()





mean <- apply(cps, 1, mean)
sd <- apply(cps, 1, sd)

ggdf <- data.frame(mean = mean, sd = sd)

ggp <- ggplot(ggdf, aes(x = log2(mean), y = sqrt(sd))) +
  geom_point(size = 0.2, alpha = 0.5) +
  xlab("log2 of mean cpm") +
  ylab("sqrt of sd of cpm")


out_name <- "edgeR_cpm"
  
pdf(paste0(out_dir, "/", out_name, ".pdf"), w=10, h=10)
print(ggp)
dev.off()



mean <- apply(log2(cps + 0.5), 1, mean)
sd <- apply(log2(cps + 0.5), 1, sd)


ggdf <- data.frame(mean = mean, sd = sd)

ggp <- ggplot(ggdf, aes(x = mean, y = sqrt(sd))) +
  geom_point(size = 0.2, alpha = 0.5) +
  xlab("mean log2(cpm + 0.5)") +
  ylab("sqrt of sd of log2(cpm + 0.5)")


out_name <- "edgeR_log2cpm"

pdf(paste0(out_dir, "/", out_name, ".pdf"), w=10, h=10)
print(ggp)
dev.off()




####################################################
### Variance vs mean plots - voom
####################################################


d <- DGEList(x, group = samps_new$tree_id)
d <- calcNormFactors(d)

# make sure a gene is expressed (CPM > 1) in more samples
cps <- cpm(d, normalized.lib.sizes = TRUE)
genes2keep <- rowSums( cps > 1 ) >= 4
d <- d[ genes2keep, ]
dim(d$counts)

cps <- cps[genes2keep, ]



### estimate dispersion for desing with trees only
design <- model.matrix(~ tree_id, data = samps_new)


out_name <- "voom_trees"

pdf(paste0(out_dir, "/", out_name, ".pdf"), w=10, h=10)
v <- voom(d, design, plot=TRUE)
dev.off()



### estimate dispersion for desing with trees and time

design <- model.matrix(~ tree_id + time_ch, data = samps_new)


out_name <- "voom_trees_time"

pdf(paste0(out_dir, "/", out_name, ".pdf"), w=10, h=10)
v <- voom(d, design, plot=TRUE)
dev.off()


####################################################
### Variance vs mean of [0,1] normalized cpms
####################################################


d <- DGEList(x, group = samps_new$tree_id)
d <- calcNormFactors(d)

# make sure a gene is expressed (CPM > 1) in more samples
cps <- cpm(d, normalized.lib.sizes = TRUE)
genes2keep <- rowSums( cps > 1 ) >= 4
d <- d[ genes2keep, ]
dim(d$counts)

cps <- cps[genes2keep, ]



mean <- apply(cps, 1, mean)
sd <- apply(normalize.counts(cps, norm.method = "01"), 1, sd)


ggdf <- data.frame(mean = mean, sd = sd)

ggp <- ggplot(ggdf, aes(x = log2(mean), y = sqrt(sd))) +
  geom_point(size = 0.2, alpha = 0.5) +
  xlab("log2 of mean cpm") +
  ylab("sqrt of sd of normalized to [0,1] cpm")


out_name <- "normalized_01"

pdf(paste0(out_dir, "/", out_name, ".pdf"), w=10, h=10)
print(ggp)
dev.off()



##############################################################################
### Variance vs mean of variance stabilized counts with DESeq2
##############################################################################

### Filtering based on edgeR because the DESeq2 filtering sugested in the vignette (dds <- dds[ rowSums(counts(dds)) > 1, ]) is too relaxed and the fitted dispersion-mean trend is bad 
### I use 4 because there are 4 different trees that we keep in the analysis and the minimal number of replicates is 6 -> min(4, 6) = 4 to have more relaxed filtering
table(samps_new$tree_id)

d <- DGEList(x, group = samps_new$tree_id)
d <- calcNormFactors(d)

cps <- cpm(d, normalized.lib.sizes = TRUE)
genes2keep <- rowSums( cps > 1 ) >= 4

table(genes2keep)


### DESeq2 pipeline

coldata <- data.frame(tree_id = samps_new$tree_id, time_ch = samps_new$time_ch)
rownames(coldata) <- colnames(x)

dds <- DESeqDataSetFromMatrix(countData = x, colData = coldata, design = ~ tree_id)

dds <- dds[genes2keep, ]

dim(dds)

# dds <- DESeq(dds) ### wrapper


dds <- estimateSizeFactors(dds)


dds <- estimateDispersions(dds, fitType = "parametric")

out_name <- "deseq2_dispersion_param"
pdf(paste0(out_dir, "/", out_name, ".pdf"), w=10, h=10)
DESeq2::plotDispEsts( dds )
dev.off()



dds <- estimateDispersions(dds, fitType = "local")

out_name <- "deseq2_dispersion_local"
pdf(paste0(out_dir, "/", out_name, ".pdf"), w=10, h=10)
DESeq2::plotDispEsts( dds )
dev.off()



# dds <- nbinomWaldTest(dds)





out_name <- "deseq2_counts"
pdf(paste0(out_dir, "/", out_name, ".pdf"), w=10, h=10)
meanSdPlot(counts(dds, normalized=TRUE))
dev.off()


out_name <- "deseq2_log"
pdf(paste0(out_dir, "/", out_name, ".pdf"), w=10, h=10)
meanSdPlot(log2(counts(dds, normalized=TRUE) + 1))
dev.off()


### varianceStabilizingTransformation


vsd_param <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE, fitType="parametric") ### returns normalized log2 counts

vsd_param_expr <- assay(vsd_param)

out_name <- "deseq2_vsd_param"
pdf(paste0(out_dir, "/", out_name, ".pdf"), w=10, h=10)
meanSdPlot(vsd_param_expr)
dev.off()



vsd_loc <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE, fitType = "local") ### returns normalized log2 counts

vsd_loc_expr <- assay(vsd_loc)

out_name <- "deseq2_vsd_loc"
pdf(paste0(out_dir, "/", out_name, ".pdf"), w=10, h=10)
meanSdPlot(vsd_loc_expr)
dev.off()



### rlog


rld_param <- DESeq2::rlog(dds, blind=TRUE, fitType = "parametric")

rld_param_expr <- assay(rld_param)


out_name <- "deseq2_rld_param"
pdf(paste0(out_dir, "/", out_name, ".pdf"), w=10, h=10)
meanSdPlot(rld_param_expr)
dev.off()




rld_loc <- DESeq2::rlog(dds, blind=TRUE, fitType = "local")

rld_loc_expr <- assay(rld_loc)


out_name <- "deseq2_rld_loc"
pdf(paste0(out_dir, "/", out_name, ".pdf"), w=10, h=10)
meanSdPlot(rld_loc_expr)
dev.off()




### vst

vsd2_param <- DESeq2::vst(dds, blind=TRUE, fitType = "parametric")

vsd2_param_expr <- assay(vsd2_param)

out_name <- "deseq2_vsd2_param"
pdf(paste0(out_dir, "/", out_name, ".pdf"), w=10, h=10)
meanSdPlot(vsd2_param_expr)
dev.off()


vsd2_loc <- DESeq2::vst(dds, blind=TRUE, fitType = "local")

vsd2_loc_expr <- assay(vsd2_loc)

out_name <- "deseq2_vsd2_loc"
pdf(paste0(out_dir, "/", out_name, ".pdf"), w=10, h=10)
meanSdPlot(vsd2_loc_expr)
dev.off()








































