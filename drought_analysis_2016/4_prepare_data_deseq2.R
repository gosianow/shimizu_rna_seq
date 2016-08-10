##############################################################################
## <<4_prepare_data_deseq2.R>>

# BioC 3.3
# Created 4 Aug 2016


##############################################################################
Sys.time()
##############################################################################

library(limma)
library(edgeR)
library(DESeq2)
library(vsn) # for meanSdPlot

##############################################################################
# Test arguments
##############################################################################

rwd='/home/Shared/data/seq/shimizu_rna_seq'
out_dir='Analysis_2016-06-22/4_Clustering/data'
dge_dir='Analysis_2016-06-22/3_DGE'
data_dir='Data_2016'

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

dir.create(out_dir, recursive = TRUE)

FDR <- 0.05

# -----------------------------------------------------------------------------
### Load data
# -----------------------------------------------------------------------------

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
samps_new$time_group <- factor(samps_new$time_group, levels = paste0("time", 1:11))
# Use time11 as reference because it has samples for all the trees
samps_new$time_group <- relevel(samps_new$time_group, ref = "time11")


### Counts
x <- read.table(paste0(data_dir, "/Data_2016-06-22/orig/Daromatica_drought_experiment_count_data_by_HTseq.csv"), sep=",", header = TRUE, as.is = TRUE)

rownames(x) <- x$gene

x <- x[, samps_new$sample_name]

### There are counts in this table that do not correspond to genes! Remove them!
genes2keep <- !grepl("_", rownames(x))
x <- x[genes2keep, ]


dim(x)
length(grep("g[0-9]+", rownames(x)))


### Load trees order data for the legends and colors in plots

trees_order <- read.table(paste0(data_dir, "/Samples/trees.xls"), header=TRUE, sep="\t", stringsAsFactors=FALSE)



# -----------------------------------------------------------------------------
### Filtering samples
# -----------------------------------------------------------------------------

### Do not consider 990, 8212, flower bud and 1st November 2008 samples because they have too few replicates 

samps2keep <- grepl("970|8266|1099|1377", samps_new$tree_id) & grepl("leaf_bud", samps_new$developmental_stage) & as.Date(samps_new$time_ch, "%Y-%m-%d") > as.Date("2008-11-01", "%Y-%m-%d")

samps_new <- samps_new[samps2keep, ]

samps_new$tree_id <- factor(samps_new$tree_id)
samps_new$time_ch <- factor(samps_new$time_ch)
samps_new$time_group <- factor(samps_new$time_group)


x <- x[, samps2keep]

trees_order <- trees_order[grepl("970|8266|1099|1377", trees_order$tree_id), ]



# -----------------------------------------------------------------------------
# Prepare the rlog transformed expression data
# -----------------------------------------------------------------------------

if(!file.exists(file.path(out_dir, "rld_param_expr.rda"))){
  
  ### Use the same filtering as for the DGE analysis
  d <- DGEList(x, group = samps_new$tree_id)
  d <- calcNormFactors(d)
  
  cps <- cpm(d, normalized.lib.sizes = TRUE)
  genes2keep <- rowSums( cps > 1 ) >= 4
  
  ### Prepare the DESeq2 object
  
  coldata <- data.frame(tree_id = samps_new$tree_id, time_ch = samps_new$time_ch)
  rownames(coldata) <- colnames(x)
  
  dds <- DESeqDataSetFromMatrix(countData = x, colData = coldata, design = ~ tree_id)
  
  dds <- dds[genes2keep, ]
  
  dds <- estimateSizeFactors(dds)
  
  
  ### rlog transformation
  
  rld_param <- DESeq2::rlog(dds, blind = TRUE, fitType = "parametric")
  
  rld_param_expr <- assay(rld_param)
  
  
  out_name <- "deseq2_rld_param.pdf"
  pdf(file.path(out_dir, out_name), w=10, h=10)
  meanSdPlot(rld_param_expr)
  dev.off()
  
  save(rld_param_expr, file = file.path(out_dir, "rld_param_expr.rda"))
  
  write.table(data.frame(gene_id = rownames(rld_param_expr), rld_param_expr, stringsAsFactors = FALSE), file.path(out_dir, "rld_param_expr.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  
}else{
  
  load(file.path(out_dir, "rld_param_expr.rda"))
  
}


# -----------------------------------------------------------------------------
# Prepare the clustering data - rlog gene expression
# -----------------------------------------------------------------------------


# ---------------------------------------------------
#  For deseq2 interaction_any_time_any_tree
# ---------------------------------------------------

de_method <- "deseq2"
de_model <- "interaction_any_time_any_tree"


out_name <- paste0(de_method, "_", de_model, "_results.txt")
res <- read.table(file.path(dge_dir, de_model, de_method, out_name), header = TRUE, sep = "\t", as.is = TRUE)

table(res$padj < FDR, useNA = "always")

de_genes <- res[res$padj < FDR & !is.na(res$padj), "gene_id"]


data_out <- data.frame(gene_id = de_genes, rld_param_expr[de_genes, ], stringsAsFactors = FALSE)


out_name <- paste0("expr_rlog_", de_method, "_", de_model, ".txt")
write.table(data_out, file.path(out_dir, out_name), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


# -----------------------------------------------------------------------------
# Prepare the clustering data - glm coefficients
# -----------------------------------------------------------------------------

### Scratch code to read the glm coeffs

# out_name <- "deseq2_dds"
# load(paste0(dge_dir, "/", out_dir_tmp, "/", out_name, ".Rdata"))
# 
# 
# out_name <- paste0(de_method, "_", de_model, "_results.txt")
# res <- read.table(paste0(dge_dir, "/", out_dir_tmp, "/", out_name), header = TRUE, sep = "\t", as.is = TRUE)
# 
# 
# table(res$padj < 0.05, useNA = "always")
# 
# # de_genes <- res[res$padj < 0.05 & !is.na(res$padj), "gene_id"]
# de_genes <- res[order(res$padj, decreasing = FALSE), "gene_id"]
# de_genes <- de_genes[1:100]
# 
# 
# betas <- coef(dds)
# colnames(betas)
# 
# mat <- betas[de_genes, grep("time11", colnames(betas))]
# 
# thr <- 3 
# mat[mat < -thr] <- -thr
# mat[mat > thr] <- thr
# 
# pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101), cluster_col=FALSE, filename = paste0(out_dir , "/time_coeffs.pdf"), width = 7, height = 20)



