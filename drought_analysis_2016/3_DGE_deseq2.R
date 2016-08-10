##############################################################################
## <<3_DGE_deseq2.R>>

# BioC 3.3
# Created 29 June 2016
# Updated 13 July 2016 


# -----------------------------------------------------------------------------
Sys.time()
# -----------------------------------------------------------------------------

library(limma)
library(edgeR)
library(DESeq2)
library(ggplot2)
library(BiocParallel)
library(vsn)
library(plyr)


# -----------------------------------------------------------------------------
# Test arguments
# -----------------------------------------------------------------------------

rwd='/home/Shared/data/seq/shimizu_rna_seq'
rcode='/home/gosia/R/shimizu_rna_seq/drought_analysis_2016'
data_dir='Data_2016'
out_dir='Analysis_2016-06-22/Results_DGE'
de_model='between_trees_any'

# -----------------------------------------------------------------------------
# Read in the arguments
# -----------------------------------------------------------------------------

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(args)

# -----------------------------------------------------------------------------

setwd(rwd)

dir.create(out_dir, recursive = TRUE)

de_method <- "deseq2"

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
### Analysis with DESeq2
# -----------------------------------------------------------------------------

# -----------------------------------------
### DGE between trees, per tree pair, no time covariate
# -----------------------------------------

# de_model <- "between_trees_per_tree_no_time"

if(de_model == "between_trees_per_tree_no_time"){
  
  
  out_dir_tmp <- paste0(de_model, "/", de_method)
  dir.create(paste0(out_dir, "/", out_dir_tmp), recursive = TRUE)
  
  
  ### Filtering based on edgeR because the DESeq2 filtering sugested in the vignette is too relaxed and the fitted dispersion-mean trend is bad
  ### I use 4 because there are 4 different trees that we keep in the analysis and the minimal number of replicates is 6 -> min(4, 6) = 4 to have more relaxed filtering
  table(samps_new$tree_id)
  
  d <- DGEList(x, group = samps_new$tree_id)
  d <- calcNormFactors(d)
  
  cps <- cpm(d, normalized.lib.sizes = TRUE)
  genes2keep <- rowSums( cps > 1 ) >= 4
  
  table(genes2keep)
  
  ### Run DESeq2 pipeline 
  
  coldata <- data.frame(tree_id = samps_new$tree_id, time_ch = samps_new$time_ch)
  rownames(coldata) <- colnames(x)
  
  table(coldata$tree_id)
  
  model.matrix(~ tree_id, coldata)
  
  
  
  dds <- DESeqDataSetFromMatrix(countData = x, colData = coldata, design = ~ tree_id)
  
  
  dds <- dds[genes2keep, ]
  
  
  dds <- DESeq(dds, parallel = TRUE, BPPARAM = MulticoreParam(workers = 10))
  
  
  
  # dds <- estimateSizeFactors(dds)
  # 
  # dds <- estimateDispersions(dds, fitType = "parametric")
  # 
  # dds <- nbinomWaldTest(dds)
  
  
  out_name <- "deseq2_dispersion"
  pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
  plotDispEsts( dds )
  dev.off()
  
  
  
  ### rlog transformation
  
  
  rld_param <- DESeq2::rlog(dds, blind = TRUE, fitType = "parametric")
  
  rld_param_expr <- assay(rld_param)
  
  
  out_name <- "deseq2_rld_param"
  pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
  meanSdPlot(rld_param_expr)
  dev.off()
  
  
  ### Get DE results for the pairs of trees defined in contrasts
  
  resultsNames(dds)
  
  contrast_names <- c("8266vs970", "8266vs1099", "8266vs1377")
  contrast_list <- list(c("tree_id", "8266", "970"), c("tree_id", "8266", "1099"), c("tree_id", "8266", "1377"))
  names(contrast_list) <- contrast_names
  
  
  for(i in 1:length(contrast_list)){
    # i = 1
    
    
    res <- results(dds, contrast = contrast_list[[i]], alpha = 0.05)
    
    head(res)
    
    summary(res)
    
    table(res$padj < 0.05)
    
    
    
    out_name <- paste0("plot_hist_pval_", contrast_names[i])
    pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
    hist(res$pvalue, breaks = 100)
    dev.off()
    
    
    out_name <- paste0("plotMA_", contrast_names[i])
    pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
    plotMA(res)
    dev.off()
    
    
    
    resMLE <- results(dds, addMLE = TRUE, alpha = 0.05)
    
    out_name <- paste0("plotMA_unshrunken_LFC_", contrast_names[i])
    pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
    plotMA(resMLE, MLE = TRUE, main="unshrunken LFC")
    dev.off()
    
    
    
    ### Choose top significant genes to plot
    top_genes <- order(res$padj, decreasing = FALSE)
    
    dir.create(paste0(out_dir, "/", out_dir_tmp, "/", "gene_expr_", contrast_names[i]), recursive = TRUE)
    
    for(j in 1:20){
      
      gene_nr <- top_genes[j]
      gene_id <- rownames(res)[gene_nr]
      adj_pvalue <- sprintf( "%.02e", res[gene_nr, "padj"])
      
      
      ### Plot raw gene expression
      
      d <- plotCounts(dds, gene = gene_id, intgroup = c("tree_id", "time_ch"), returnData=TRUE)
      
      d$time_ch <- as.Date(d$time_ch, "%Y-%m-%d")
      d$tree_id <- factor(d$tree_id, levels = trees_order$tree_id)
      
      
      ggp <- ggplot(d, aes(x = time_ch, y = count, group = tree_id, color = tree_id)) +
        geom_line(size = 0.2, linetype = 2) + 
        geom_point(size = 2) +
        xlab("Time") +
        ylab("Gene expression in normalized counts") +
        ggtitle(paste0(gene_id, " adj_pavalue = ", adj_pvalue)) + 
        theme_bw() +
        scale_x_date(date_labels = "%d %b %Y", date_breaks = "1 month") +
        scale_color_manual(values = trees_order$tree_col)
      
      
      out_name <- paste0("gene_expr_", contrast_names[i], "/", "gene_expr_raw_", contrast_names[i], "_", j, "_", gene_id)
      pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=12, h=5)
      print(ggp)
      dev.off()
      
      
      
      ### Plot rlog transformed gene expression
      
      d <- coldata
      d$count <- rld_param_expr[gene_id, ]
      
      d$time_ch <- as.Date(d$time_ch, "%Y-%m-%d")
      d$tree_id <- factor(d$tree_id, levels = trees_order$tree_id)
      
      
      ggp <- ggplot(d, aes(x = time_ch, y = count, group = tree_id, color = tree_id)) +
        geom_line(size = 0.2, linetype = 2) + 
        geom_point(size = 2) +
        xlab("Time") +
        ylab("Gene expression in rlog transformed counts") +
        ggtitle(paste0(gene_id, " adj_pavalue = ", adj_pvalue)) + 
        theme_bw() +
        scale_x_date(date_labels = "%d %b %Y", date_breaks = "1 month") +
        scale_color_manual(values = trees_order$tree_col)
      
      
      out_name <- paste0("gene_expr_", contrast_names[i], "/", "gene_expr_rlog_", contrast_names[i], "_", j, "_", gene_id)
      pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=12, h=5)
      print(ggp)
      dev.off()
      
      
      
    }
    
    
    ### Prepare table with results to save
    
    res_out <- as.data.frame(res)
    
    res_out$gene_id <- rownames(res_out)
    
    
    out_name <- paste0(de_method, "_", de_model, "_", contrast_names[i], "_results")
    write.table(res_out, paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
    
    
    
  }
  
  
  
}








# -----------------------------------------
### DGE between trees, per tree pair, with time covariate
# -----------------------------------------

# de_model <- "between_trees_per_tree"

if(de_model == "between_trees_per_tree"){
  
  out_dir_tmp <- paste0(de_model, "/", de_method)
  dir.create(paste0(out_dir, "/", out_dir_tmp), recursive = TRUE)
  
  
  ### Filtering based on edgeR because the DESeq2 filtering sugested in the vignette is too relaxed and the fitted dispersion-mean trend is bad
  ### I use 4 because there are 4 different trees that we keep in the analysis and the minimal number of replicates is 6 -> min(4, 6) = 4 to have more relaxed filtering
  table(samps_new$tree_id)
  
  d <- DGEList(x, group = samps_new$tree_id)
  d <- calcNormFactors(d)
  
  cps <- cpm(d, normalized.lib.sizes = TRUE)
  genes2keep <- rowSums( cps > 1 ) >= 4
  
  table(genes2keep)
  
  ### Run DESeq2 pipeline 
  
  coldata <- data.frame(tree_id = samps_new$tree_id, time_ch = samps_new$time_ch, time_group = samps_new$time_group)
  rownames(coldata) <- colnames(x)
  
  table(coldata$tree_id)
  
  model.matrix(~ time_group + tree_id, coldata)
  
  
  dds <- DESeqDataSetFromMatrix(countData = x, colData = coldata, design = ~ time_group + tree_id)
  
  
  dds <- dds[genes2keep, ]
  
  
  dds <- DESeq(dds, parallel = TRUE, BPPARAM = MulticoreParam(workers = 10))
  
  
  
  # dds <- estimateSizeFactors(dds)
  # 
  # dds <- estimateDispersions(dds, fitType = "parametric")
  # 
  # dds <- nbinomWaldTest(dds)
  
  ### Plot dispersion 
  out_name <- "deseq2_dispersion"
  pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
  plotDispEsts( dds )
  dev.off()
  
  
  
  ### rlog transformation
  
  rld_param <- DESeq2::rlog(dds, blind = TRUE, fitType = "parametric")
  
  rld_param_expr <- assay(rld_param)
  
  out_name <- "deseq2_rld_param"
  pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
  meanSdPlot(rld_param_expr)
  dev.off()
  
  
  ### Get DE results for the pairs of trees
  
  resultsNames(dds)
  
  contrast_names <- c("8266vs970", "8266vs1099", "8266vs1377")
  contrast_list <- list(c("tree_id", "8266", "970"), c("tree_id", "8266", "1099"), c("tree_id", "8266", "1377"))
  names(contrast_list) <- contrast_names
  
  
  for(i in 1:length(contrast_list)){
    # i = 1
    
    
    res <- results(dds, contrast = contrast_list[[i]], alpha = 0.05)
    
    head(res)
    
    summary(res)
    
    table(res$padj < 0.05)
    
    
    
    out_name <- paste0("plot_hist_pval_", contrast_names[i])
    pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
    hist(res$pvalue, breaks = 100)
    dev.off()
    
    
    out_name <- paste0("plotMA_", contrast_names[i])
    pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
    plotMA(res)
    dev.off()
    
    
    
    resMLE <- results(dds, addMLE=TRUE, alpha = 0.05)
    
    out_name <- paste0("plotMA_unshrunken_LFC_", contrast_names[i])
    pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
    plotMA(resMLE, MLE = TRUE, main="unshrunken LFC")
    dev.off()
    
    
    
    ### Choose top significant genes to plot
    top_genes <- order(res$padj, decreasing = FALSE)
    
    dir.create(paste0(out_dir, "/", out_dir_tmp, "/", "gene_expr_", contrast_names[i]), recursive = TRUE)
    
    for(j in 1:20){
      
      gene_nr <- top_genes[j]
      gene_id <- rownames(res)[gene_nr]
      adj_pvalue <- sprintf( "%.02e", res[gene_nr, "padj"])
      
      
      ### Plot raw gene expression
      
      d <- plotCounts(dds, gene = gene_id, intgroup = c("tree_id", "time_ch"), returnData=TRUE)
      
      d$time_ch <- as.Date(d$time_ch, "%Y-%m-%d")
      d$tree_id <- factor(d$tree_id, levels = trees_order$tree_id)
      
      
      ggp <- ggplot(d, aes(x = time_ch, y = count, group = tree_id, color = tree_id)) +
        geom_line(size = 0.2, linetype = 2) + 
        geom_point(size = 2) +
        xlab("Time") +
        ylab("Gene expression in normalized counts") +
        ggtitle(paste0(gene_id, " adj_pavalue = ", adj_pvalue)) + 
        theme_bw() +
        scale_x_date(date_labels = "%d %b %Y", date_breaks = "1 month") +
        scale_color_manual(values = trees_order$tree_col)
      
      
      out_name <- paste0("gene_expr_", contrast_names[i], "/", "gene_expr_raw_", contrast_names[i], "_", j, "_", gene_id)
      pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=12, h=5)
      print(ggp)
      dev.off()
      
      
      
      ### Plot rlog transformed gene expression
      
      d <- coldata
      d$count <- rld_param_expr[gene_id, ]
      
      d$time_ch <- as.Date(d$time_ch, "%Y-%m-%d")
      d$tree_id <- factor(d$tree_id, levels = trees_order$tree_id)
      
      
      ggp <- ggplot(d, aes(x = time_ch, y = count, group = tree_id, color = tree_id)) +
        geom_line(size = 0.2, linetype = 2) + 
        geom_point(size = 2) +
        xlab("Time") +
        ylab("Gene expression in rlog transformed counts") +
        ggtitle(paste0(gene_id, " adj_pavalue = ", adj_pvalue)) + 
        theme_bw() +
        scale_x_date(date_labels = "%d %b %Y", date_breaks = "1 month") +
        scale_color_manual(values = trees_order$tree_col)
      
      
      out_name <- paste0("gene_expr_", contrast_names[i], "/", "gene_expr_rlog_", contrast_names[i], "_", j, "_", gene_id)
      pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=12, h=5)
      print(ggp)
      dev.off()
      
      
      
    }
    
    
    ### Prepare table with results to save
    
    res_out <- as.data.frame(res)
    
    res_out$gene_id <- rownames(res_out)
    
    out_name <- paste0(de_method, "_", de_model, "_", contrast_names[i], "_results")
    write.table(res_out, paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
    
    
    
  }
  
  
  
  
}




# -----------------------------------------
### DGE interaction between any time and any tree
# -----------------------------------------



if(de_model == "interaction_any_time_any_tree"){
  
  # de_model <- "interaction_any_time_any_tree"
  
  out_dir_tmp <- paste0(de_model, "/", de_method)
  dir.create(paste0(out_dir, "/", out_dir_tmp), recursive = TRUE)
  
  
  ### Filtering based on edgeR because the DESeq2 filtering sugested in the vignette is too relaxed and the fitted dispersion-mean trend is bad
  ### I use 4 because there are 4 different trees that we keep in the analysis and the minimal number of replicates is 6 -> min(4, 6) = 4 to have more relaxed filtering
  table(samps_new$tree_id)
  
  d <- DGEList(x, group = samps_new$tree_id)
  d <- calcNormFactors(d)
  
  cps <- cpm(d, normalized.lib.sizes = TRUE)
  genes2keep <- rowSums( cps > 1 ) >= 4
  
  table(genes2keep)
  
  ### Run DESeq2 pipeline 
  
  coldata <- data.frame(tree_id = samps_new$tree_id, time_ch = samps_new$time_ch, time_group = samps_new$time_group)
  rownames(coldata) <- colnames(x)
  
  
  dds <- DESeqDataSetFromMatrix(countData = x, colData = coldata, design = ~ tree_id)
  
  dds <- dds[genes2keep, ]
  
  
  
  ###### Approach 1 with specified model matrices
  
  ### This matrix contains columns with only 0s because there are combinations of factors with no corresponding samples
  modelMatrix <- model.matrix(~ time_group + tree_id + time_group:tree_id, coldata)
  
  qr(modelMatrix)$rank
  ncol(modelMatrix)
  
  ### Remove zero columns
  modelMatrix_full <- modelMatrix[, colSums(modelMatrix) > 0]
  
  qr(modelMatrix_full)$rank
  ncol(modelMatrix_full)
  
  
  modelMatrix_reduced <- model.matrix(~ time_group + tree_id, coldata)
  
  qr(modelMatrix_reduced)$rank
  ncol(modelMatrix_reduced)
  
  
  ### Use LR to test for the effect of any of the time-tree interactions
  dds <- DESeq(dds, test="LRT", full = modelMatrix_full, reduced = modelMatrix_reduced, parallel = TRUE, BPPARAM = MulticoreParam(workers = 10))
  
  
  ###### Approach 2 with defined model formulas - Does NOT work becaue "Levels or combinations of levels without any samples have resulted in column(s) of zeros in the model matrix."
  
  # design(dds) <- ~ time_group + tree_id + time_group:tree_id
  # dds <- DESeq(dds, test="LRT", reduced = ~ time_group + tree_id, parallel = TRUE, BPPARAM = MulticoreParam(workers = 10))
  
  
  ### Plot dispersion 
  out_name <- "deseq2_dispersion"
  pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
  plotDispEsts( dds )
  dev.off()
  
  ### rlog transformation
  
  rld_param <- DESeq2::rlog(dds, blind = TRUE, fitType = "parametric")
  
  rld_param_expr <- assay(rld_param)
  
  out_name <- "deseq2_rld_param"
  pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
  meanSdPlot(rld_param_expr)
  dev.off()
  
  ### Get results
  resultsNames(dds)
  
  res <- results(dds, alpha = 0.05)
  
  head(res)
  
  summary(res)
  
  table(res$padj < 0.05)
  
  
  
  out_name <- paste0("plot_hist_pval")
  pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
  hist(res$pvalue, breaks = 100)
  dev.off()
  
  
  
  ### Choose top significant genes to plot
  top_genes <- order(res$padj, decreasing = FALSE)
  
  dir.create(paste0(out_dir, "/", out_dir_tmp, "/", "gene_expr"), recursive = TRUE)
  
  for(j in 1:20){
    
    gene_nr <- top_genes[j]
    gene_id <- rownames(res)[gene_nr]
    adj_pvalue <- sprintf( "%.02e", res[gene_nr, "padj"])
    
    
    ### Plot raw gene expression
    
    d <- plotCounts(dds, gene = gene_id, intgroup = c("tree_id", "time_ch"), returnData=TRUE)
    
    d$time_ch <- as.Date(d$time_ch, "%Y-%m-%d")
    d$tree_id <- factor(d$tree_id, levels = trees_order$tree_id)
    
    
    ggp <- ggplot(d, aes(x = time_ch, y = count, group = tree_id, color = tree_id)) +
      geom_line(size = 0.2, linetype = 2) + 
      geom_point(size = 2) +
      xlab("Time") +
      ylab("Gene expression in normalized counts") +
      ggtitle(paste0(gene_id, " adj_pavalue = ", adj_pvalue)) + 
      theme_bw() +
      scale_x_date(date_labels = "%d %b %Y", date_breaks = "1 month") +
      scale_color_manual(values = trees_order$tree_col)
    
    
    out_name <- paste0("gene_expr", "/", "gene_expr_raw_", j, "_", gene_id)
    pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=12, h=5)
    print(ggp)
    dev.off()
    
    
    
    ### Plot rlog transformed gene expression
    
    d <- coldata
    d$count <- rld_param_expr[gene_id, ]
    
    d$time_ch <- as.Date(d$time_ch, "%Y-%m-%d")
    d$tree_id <- factor(d$tree_id, levels = trees_order$tree_id)
    
    
    ggp <- ggplot(d, aes(x = time_ch, y = count, group = tree_id, color = tree_id)) +
      geom_line(size = 0.2, linetype = 2) + 
      geom_point(size = 2) +
      xlab("Time") +
      ylab("Gene expression in rlog transformed counts") +
      ggtitle(paste0(gene_id, " adj_pavalue = ", adj_pvalue)) + 
      theme_bw() +
      scale_x_date(date_labels = "%d %b %Y", date_breaks = "1 month") +
      scale_color_manual(values = trees_order$tree_col)
    
    
    out_name <- paste0("gene_expr", "/", "gene_expr_rlog_", j, "_", gene_id)
    pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=12, h=5)
    print(ggp)
    dev.off()
    
    
    
  }
  
  
  ### Prepare table with results to save
  
  res_out <- as.data.frame(res)
  
  res_out$gene_id <- rownames(res_out)
  
  out_name <- paste0(de_method, "_", de_model, "_results")
  write.table(res_out, paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  
  
  
  
}



# -----------------------------------------
### DGE interaction between time, per tree - genes that are DE between pair of trees at a given time
### Corresponds to full ~ tree_id + time_group:tree_id
# -----------------------------------------



if(de_model == "interaction_per_time_per_tree"){
  
  # de_model <- "interaction_per_time_per_tree"
  
  out_dir_tmp <- paste0(de_model, "/", de_method)
  dir.create(paste0(out_dir, "/", out_dir_tmp), recursive = TRUE)
  
  
  ### Filtering based on edgeR because the DESeq2 filtering sugested in the vignette is too relaxed and the fitted dispersion-mean trend is bad
  ### I use 4 because there are 4 different trees that we keep in the analysis and the minimal number of replicates is 6 -> min(4, 6) = 4 to have more relaxed filtering
  table(samps_new$tree_id)
  
  d <- DGEList(x, group = samps_new$tree_id)
  d <- calcNormFactors(d)
  
  cps <- cpm(d, normalized.lib.sizes = TRUE)
  genes2keep <- rowSums( cps > 1 ) >= 4
  
  table(genes2keep)
  
  ### Run DESeq2 pipeline 
  
  coldata <- data.frame(tree_id = samps_new$tree_id, time_ch = samps_new$time_ch, time_group = samps_new$time_group)
  rownames(coldata) <- colnames(x)
  
  
  dds <- DESeqDataSetFromMatrix(countData = x, colData = coldata, design = ~ tree_id)
  
  dds <- dds[genes2keep, ]
  
  
  
  ###### Approach 1 with specified model matrices
  
  ### This matrix contains columns with only 0s because there are combinations of factors with no corresponding samples
  modelMatrix <- model.matrix(~ time_group + tree_id + time_group:tree_id, coldata)
  
  qr(modelMatrix)$rank
  ncol(modelMatrix)
  
  ### Remove zero columns
  modelMatrix_full <- modelMatrix[, colSums(modelMatrix) > 0]
  
  qr(modelMatrix_full)$rank
  ncol(modelMatrix_full)
  

  dds <- DESeq(dds, full = modelMatrix_full, parallel = TRUE, BPPARAM = MulticoreParam(workers = 10))
  
  
  ### Plot dispersion 
  out_name <- "deseq2_dispersion"
  pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
  plotDispEsts( dds )
  dev.off()
  
  ### rlog transformation
  
  rld_param <- DESeq2::rlog(dds, blind = TRUE, fitType = "parametric")
  
  rld_param_expr <- assay(rld_param)
  
  out_name <- "deseq2_rld_param"
  pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
  meanSdPlot(rld_param_expr)
  dev.off()
  
  
  
  ### Get DE results for the pairs of trees
  
  rn <- resultsNames(dds)
  cl <- list(rn[grep(".tree_id970", rn)], rn[grep(".tree_id1099", rn)], rn[grep(".tree_id1377", rn)])
  cl_time <- unlist(cl)
  cl_time <- strsplit2(cl_time, "time_group|.tree")[,2]
  
  contrast_names <- paste0(rep(c("8266vs970", "8266vs1099", "8266vs1377"), times = sapply(cl, length)), "_", cl_time)
  
  
  ### each contrast is for a tree effect at a given time 
  contrast_list <- data.frame(c1 = rep(c("tree_id970", "tree_id1099", "tree_id1377"), times = sapply(cl, length)), c2 = unlist(cl), stringsAsFactors = FALSE)

  
  
  for(i in 1:nrow(contrast_list)){
    # i = 1
  
    res <- results(dds, contrast = list(as.character(contrast_list[i, ])), alpha = 0.05)
    
    head(res)
    
    summary(res)
    
    table(res$padj < 0.05)
    

    out_name <- paste0("plot_hist_pval_", contrast_names[i])
    pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
    hist(res$pvalue, breaks = 100)
    dev.off()
    
    
    ### Choose top significant genes to plot
    top_genes <- order(res$padj, decreasing = FALSE)
    
    dir.create(paste0(out_dir, "/", out_dir_tmp, "/", "gene_expr_", contrast_names[i]), recursive = TRUE)
    
    for(j in 1:20){
      
      gene_nr <- top_genes[j]
      gene_id <- rownames(res)[gene_nr]
      adj_pvalue <- sprintf( "%.02e", res[gene_nr, "padj"])
      
      
      ### Plot raw gene expression
      
      d <- plotCounts(dds, gene = gene_id, intgroup = c("tree_id", "time_ch"), returnData=TRUE)
      
      d$time_ch <- as.Date(d$time_ch, "%Y-%m-%d")
      d$tree_id <- factor(d$tree_id, levels = trees_order$tree_id)
      
      
      ggp <- ggplot(d, aes(x = time_ch, y = count, group = tree_id, color = tree_id)) +
        geom_line(size = 0.2, linetype = 2) + 
        geom_point(size = 2) +
        xlab("Time") +
        ylab("Gene expression in normalized counts") +
        ggtitle(paste0(gene_id, " adj_pavalue = ", adj_pvalue)) + 
        theme_bw() +
        scale_x_date(date_labels = "%d %b %Y", date_breaks = "1 month") +
        scale_color_manual(values = trees_order$tree_col)
      
      
      out_name <- paste0("gene_expr_", contrast_names[i], "/", "gene_expr_raw_", contrast_names[i], "_", j, "_", gene_id)
      pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=12, h=5)
      print(ggp)
      dev.off()
      
      
      
      ### Plot rlog transformed gene expression

      d <- coldata
      d$count <- rld_param_expr[gene_id, ]

      d$time_ch <- as.Date(d$time_ch, "%Y-%m-%d")
      d$tree_id <- factor(d$tree_id, levels = trees_order$tree_id)


      ggp <- ggplot(d, aes(x = time_ch, y = count, group = tree_id, color = tree_id)) +
        geom_line(size = 0.2, linetype = 2) +
        geom_point(size = 2) +
        xlab("Time") +
        ylab("Gene expression in rlog transformed counts") +
        ggtitle(paste0(gene_id, " adj_pavalue = ", adj_pvalue)) +
        theme_bw() +
        scale_x_date(date_labels = "%d %b %Y", date_breaks = "1 month") +
        scale_color_manual(values = trees_order$tree_col)


      out_name <- paste0("gene_expr_", contrast_names[i], "/", "gene_expr_rlog_", contrast_names[i], "_", j, "_", gene_id)
      pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=12, h=5)
      print(ggp)
      dev.off()
      
      
      
    }
    
    
    ### Prepare table with results to save
    
    res_out <- as.data.frame(res)
    
    res_out$gene_id <- rownames(res_out)
    
    out_name <- paste0(de_method, "_", de_model, "_", contrast_names[i], "_results")
    write.table(res_out, paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
    
    
    
  }
  
  
  
  
}





# -----------------------------------------
### DGE between trees
# -----------------------------------------



if(de_model == "between_trees_any_tree"){
  
  # de_model <- "between_trees_any_tree"
  
  out_dir_tmp <- paste0(de_model, "/", de_method)
  dir.create(paste0(out_dir, "/", out_dir_tmp), recursive = TRUE)
  
  
  ### Filtering based on edgeR because the DESeq2 filtering sugested in the vignette is too relaxed and the fitted dispersion-mean trend is bad
  ### I use 4 because there are 4 different trees that we keep in the analysis and the minimal number of replicates is 6 -> min(4, 6) = 4 to have more relaxed filtering
  table(samps_new$tree_id)
  
  d <- DGEList(x, group = samps_new$tree_id)
  d <- calcNormFactors(d)
  
  cps <- cpm(d, normalized.lib.sizes = TRUE)
  genes2keep <- rowSums( cps > 1 ) >= 4
  
  table(genes2keep)
  
  ### Run DESeq2 pipeline 
  
  coldata <- data.frame(tree_id = samps_new$tree_id, time_ch = samps_new$time_ch, time_group = samps_new$time_group)
  rownames(coldata) <- colnames(x)
  
  
  dds <- DESeqDataSetFromMatrix(countData = x, colData = coldata, design = ~ tree_id)
  
  dds <- dds[genes2keep, ]
  
  
  ###### Approach 2 with defined model formulas
  
  design(dds) <- ~ time_group + tree_id
  dds <- DESeq(dds, test = "LRT", reduced = ~ time_group, parallel = TRUE, BPPARAM = MulticoreParam(workers = 10))
  
  
  ### Plot dispersion 
  out_name <- "deseq2_dispersion"
  pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
  plotDispEsts( dds )
  dev.off()
  
  ### rlog transformation
  
  rld_param <- DESeq2::rlog(dds, blind = TRUE, fitType = "parametric")
  
  rld_param_expr <- assay(rld_param)
  
  out_name <- "deseq2_rld_param"
  pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
  meanSdPlot(rld_param_expr)
  dev.off()
  
  ### Get results
  resultsNames(dds)
  
  res <- results(dds, alpha = 0.05)
  
  head(res)
  
  summary(res)
  
  table(res$padj < 0.05)
  
  
  
  out_name <- paste0("plot_hist_pval")
  pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
  hist(res$pvalue, breaks = 100)
  dev.off()
  
  
  
  ### Choose top significant genes to plot
  top_genes <- order(res$padj, decreasing = FALSE)
  
  dir.create(paste0(out_dir, "/", out_dir_tmp, "/", "gene_expr"), recursive = TRUE)
  
  for(j in 1:20){
    
    gene_nr <- top_genes[j]
    gene_id <- rownames(res)[gene_nr]
    adj_pvalue <- sprintf( "%.02e", res[gene_nr, "padj"])
    
    
    ### Plot raw gene expression
    
    d <- plotCounts(dds, gene = gene_id, intgroup = c("tree_id", "time_ch"), returnData=TRUE)
    
    d$time_ch <- as.Date(d$time_ch, "%Y-%m-%d")
    d$tree_id <- factor(d$tree_id, levels = trees_order$tree_id)
    
    
    ggp <- ggplot(d, aes(x = time_ch, y = count, group = tree_id, color = tree_id)) +
      geom_line(size = 0.2, linetype = 2) + 
      geom_point(size = 2) +
      xlab("Time") +
      ylab("Gene expression in normalized counts") +
      ggtitle(paste0(gene_id, " adj_pavalue = ", adj_pvalue)) + 
      theme_bw() +
      scale_x_date(date_labels = "%d %b %Y", date_breaks = "1 month") +
      scale_color_manual(values = trees_order$tree_col)
    
    
    out_name <- paste0("gene_expr", "/", "gene_expr_raw_", j, "_", gene_id)
    pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=12, h=5)
    print(ggp)
    dev.off()
    
    
    
    ### Plot rlog transformed gene expression
    
    d <- coldata
    d$count <- rld_param_expr[gene_id, ]
    
    d$time_ch <- as.Date(d$time_ch, "%Y-%m-%d")
    d$tree_id <- factor(d$tree_id, levels = trees_order$tree_id)
    
    
    ggp <- ggplot(d, aes(x = time_ch, y = count, group = tree_id, color = tree_id)) +
      geom_line(size = 0.2, linetype = 2) + 
      geom_point(size = 2) +
      xlab("Time") +
      ylab("Gene expression in rlog transformed counts") +
      ggtitle(paste0(gene_id, " adj_pavalue = ", adj_pvalue)) + 
      theme_bw() +
      scale_x_date(date_labels = "%d %b %Y", date_breaks = "1 month") +
      scale_color_manual(values = trees_order$tree_col)
    
    
    out_name <- paste0("gene_expr", "/", "gene_expr_rlog_", j, "_", gene_id)
    pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=12, h=5)
    print(ggp)
    dev.off()
    
    
    
  }
  
  
  ### Prepare table with results to save
  
  res_out <- as.data.frame(res)
  
  res_out$gene_id <- rownames(res_out)
  
  out_name <- paste0(de_method, "_", de_model, "_results")
  write.table(res_out, paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  
  
}



# -----------------------------------------
### DGE between time
# -----------------------------------------



if(de_model == "between_time_all_trees"){
  
  # de_model <- "between_time_all_trees"
  
  out_dir_tmp <- paste0(de_model, "/", de_method)
  dir.create(paste0(out_dir, "/", out_dir_tmp), recursive = TRUE)
  
  
  ### Filtering based on edgeR because the DESeq2 filtering sugested in the vignette is too relaxed and the fitted dispersion-mean trend is bad
  ### I use 4 because there are 4 different trees that we keep in the analysis and the minimal number of replicates is 6 -> min(4, 6) = 4 to have more relaxed filtering
  table(samps_new$tree_id)
  
  d <- DGEList(x, group = samps_new$tree_id)
  d <- calcNormFactors(d)
  
  cps <- cpm(d, normalized.lib.sizes = TRUE)
  genes2keep <- rowSums( cps > 1 ) >= 4
  
  table(genes2keep)
  
  ### Run DESeq2 pipeline 
  
  coldata <- data.frame(tree_id = samps_new$tree_id, time_ch = samps_new$time_ch, time_group = samps_new$time_group)
  rownames(coldata) <- colnames(x)
  
  
  dds <- DESeqDataSetFromMatrix(countData = x, colData = coldata, design = ~ tree_id)
  
  dds <- dds[genes2keep, ]
  
  
  ###### Approach 2 with defined model formulas
  
  design(dds) <- ~ time_group + tree_id
  dds <- DESeq(dds, test = "LRT", reduced = ~ tree_id, parallel = TRUE, BPPARAM = MulticoreParam(workers = 10))
  
  
  out_name <- "deseq2_dds"
  save(dds, file = paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".Rdata"))
  
  
  ### Plot dispersion 
  out_name <- "deseq2_dispersion"
  pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
  plotDispEsts( dds )
  dev.off()
  
  ### rlog transformation
  
  rld_param <- DESeq2::rlog(dds, blind = TRUE, fitType = "parametric")
  
  rld_param_expr <- assay(rld_param)
  
  out_name <- "deseq2_rld_param"
  pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
  meanSdPlot(rld_param_expr)
  dev.off()
  
  ### Get results
  
  resultsNames(dds)
  
  res <- results(dds, alpha = 0.05)
  
  head(res)
  
  summary(res)
  
  table(res$padj < 0.05)
  
  
  
  out_name <- paste0("plot_hist_pval")
  pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
  hist(res$pvalue, breaks = 100)
  dev.off()
  
  
  
  ### Choose top significant genes to plot
  top_genes <- order(res$padj, decreasing = FALSE)
  
  dir.create(paste0(out_dir, "/", out_dir_tmp, "/", "gene_expr"), recursive = TRUE)
  
  for(j in 1:20){
    
    gene_nr <- top_genes[j]
    gene_id <- rownames(res)[gene_nr]
    adj_pvalue <- sprintf( "%.02e", res[gene_nr, "padj"])
    
    
    ### Plot raw gene expression
    
    d <- plotCounts(dds, gene = gene_id, intgroup = c("tree_id", "time_ch"), returnData=TRUE)
    
    d$time_ch <- as.Date(d$time_ch, "%Y-%m-%d")
    d$tree_id <- factor(d$tree_id, levels = trees_order$tree_id)
    
    
    ggp <- ggplot(d, aes(x = time_ch, y = count, group = tree_id, color = tree_id)) +
      geom_line(size = 0.2, linetype = 2) + 
      geom_point(size = 2) +
      xlab("Time") +
      ylab("Gene expression in normalized counts") +
      ggtitle(paste0(gene_id, " adj_pavalue = ", adj_pvalue)) + 
      theme_bw() +
      scale_x_date(date_labels = "%d %b %Y", date_breaks = "1 month") +
      scale_color_manual(values = trees_order$tree_col)
    
    
    out_name <- paste0("gene_expr", "/", "gene_expr_raw_", j, "_", gene_id)
    pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=12, h=5)
    print(ggp)
    dev.off()
    
    
    
    ### Plot rlog transformed gene expression
    
    d <- coldata
    d$count <- rld_param_expr[gene_id, ]
    
    d$time_ch <- as.Date(d$time_ch, "%Y-%m-%d")
    d$tree_id <- factor(d$tree_id, levels = trees_order$tree_id)
    
    
    ggp <- ggplot(d, aes(x = time_ch, y = count, group = tree_id, color = tree_id)) +
      geom_line(size = 0.2, linetype = 2) + 
      geom_point(size = 2) +
      xlab("Time") +
      ylab("Gene expression in rlog transformed counts") +
      ggtitle(paste0(gene_id, " adj_pavalue = ", adj_pvalue)) + 
      theme_bw() +
      scale_x_date(date_labels = "%d %b %Y", date_breaks = "1 month") +
      scale_color_manual(values = trees_order$tree_col)
    
    
    out_name <- paste0("gene_expr", "/", "gene_expr_rlog_", j, "_", gene_id)
    pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=12, h=5)
    print(ggp)
    dev.off()
    
    
    
  }
  
  
  ### Prepare table with results to save
  
  res_out <- as.data.frame(res)
  
  res_out$gene_id <- rownames(res_out)
  
  out_name <- paste0(de_method, "_", de_model, "_results")
  write.table(res_out, paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  
  
  
  
}






# -----------------------------------------
### DGE intercept
# -----------------------------------------



if(de_model == "intercept"){
  
  # de_model <- "intercept"
  
  out_dir_tmp <- paste0(de_model, "/", de_method)
  dir.create(paste0(out_dir, "/", out_dir_tmp), recursive = TRUE)
  
  
  ### Filtering based on edgeR because the DESeq2 filtering sugested in the vignette is too relaxed and the fitted dispersion-mean trend is bad
  ### I use 4 because there are 4 different trees that we keep in the analysis and the minimal number of replicates is 6 -> min(4, 6) = 4 to have more relaxed filtering
  table(samps_new$tree_id)
  
  d <- DGEList(x, group = samps_new$tree_id)
  d <- calcNormFactors(d)
  
  cps <- cpm(d, normalized.lib.sizes = TRUE)
  genes2keep <- rowSums( cps > 1 ) >= 4
  
  table(genes2keep)
  
  ### Run DESeq2 pipeline 
  
  coldata <- data.frame(tree_id = samps_new$tree_id, time_ch = samps_new$time_ch, time_group = samps_new$time_group)
  rownames(coldata) <- colnames(x)
  
  
  dds <- DESeqDataSetFromMatrix(countData = x, colData = coldata, design = ~ tree_id)
  
  dds <- dds[genes2keep, ]
  
  
  ###### Approach 2 with defined model formulas 
  
  design(dds) <- ~ time_group + tree_id
  dds <- DESeq(dds, test = "LRT", reduced = ~ 1, parallel = TRUE, BPPARAM = MulticoreParam(workers = 10))
  
  
  ### Plot dispersion 
  out_name <- "deseq2_dispersion"
  pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
  plotDispEsts( dds )
  dev.off()
  
  ### rlog transformation
  
  rld_param <- DESeq2::rlog(dds, blind = TRUE, fitType = "parametric")
  
  rld_param_expr <- assay(rld_param)
  
  out_name <- "deseq2_rld_param"
  pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
  meanSdPlot(rld_param_expr)
  dev.off()
  
  ### Get results
  
  
  resultsNames(dds)
  
  res <- results(dds, alpha = 0.05)
  
  head(res)
  
  summary(res)
  
  table(res$padj < 0.05)
  
  
  
  out_name <- paste0("plot_hist_pval")
  pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
  hist(res$pvalue, breaks = 100)
  dev.off()
  
  
  
  ### Choose top significant genes to plot
  top_genes <- order(res$padj, decreasing = FALSE)
  
  dir.create(paste0(out_dir, "/", out_dir_tmp, "/", "gene_expr"), recursive = TRUE)
  
  for(j in 1:20){
    
    gene_nr <- top_genes[j]
    gene_id <- rownames(res)[gene_nr]
    adj_pvalue <- sprintf( "%.02e", res[gene_nr, "padj"])
    
    
    ### Plot raw gene expression
    
    d <- plotCounts(dds, gene = gene_id, intgroup = c("tree_id", "time_ch"), returnData=TRUE)
    
    d$time_ch <- as.Date(d$time_ch, "%Y-%m-%d")
    d$tree_id <- factor(d$tree_id, levels = trees_order$tree_id)
    
    
    ggp <- ggplot(d, aes(x = time_ch, y = count, group = tree_id, color = tree_id)) +
      geom_line(size = 0.2, linetype = 2) + 
      geom_point(size = 2) +
      xlab("Time") +
      ylab("Gene expression in normalized counts") +
      ggtitle(paste0(gene_id, " adj_pavalue = ", adj_pvalue)) + 
      theme_bw() +
      scale_x_date(date_labels = "%d %b %Y", date_breaks = "1 month") +
      scale_color_manual(values = trees_order$tree_col)
    
    
    out_name <- paste0("gene_expr", "/", "gene_expr_raw_", j, "_", gene_id)
    pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=12, h=5)
    print(ggp)
    dev.off()
    
    
    
    ### Plot rlog transformed gene expression
    
    d <- coldata
    d$count <- rld_param_expr[gene_id, ]
    
    d$time_ch <- as.Date(d$time_ch, "%Y-%m-%d")
    d$tree_id <- factor(d$tree_id, levels = trees_order$tree_id)
    
    
    ggp <- ggplot(d, aes(x = time_ch, y = count, group = tree_id, color = tree_id)) +
      geom_line(size = 0.2, linetype = 2) + 
      geom_point(size = 2) +
      xlab("Time") +
      ylab("Gene expression in rlog transformed counts") +
      ggtitle(paste0(gene_id, " adj_pavalue = ", adj_pvalue)) + 
      theme_bw() +
      scale_x_date(date_labels = "%d %b %Y", date_breaks = "1 month") +
      scale_color_manual(values = trees_order$tree_col)
    
    
    out_name <- paste0("gene_expr", "/", "gene_expr_rlog_", j, "_", gene_id)
    pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=12, h=5)
    print(ggp)
    dev.off()
    
    
    
  }
  
  
  ### Prepare table with results to save
  
  res_out <- as.data.frame(res)
  
  res_out$gene_id <- rownames(res_out)
  
  out_name <- paste0(de_method, "_", de_model, "_results")
  write.table(res_out, paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  
  
  
  
}




# -----------------------------------------
### DGE between conditions: drought vs. control
# -----------------------------------------

### Create new tree variable as sugessted in section 3.12.1 Linear combinations in DESeq2 vignette

samps_new$tree_id_new <- samps_new$tree_id

samps_new$tree_id_new <- mapvalues(samps_new$tree_id_new, from = c("8266", "970" , "1099", "1377"), to = c("1", "2", "1", "2"))

samps_new$condition <- factor(samps_new$drough_control, levels = c("control", "drought"))




if(de_model == "between_conditions"){
  
  # de_model <- "between_conditions"
  
  out_dir_tmp <- paste0(de_model, "/", de_method)
  dir.create(paste0(out_dir, "/", out_dir_tmp), recursive = TRUE)
  
  
  ### Filtering based on edgeR because the DESeq2 filtering sugested in the vignette is too relaxed and the fitted dispersion-mean trend is bad
  ### I use 4 because there are 4 different trees that we keep in the analysis and the minimal number of replicates is 6 -> min(4, 6) = 4 to have more relaxed filtering
  table(samps_new$tree_id)
  
  d <- DGEList(x, group = samps_new$tree_id)
  d <- calcNormFactors(d)
  
  cps <- cpm(d, normalized.lib.sizes = TRUE)
  genes2keep <- rowSums( cps > 1 ) >= 4
  
  table(genes2keep)
  
  ### Run DESeq2 pipeline 
  
  coldata <- data.frame(tree_id = samps_new$tree_id, time_ch = samps_new$time_ch, tree_id_new = samps_new$tree_id_new, time_group = samps_new$time_group, condition = samps_new$condition)
  rownames(coldata) <- colnames(x)
  
  
  dds <- DESeqDataSetFromMatrix(countData = x, colData = coldata, design = ~ 1)
  
  dds <- dds[genes2keep, ]
  

  design(dds) <- ~ condition + tree_id_new + time_group
  
  dds <- DESeq(dds, parallel = TRUE, BPPARAM = MulticoreParam(workers = 10))
  
  
  ### Plot dispersion 
  out_name <- "deseq2_dispersion"
  pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
  plotDispEsts( dds )
  dev.off()
  
  
  
  ### rlog transformation
  
  rld_param <- DESeq2::rlog(dds, blind = TRUE, fitType = "parametric")
  
  rld_param_expr <- assay(rld_param)
  
  out_name <- "deseq2_rld_param"
  pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
  meanSdPlot(rld_param_expr)
  dev.off()
  
  
  
  ### Get results
  resultsNames(dds)
  
  res <- results(dds, contrast = c("condition", "control", "drought"), alpha = 0.05)
  
  head(res)
  
  summary(res)
  
  table(res$padj < 0.05)
  
  
  
  out_name <- paste0("plot_hist_pval")
  pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=10, h=10)
  hist(res$pvalue, breaks = 100)
  dev.off()
  
  
  
  ### Choose top significant genes to plot
  top_genes <- order(res$padj, decreasing = FALSE)
  
  dir.create(paste0(out_dir, "/", out_dir_tmp, "/", "gene_expr"), recursive = TRUE)
  
  for(j in 1:20){
    
    gene_nr <- top_genes[j]
    gene_id <- rownames(res)[gene_nr]
    adj_pvalue <- sprintf( "%.02e", res[gene_nr, "padj"])
    
    
    ### Plot raw gene expression
    
    d <- plotCounts(dds, gene = gene_id, intgroup = c("tree_id", "time_ch"), returnData=TRUE)
    
    d$time_ch <- as.Date(d$time_ch, "%Y-%m-%d")
    d$tree_id <- factor(d$tree_id, levels = trees_order$tree_id)
    
    
    ggp <- ggplot(d, aes(x = time_ch, y = count, group = tree_id, color = tree_id)) +
      geom_line(size = 0.2, linetype = 2) + 
      geom_point(size = 2) +
      xlab("Time") +
      ylab("Gene expression in normalized counts") +
      ggtitle(paste0(gene_id, " adj_pavalue = ", adj_pvalue)) + 
      theme_bw() +
      scale_x_date(date_labels = "%d %b %Y", date_breaks = "1 month") +
      scale_color_manual(values = trees_order$tree_col)
    
    
    out_name <- paste0("gene_expr", "/", "gene_expr_raw_", j, "_", gene_id)
    pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=12, h=5)
    print(ggp)
    dev.off()
    
    
    
    ### Plot rlog transformed gene expression
    
    # d <- coldata
    # d$count <- rld_param_expr[gene_id, ]
    # 
    # d$time_ch <- as.Date(d$time_ch, "%Y-%m-%d")
    # d$tree_id <- factor(d$tree_id, levels = trees_order$tree_id)
    # 
    # 
    # ggp <- ggplot(d, aes(x = time_ch, y = count, group = tree_id, color = tree_id)) +
    #   geom_line(size = 0.2, linetype = 2) + 
    #   geom_point(size = 2) +
    #   xlab("Time") +
    #   ylab("Gene expression in rlog transformed counts") +
    #   ggtitle(paste0(gene_id, " adj_pavalue = ", adj_pvalue)) + 
    #   theme_bw() +
    #   scale_x_date(date_labels = "%d %b %Y", date_breaks = "1 month") +
    #   scale_color_manual(values = trees_order$tree_col)
    # 
    # 
    # out_name <- paste0("gene_expr", "/", "gene_expr_rlog_", j, "_", gene_id)
    # pdf(paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".pdf"), w=12, h=5)
    # print(ggp)
    # dev.off()
    
    
    
  }
  
  
  ### Prepare table with results to save
  
  res_out <- as.data.frame(res)
  
  res_out$gene_id <- rownames(res_out)
  
  out_name <- paste0(de_method, "_", de_model, "_results")
  write.table(res_out, paste0(out_dir, "/", out_dir_tmp, "/", out_name, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  
  
}






###########################
# 3_DGE_deseq2 done!
###########################