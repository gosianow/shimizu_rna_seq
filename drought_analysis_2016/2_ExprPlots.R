##############################################################################
## <<2_ExprPlots.R>>

# BioC 3.3
# Created 23 June 2016
# Updated 13 July 2016

### plots of expression for flowering genes

##############################################################################
Sys.time()
##############################################################################

library(edgeR)
library(ggplot2)
library(pheatmap)
library(limma)
library(DESeq2)
library(vsn)

##############################################################################
# Test arguments
##############################################################################

# rwd='/home/Shared/data/seq/shimizu_rna_seq'
# rcode='/home/gosia/R/shimizu_rna_seq/drought_analysis_2016'
# data_dir='Data_2016'
# analysis_dir='Analysis_2016-06-22'
# out_dir='Plots_of_flowering_genes'
# normalize_counts_function <- "normalize_counts.R"

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

### Samples info
samps_new <- read.table(paste0(data_dir, "/Samples/new_samps_interpolation.xls"), header=TRUE, sep="\t", stringsAsFactors=FALSE)

samps_new <- samps_new[order(samps_new$drough_control, samps_new$tree_id, samps_new$time_nr), ]

samps2keep <- samps_new$developmental_stage == "leaf_bud"

samps_new <- samps_new[samps2keep, ]

### Convert into factor bcs it is needed for design matrix
samps_new$tree_id <- factor(samps_new$tree_id)
samps_new$tree_id <- relevel(samps_new$tree_id, ref="8266")



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



### Load trees order data for the legends and colors in plots

trees_order <- read.table(paste0(data_dir, "/Samples/trees.xls"), header=TRUE, sep="\t", stringsAsFactors=FALSE)


### Counts
x <- read.table(paste0(data_dir, "/Data_2016-06-22/orig/Daromatica_drought_experiment_count_data_by_HTseq.csv"), sep=",", header = TRUE, as.is = TRUE)

rownames(x) <- x$gene

x <- x[, samps_new$sample_name]


### Gene description

gene_description <- read.table(paste0(data_dir, "/Data_2016-06-22/genes_description_full.xls"), sep="\t", header = TRUE, as.is = TRUE)


### Flowering-related genes

flowering_genes <- gene_description[!is.na(gene_description$flowering), ]

# extract more relevant A.t. gene symbol from at_description_flowering
flowering_genes$at_description_flowering_short <- strsplit2(flowering_genes$at_description_flowering, " \\(|/")[, 1]


# flowering_genes[grep("FT", flowering_genes$at_symbol), ]
# 
# flowering_genes[grep("SVP", flowering_genes$at_symbol), ]


### Drought-related genes

drought_up_genes <- gene_description[!is.na(gene_description$drought_regulation) & gene_description$drought_regulation == "Up", ]


drought_down_genes <- gene_description[!is.na(gene_description$drought_regulation) & gene_description$drought_regulation == "Down", ]


#############################################################################
### plots of expression for interesting gene sets from raw data
#############################################################################



plot_expr_wrapper <- function(expr, gene_set, samps_new, trees_order, out_dir_plots, out_prefix, out_suffix, ylabel){
  
  dir.create(out_dir_plots, recursive = TRUE, showWarnings = FALSE)

  for(i in 1:nrow(gene_set)){
    # i = 1

    gene_id <- gene_set$sl_id[i]
    
    if(gene_id %in% rownames(expr)){
      
      cat(paste(i, ", "))
      
      ggdf <- data.frame(expression = expr[gene_id, ], tree_legend = samps_new$tree_legend, time_date = as.Date(samps_new$time_ch, "%Y-%m-%d"), stringsAsFactors = FALSE)
      ggdf$tree_legend <- factor(ggdf$tree_legend, levels = trees_order$tree_legend)
      
      ggp <- ggplot(ggdf, aes(x = time_date, y = expression, group = tree_legend, color = tree_legend)) +
        geom_line(size = 0.2, linetype = 2) + 
        geom_point(size = 2) +
        xlab("Time") +
        ylab(ylabel) +
        ggtitle(paste0(gene_set$at_id[i], " - ", gene_id, "\n", gene_set$at_description[i])) + 
        theme_bw() +
        scale_x_date(date_labels = "%d %b %Y", date_breaks = "1 month") +
        scale_color_manual(values = trees_order$tree_col)
      
      
      pdf(paste0(out_dir_plots , "/", out_prefix, gene_set$at_description_short[i], "_" ,flowering_genes$at_id[i], "_", gene_set$sl_id[i], out_suffix, ".pdf"), h=5, w=12)
      print(ggp)
      dev.off()
      
    }
    
  }
  
}



### Plots of expression for flowering genes in CPMs
### Calculate CPMs

d <- DGEList(x, group = samps_new$tree_id)
d <- calcNormFactors(d)

dcpm <- cpm(d, normalized.lib.sizes=TRUE)
genes2keep <- rowSums( dcpm > 1 ) > 3
expr <- dcpm[genes2keep, ]


### plots for flowering genes

gene_set <- flowering_genes[, c("sl_id", "at_id", "at_description_flowering", "at_description_flowering_short")]
colnames(gene_set) <- c("sl_id", "at_id", "at_description", "at_description_short")

plot_expr_wrapper(expr = expr, gene_set = gene_set, samps_new = samps_new, trees_order = trees_order, out_dir_plots = paste0(out_dir, "/", "flowering_genes"), out_prefix = "flowering_genes_", out_suffix = "_cpm", ylabel = "Gene expression in cpm")



### plots for flowering genes: FT and SVP

gene_set <- flowering_genes[, c("sl_id", "at_id", "at_description_flowering", "at_description_flowering_short")]
colnames(gene_set) <- c("sl_id", "at_id", "at_description", "at_description_short")

gene_set <- gene_set[gene_set$at_description_short %in% c("FT", "SVP"), ]


plot_expr_wrapper(expr = expr, gene_set = gene_set, samps_new = samps_new, trees_order = trees_order, out_dir_plots = paste0(out_dir, "/", "flowering_genes_ft_svp"), out_prefix = "flowering_genes_", out_suffix = "_cpm", ylabel = "Gene expression in cpm")




### Plots of expression for flowering genes in normalized counts from DESeq2
### Do filtering based on edgeR CPM

d <- DGEList(x, group = samps_new$tree_id)
d <- calcNormFactors(d)

dcpm <- cpm(d, normalized.lib.sizes=TRUE)
genes2keep <- rowSums( dcpm > 1 ) > 3

### Create DESeq2 object

coldata <- data.frame(tree_id = samps_new$tree_id, time_ch = samps_new$time_ch)
rownames(coldata) <- colnames(x)

dds <- DESeqDataSetFromMatrix(countData = x, colData = coldata, design = ~ tree_id)

dds <- dds[genes2keep, ]

dds <- estimateSizeFactors(dds)

expr <- counts(dds, normalized = TRUE)


### plots for flowering genes

gene_set <- flowering_genes[, c("sl_id", "at_id", "at_description_flowering", "at_description_flowering_short")]
colnames(gene_set) <- c("sl_id", "at_id", "at_description", "at_description_short")

plot_expr_wrapper(expr = expr, gene_set = gene_set, samps_new = samps_new, trees_order = trees_order, out_dir_plots = paste0(out_dir, "/", "flowering_genes"), out_prefix = "flowering_genes_", out_suffix = "_deseq2_norm_counts", ylabel = "Gene expression in DESeq2 normalized counts")



### plots for flowering genes: FT and SVP

gene_set <- flowering_genes[, c("sl_id", "at_id", "at_description_flowering", "at_description_flowering_short")]
colnames(gene_set) <- c("sl_id", "at_id", "at_description", "at_description_short")

gene_set <- gene_set[gene_set$at_description_short %in% c("FT", "SVP"), ]


plot_expr_wrapper(expr = expr, gene_set = gene_set, samps_new = samps_new, trees_order = trees_order, out_dir_plots = paste0(out_dir, "/", "flowering_genes_ft_svp"), out_prefix = "flowering_genes_", out_suffix = "_deseq2_norm_counts", ylabel = "Gene expression in DESeq2 normalized counts")




### Plots of expression for flowering genes after the rlog transform from DESeq2

### rlog transformation

rld_param <- DESeq2::rlog(dds, blind = TRUE, fitType = "parametric")

out_name <- "deseq2_rld_param"
pdf(paste0(out_dir, "/", out_name, ".pdf"), w=10, h=10)
meanSdPlot(rld_param_expr)
dev.off()

expr <- assay(rld_param)


### plots for flowering genes

gene_set <- flowering_genes[, c("sl_id", "at_id", "at_description_flowering", "at_description_flowering_short")]
colnames(gene_set) <- c("sl_id", "at_id", "at_description", "at_description_short")

plot_expr_wrapper(expr = expr, gene_set = gene_set, samps_new = samps_new, trees_order = trees_order, out_dir_plots = paste0(out_dir, "/", "flowering_genes"), out_prefix = "flowering_genes_", out_suffix = "_deseq2_rlog", ylabel = "Gene expression after the rlog transformation")



### plots for flowering genes: FT and SVP

gene_set <- flowering_genes[, c("sl_id", "at_id", "at_description_flowering", "at_description_flowering_short")]
colnames(gene_set) <- c("sl_id", "at_id", "at_description", "at_description_short")

gene_set <- gene_set[gene_set$at_description_short %in% c("FT", "SVP"), ]


plot_expr_wrapper(expr = expr, gene_set = gene_set, samps_new = samps_new, trees_order = trees_order, out_dir_plots = paste0(out_dir, "/", "flowering_genes_ft_svp"), out_prefix = "flowering_genes_", out_suffix = "_deseq2_rlog", ylabel = "Gene expression after the rlog transformation")





#############################################################################
### Plot gene expression in a heatmap, per tree (1099, 1377, 970, 8266)
#############################################################################

### Use rlog transformed expression from DESeq2

expr <- assay(rld_param)

which_in <- flowering_genes$sl_id %in% rownames(expr)

expr_fl <- expr[flowering_genes$sl_id[which_in], ]



expr <- normalize.counts(expr_fl, norm.method = "mean0")
expr[expr < -2 ] <- -2
expr[expr > 2 ] <- 2


trees <- c("1099", "1377", "970", "8266")

for(i in 1:length(trees)){
  # i = 1
  
  samps_new_tmp <- samps_new[samps_new$tree_id == trees[i], ]
  expr_tmp <- expr[, samps_new_tmp$sample_name]
  
  annotation_col <- samps_new_tmp[, "water_potential", drop = FALSE]
  rownames(annotation_col) <- colnames(expr_tmp)
  annotation_colors = list(water_potential = samps_new_tmp$water_potential_col)
  
  labels_row <- flowering_genes$at_description_flowering_short
  labels_col <- samps_new_tmp$sample_name_short
  
  cluster_rows <- hclust(dist(expr_tmp), method = "average")
  
  
  pheatmap(mat = expr_tmp, cluster_cols = FALSE, cluster_rows = cluster_rows, border_color = NA, annotation_col = annotation_col, annotation_colors = annotation_colors, annotation_legend = FALSE, labels_col = labels_col, labels_row = labels_row, breaks = seq(from = -2, to = 2, length.out = 101), legend_breaks = seq(from = -2, to = 2, by = 1), fontsize_row = 6, fontsize_col = 6, fontsize = 6, filename = paste0(out_dir , "/flowering_genes_heatmap_", trees[i] , "_deseq2_rlog_mean0.pdf"), width = 7, height = 20)
  
  
}


### Comment: when genes are not DE then using norm makes the unrelevant differences look relevant...

# expr <- normalize.counts(expr_fl, norm.method = "norm")
# expr[expr < -3 ] <- -3
# expr[expr > 3 ] <- 3

# trees <- c("1099", "1377", "970", "8266")

# for(i in 1:length(trees)){
#   # i = 1
  
#   samps_new_tmp <- samps_new[samps_new$tree_id == trees[i], ]
#   expr_tmp <- expr[, samps_new_tmp$sample_name]
  
#   annotation_col <- samps_new_tmp[, "water_potential", drop = FALSE]
#   rownames(annotation_col) <- colnames(expr_tmp)
#   annotation_colors = list(water_potential = samps_new_tmp$water_potential_col)
  
#   labels_row <- flowering_genes$at_description_flowering_short
#   labels_col <- samps_new_tmp$sample_name_short
  
#   cluster_rows <- hclust(dist(expr_tmp), method = "average")
  
#   pheatmap(mat = expr_tmp, cluster_cols = FALSE, cluster_rows = cluster_rows, border_color = NA, annotation_col = annotation_col, annotation_colors = annotation_colors, annotation_legend = FALSE, labels_col = labels_col, labels_row = labels_row, breaks = seq(from = -3, to = 3, length.out = 101), legend_breaks = seq(from = -3, to = 3, by = 1), fontsize_row = 6, fontsize_col = 6, fontsize = 6, filename = paste0(out_dir , "/flowering_genes_heatmap_", trees[i] , "_deseq2_rlog_norm.pdf"), width = 7, height = 20)
  
  
# }




#############################################################################
### Plot variables interpolated at the time of taken samples
#############################################################################

plot_vars <- c("water_potential", "soil_moisture")

for(v in plot_vars){
  
  ggdf <- data.frame(expression = samps_new[, v], tree_legend = samps_new$tree_legend, time_date = as.Date(samps_new$time_ch, "%Y-%m-%d"), stringsAsFactors = FALSE)
  ggdf$tree_legend <- factor(ggdf$tree_legend, levels = trees_order$tree_legend)
  
  ggp <- ggplot(ggdf, aes(x = time_date, y = expression, group = tree_legend, color = tree_legend)) +
    geom_line(size = 0.5) + 
    geom_point(size = 2) +
    xlab("Time") +
    ylab(v) +
    theme_bw() +
    scale_x_date(date_labels = "%d %b %Y", date_breaks = "1 month") +
    scale_color_manual(values = trees_order$tree_col)
  
  
  pdf(paste0(out_dir , "/variables_", v ,".pdf"), h=5, w=12)
  print(ggp)
  dev.off()
  
  
}




