##############################################################################
## <<2_ExprPlots.R>>

# BioC 3.3
# Created 23 June 2016


### plots of expression for flowering genes

##############################################################################
Sys.time()
##############################################################################

library(edgeR)
library(ggplot2)

##############################################################################
# Test arguments
##############################################################################

# rwd='/home/Shared/data/seq/shimizu_rna_seq'
# rcode='/home/gosia/R/shimizu_rna_seq/drought_analysis_2016'
# data_dir='Data_2016'
# analysis_dir='Analysis_2016-06-22'
# out_dir='Plots_of_flowering_genes'


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

samps2keep <- new.samps$developmental_stage == "leaf_bud"

new.samps <- new.samps[samps2keep, ]


###

trees.order <- read.table(paste0(data_dir, "/Samples/trees.xls"), header=TRUE, sep="\t", stringsAsFactors=FALSE)


### Counts
x <- read.table(paste0(data_dir, "/Data_2016-06-22/orig/Daromatica_drought_experiment_count_data_by_HTseq.csv"), sep=",", header = TRUE, as.is = TRUE)

rownames(x) <- x$gene

x <- x[, new.samps$sample_name]


### Gene description

genes.full.description <- read.table(paste0(data_dir, "/Data_2016-06-22/genes_description_full.xls"), sep="\t", header = TRUE, as.is = TRUE)

flowering_genes <- genes.full.description[!is.na(genes.full.description$Flowering), ]


# flowering_genes[grep("FT", flowering_genes$AT_symbol), ]
# 
# flowering_genes[grep("SVP", flowering_genes$AT_symbol), ]


####################################################
### plots of expression for flowering genes from raw data
####################################################

### Calculate CPMs

d.org <- DGEList(x, group = new.samps$tree_ID)
d.org <- calcNormFactors(d.org)

d.cpm <- cpm(d.org, normalized.lib.sizes=TRUE)

d.cpm.l <- log(d.cpm + min(d.cpm[d.cpm != 0]))



for(i in 1:nrow(flowering_genes)){
  # i = 1
  cat(paste(i, ", "))
  
  gene_id <- flowering_genes$ID[i]
  
  if(gene_id %in% rownames(x)){
    
    ggdf <- data.frame(expression = d.cpm[gene_id, ], tree_legend = new.samps$tree_legend, time_date = as.Date(new.samps$time_ch, "%Y-%m-%d"), stringsAsFactors = FALSE)
    ggdf$tree_legend <- factor(ggdf$tree_legend, levels = trees.order$tree_legend)
    
    ggp <- ggplot(ggdf, aes(x = time_date, y = expression, group = tree_legend, color = tree_legend)) +
      geom_line(size = 1) + 
      geom_point() +
      xlab("Time") +
      ylab("Gene expression in cpm") +
      ggtitle(paste0(gene_id, "\n", flowering_genes$AT_NAME[i])) + 
      theme_bw() +
      scale_x_date(date_labels = "%d %b %Y", date_breaks = "1 month") +
      scale_color_manual(values = trees.order$tree_col)
    
    
    pdf(paste0(out_dir , "/Flowering_gene_", flowering_genes$AT_ID[i], "_", flowering_genes$ID[i], ".pdf"), h=5, w=12)
    print(ggp)
    dev.off()
    
  }
  
}



### Plot variables interpolated at the time of takes samples

plot.vars <- c("Water.Potential", "Soil.Moisture")

for(v in plot.vars){
  
  ggdf <- data.frame(expression = new.samps[, v], tree_legend = new.samps$tree_legend, time_date = as.Date(new.samps$time_ch, "%Y-%m-%d"), stringsAsFactors = FALSE)
  ggdf$tree_legend <- factor(ggdf$tree_legend, levels = trees.order$tree_legend)
  
  ggp <- ggplot(ggdf, aes(x = time_date, y = expression, group = tree_legend, color = tree_legend)) +
    geom_line(size = 1) + 
    geom_point() +
    xlab("Time") +
    ylab(v) +
    theme_bw() +
    scale_x_date(date_labels = "%d %b %Y", date_breaks = "1 month") +
    scale_color_manual(values = trees.order$tree_col)
  
  
  pdf(paste0(out_dir , "/Variables_", v ,".pdf"), h=5, w=12)
  print(ggp)
  dev.off()
  
  
}




