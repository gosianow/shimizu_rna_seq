#####################################################################################################
# OLD data analysis:

# plots of expression
# GLM - simple with run.edgeR function
# run.edgeR.robust 

#####################################################################################################

#####################################################################################################

### Load data 

#####################################################################################################

setwd("/home/Shared/data/seq/Shimizu_RNA_seq/")

### samples & factors
new.samps <- read.table("Samples_out/new_samps_interpolation_November.csv", header=TRUE, sep=";", stringsAsFactors=FALSE)

rownames(new.samps) <- new.samps$sample_name
new.samps$tree_ID <- as.factor(new.samps$tree_ID)

# define param for ploting: colors, legend
new.samps$tree_col <- new.samps$tree_ID
levels(new.samps$tree_col) <- c("orange", "cyan", "green", "blue", "red", "magenta")
new.samps$tree_col <- as.character(new.samps$tree_col)
new.samps$tree_col[new.samps$sample_name=="E7_8266_20090416"] <- "darkmagenta"
new.samps$tree_legend <- paste(new.samps$tree_ID, new.samps$drough.control, new.samps$developmental_stage, sep="-")

trees.order <- data.frame(legend=c("990-control-leaf_bud", "1099-control-leaf_bud", "1377-control-leaf_bud", "970-drought-leaf_bud", "8212-drought-leaf_bud", "8266-drought-leaf_bud", "8266-drought-flower_bud"), color=c("cyan", "green", "blue", "orange", "red", "magenta", "darkmagenta"), pch=c(16, 16, 16, 18, 18, 18, 18), cex=c(1, 1, 1, 1, 1, 1, 2)) 
trees.order$color <- as.character(trees.order$color)


# colors reprezening WP level 
colors.wp <- new.samps[, c("sample_name", "Water.Potential")]
colors.wp$Water.Potential <- 100*round(colors.wp$Water.Potential, digits=2)
colors.palette <- as.character(colorRampPalette(c("red","blue"))(diff(range(na.omit(colors.wp$Water.Potential)))))
colors.wp$colors <- colors.palette[colors.wp$Water.Potential - min(na.omit(colors.wp$Water.Potential))+1]
colors.wp$colors[is.na(colors.wp$colors)] <- "#838B8B" 
rownames(colors.wp) <- colors.wp$sample_name
colors.wp

new.samps$short.name <- paste(new.samps$tree_ID, paste(substr(new.samps$year.month.day, 7,8), substr(new.samps$year.month.day, 4,5), substr(new.samps$year.month.day, 1,2), sep="."), sep="_")
new.samps$short.name


### counts
x <- read.table("Data/raw_data-clean.csv", sep=",", header=T, row.names=1)
x <- x[, new.samps$sample_name]

flowered.samps <- c("E7_8266_20090416") # flowered sample
control.samps <- rownames(new.samps[new.samps$drough.control=="control",])
November.samps <- c("N9_970_20090114","M1_1099_20090114","I5_8212_20090119","F8_8266_20090119")

### list of all unique days in 2008 and 2009
days <- read.table("Data/Unique_days_short2.csv", sep=";")
all.days <- data.frame(dataF=days, numericF=as.numeric(strptime(days[,], "%d.%m.%y")))
library(stringr)
month.days <- all.days[str_sub(all.days[,1], 1, 2)=="01",]


###### files with control genes

AT.id <- read.table("Data/genes_descr_control/best_hit_blast_result_Sl_predicted_exons.txt", sep=",", stringsAsFactors=FALSE)

genes.description <- read.table("Data/genes_descr_control/genes_description.xls", sep="\t", header=T)

### Flowering genes - Table S4
Athaliana.flowering.genes <- read.table("Data/genes_descr_control/gene_list_flowering_related_MolEcol_Athaliana_S4.csv", sep=",", header=T, stringsAsFactors=FALSE, skip=1)
Athaliana.flowering.genes[,1] <- gsub(" ", "", Athaliana.flowering.genes[,1])

genes.full.description <- read.table("Data/genes_descr_control/genes_description_full.csv", header=T, sep=";", stringsAsFactors=F)

UP.genes <- read.table("Data/genes_descr_control/mDr_Day10_drought_up_regulated_genes_Harb_etal.csv", sep=";", header=T, stringsAsFactors=FALSE)
DOWN.genes <- read.table("Data/genes_descr_control/mDr_Day10_drought_down_regulated_genes_Harb_etal.csv", sep=";", header=T, stringsAsFactors=FALSE)
Drought.genes <- rbind(UP.genes[,1:3], DOWN.genes[,1:3])
names(Drought.genes) <- c("AT_ID", "Drought_regulation", "Description3")

genes.full.description <- merge(genes.full.description, Drought.genes, all.x=TRUE)
rownames(genes.full.description) <- genes.full.description$ID


# save.image("Shimizu_workspace.Rdata")


#####################################################################################################

### plots of expression for flowering genes

#####################################################################################################


#AT.genes <- Athaliana.flowering.genes.FT.SVP
AT.genes <- Athaliana.flowering.genes

#elim.samps=c(flowered.samps)
elim.samps=NULL
x <- x[,!names(x) %in% elim.samps]
new.samps <- new.samps[!new.samps$sample_name %in% elim.samps, ]

library(edgeR)
d.org <- DGEList(x, group=new.samps$tree_ID)
d.org <- calcNormFactors(d.org)

d.cpm <- cpm(d.org, normalized.lib.sizes=TRUE)
d.cpm.l <- log(d.cpm + min(d.cpm[d.cpm != 0]))


dir.create("Plots_of_flowering_genes/", showWarnings=F, recursive=T)


pdf(paste0("Plots_of_flowering_genes/" , "Flowering_genes_from_S4_table_Nov_fl" ,".pdf"), h=5, w=10)
info.table <- NULL

for(g in 1:nrow(AT.genes)){
  # g=1
  cat(paste(g, ", "))
  genes <- AT.id[AT.id[,2] == AT.genes[g, 1], 1]
  
  if(length(genes!=0)){
    for(j in 1:length(genes)){
      # j=1
      
      info.table <- rbind(info.table, data.frame(AT.ID=AT.genes[g, 1], AT.description=AT.genes[g, 2] , ID=genes[j]))
      
      plot(0, type="n", main=paste0(AT.genes[g, 1] ," - ", AT.genes[g, 2] , "\n",  genes[j]) ,xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(d.cpm[genes[j], ])), max(na.omit(d.cpm[genes[j], ]))), xlab="Time", ylab="Gene Expression in cpm", xaxt = "n")
      axis(side=1, at=month.days[,2], labels=month.days[,1])
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[246])
      for(t in trees.order$legend){
        # t=trees.order$legend[1]
        lines(new.samps$time_nr[new.samps$tree_legend == t], d.cpm[genes[j], new.samps$tree_legend == t] , col=trees.order$color[trees.order$legend==t], type="b", pch=trees.order$pch[trees.order$legend==t], cex=trees.order$cex[trees.order$legend==t]) 
      }
      legend("topright", legend = trees.order$legend, col=trees.order$color, cex=0.5, text.col=trees.order$color)
      
    }
    
  }
}

# write.table(info.table,"Plots_of_flowering_genes/info_table.xls", quote=FALSE, sep="\t", row.names=FALSE)

dev.off()



pdf(paste0("Plots_of_flowering_genes/" , "Variables" ,".pdf"), h=5, w=10)

plot.vars <- names(new.samps)[20:34]

for(v in plot.vars){
  plot(0, type="n", xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(new.samps[,v])), max(na.omit(new.samps[,v]))), xlab="Time", ylab=v, xaxt = "n")
  for(t in trees.order$legend){
    lines(new.samps$time_nr[new.samps$tree_legend==t], new.samps[new.samps$tree_legend==t, v], col=trees.order$color[trees.order$legend==t], type="b", pch=trees.order$pch[trees.order$legend==t], cex=trees.order$cex[trees.order$legend==t]) 
  }
  axis(side=1, at=month.days[,2], labels=month.days[,1])
  legend("topright", legend = trees.order$legend, col=trees.order$color, cex=0.5, text.col=trees.order$color)
  
}
dev.off()


### log cpm

pdf(paste0("Plots_of_flowering_genes/" , "Flowering_genes_from_S4_table_Nov_fl_log" ,".pdf"), h=5, w=10)

info.table <- NULL

for(g in 1:nrow(AT.genes)){
  # g=1
  cat(paste(g, ", "))
  genes <- AT.id[AT.id[,2] == AT.genes[g, 1], 1]
  
  if(length(genes!=0)){
    for(j in 1:length(genes)){
      # j=1
      
      info.table <- rbind(info.table, data.frame(AT.ID=AT.genes[g, 1], AT.description=AT.genes[g, 2] , ID=genes[j]))
      
      plot(0, type="n", main=paste0(AT.genes[g, 1] ," - ", AT.genes[g, 2] , "\n",  genes[j]) ,xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(d.cpm.l[genes[j], ])), max(na.omit(d.cpm.l[genes[j], ]))), xlab="Time", ylab="Gene Expression in cpm", xaxt = "n")
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[246])
      for(t in trees.order$legend){
        # t=trees.order$legend[1]
        lines(new.samps$time_nr[new.samps$tree_legend ==t], d.cpm.l[genes[j], new.samps$tree_legend ==t] , col=trees.order$color[trees.order$legend==t], type="b", pch=trees.order$pch[trees.order$legend==t], cex=trees.order$cex[trees.order$legend==t]) 
      }
      axis(side=1, at=month.days[,2], labels=month.days[,1])
      legend("topright", legend = trees.order$legend, col=trees.order$color, cex=0.5, text.col=trees.order$color)
      
    }
    
  }
}

plot.vars <- names(new.samps)[20:34]

for(v in plot.vars){
  plot(0, type="n", xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(new.samps[,v])), max(na.omit(new.samps[,v]))), xlab="Time", ylab=v, xaxt = "n")
  for(t in trees.order$legend){
    lines(new.samps$time_nr[new.samps$tree_legend==t], new.samps[new.samps$tree_legend==t, v], col=trees.order$color[trees.order$legend==t], type="b", pch=trees.order$pch[trees.order$legend==t], cex=trees.order$cex[trees.order$legend==t]) 
  }
  axis(side=1, at=month.days[,2], labels=month.days[,1])
  legend("topright", legend = trees.order$legend, col=trees.order$color, cex=0.5, text.col=trees.order$color)
  
}

dev.off()


####################################
### expression of SANGER genes
####################################

new.samps$sample_name_short <- substr(new.samps$sample_name, 1, nchar(new.samps$sample_name)-9)

sanger.lib.size <- read.table("Data/new_data_29Nov_SANGER/sanger_map_count_clean.csv",sep=",",nrow=1,stringsAsFactors=FALSE, row.names=1, header=T)
sanger.lib.size <- as.numeric(sanger.lib.size[new.samps$sample_name_short])
names(sanger.lib.size) <- new.samps$sample_name

sanger.x <- read.table("Data/new_data_29Nov_SANGER/sanger_map_count_clean.csv", sep=",", row.names=1, header=T)
sanger.x <- sanger.x[-1,]

sanger.x <- sanger.x[, new.samps$sample_name_short]
names(sanger.x) <- new.samps$sample_name

AT.genes <- Athaliana.flowering.genes

elim.samps=NULL
sanger.x <- sanger.x[,!names(sanger.x) %in% elim.samps]
new.samps <- new.samps[!new.samps$sample_name %in% elim.samps, ]

library(edgeR)
sanger.d.org <- DGEList(sanger.x, group=new.samps$tree_ID, lib.size=sanger.lib.size)
sanger.d.org <- calcNormFactors(sanger.d.org)

sanger.d.cpm <- cpm(sanger.d.org, normalized.lib.sizes=TRUE)
sanger.d.cpm.l <- log(sanger.d.cpm + min(sanger.d.cpm[sanger.d.cpm != 0]))

sanger.genes <- row.names(sanger.d.cpm)

pdf(paste0("Plots_of_flowering_genes/" , "Sanger_genes" ,".pdf"), h=5, w=10)

for(g in 1:nrow(sanger.d.cpm)){
# g=1
      plot(0, type="n", main=paste0(sanger.genes[g]) ,xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(sanger.d.cpm[g, ])), max(na.omit(sanger.d.cpm[g, ]))), xlab="Time", ylab="Gene Expression in cpm")
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[246])
      for(t in trees.order$legend){
        # t=trees.order$legend[1]
        lines(new.samps$time_nr[new.samps$tree_legend == t], sanger.d.cpm[g, new.samps$tree_legend == t] , col=trees.order$color[trees.order$legend==t], type="b", pch=trees.order$pch[trees.order$legend==t], cex=trees.order$cex[trees.order$legend==t]) 
      }
      legend("topleft", legend = trees.order$legend, col=trees.order$color, cex=0.5, text.col=trees.order$color)
      
}

dev.off()



#####################################################################################################

### GLM - models with one covariate 

#####################################################################################################

####################################
### run.edgeR()
####################################

# function

library(edgeR)

model.formula = ~ FT_gene
varialbs=c("FT_gene")
varialbs.plot = "FT_gene"


elim.samps = c(November.samps)
plots.path="Plots_RUN_edgeR"
plot.name=""
fit.clr = 1
LRTcoef=2
FDR=0.1


run.edgeR <- function(x, new.samps, model.formula, varialbs, elim.samps, FDR=0.1, plots.path="Plots_RUN_edgeR", plot.name="", varialbs.plot, fit.clr = 1, LRTcoef=2, trees.order=trees.order, genes.full.description=genes.full.description, AT.id=AT.id){ 
  
  dir.create(plots.path, recursive=T, showWarnings = FALSE)
  
  library(stringr)
  model.char <- str_replace_all(as.character(model.formula), " ", "") # can be gsub()
  model.char <- str_replace_all(model.char[2], "\\.", "_")
  
  
  x.tmp <- x[, new.samps$sample_name]  
  x.tmp <- x.tmp[,!names(x.tmp) %in% elim.samps]
  new.samps.tmp <- new.samps[!new.samps$sample_name %in% elim.samps, ]
  
#   for(i in varialbs){
#     x.tmp <- x.tmp[,!is.na(new.samps.tmp[, i])]
#     new.samps.tmp <- new.samps.tmp[!is.na(new.samps.tmp[, i]), ]
#   } 
  

  d <- DGEList(x.tmp, group=new.samps.tmp$tree_ID)
  d <- calcNormFactors(d)
  
  ### make sure a gene is expressed (CPM > 1) in more than 2 samples
  d.cpm <- cpm(d, normalized.lib.sizes=TRUE)
  
eps <- min(d.cpm[d.cpm!=0])

  new.samps.tmp$FT_gene <- log(d.cpm["GID031739_3527284",]+eps)
  new.samps.tmp$SVP_gene <- log(d.cpm["GID037469_3529277",]+eps)
  new.samps.tmp$FUL_gene <- log(d.cpm["GID056430_3537430",]+eps)

  d <- d[ rowSums(d.cpm>1) > 10, ]
  # cat(paste("*", nrow(d), "genes expressed in more than two samples \n"))
  
  
  # design model matrix
  design <- model.matrix(model.formula, data=new.samps.tmp)
  
  # estimate dispersion  
  d <- estimateGLMCommonDisp(d,design)
  d <- estimateGLMTrendedDisp(d,design)
  d <- estimateGLMTagwiseDisp(d,design)
  
  # glmFit
  fit <- glmFit(d,design)
  # glmLRT
  lrt <- glmLRT(fit, coef=LRTcoef)
  
  ## plot TOP genes // expression vs meteo var
  top.tags <- topTags(lrt, n=nrow(lrt$table))
  top.tags <- top.tags$table
  top.genes <- rownames(top.tags[top.tags$FDR < FDR,])
  length(top.genes)
  #fit$coefficients[top.genes,]

write.table( merge(fit$coefficients[top.genes,], top.tags[top.genes,], by=0) ,paste(plots.path, "/", model.char, "_top_fitting", plot.name,".txt", sep=""), sep="\t", row.names=F, quote=F)

  
  d.fit <- DGEList(fit$fitted.values, group=new.samps.tmp$tree_ID[match(names(x.tmp), new.samps.tmp$sample_name)])
  d.fit <- calcNormFactors(d.fit)
  
  d.cpm <- cpm(d, normalized.lib.sizes=TRUE)
  d.fit.cpm <- cpm(d.fit, normalized.lib.sizes=TRUE)
  

  Top.fitting.list <- list()
  UP.coeffs <- matrix(NA, length(top.genes), ncol(design))
  rownames(UP.coeffs) <- top.genes  
  DOWN.coeffs <- matrix(NA, length(top.genes), ncol(design))
  rownames(DOWN.coeffs) <- top.genes
  
  
  pdf( paste(plots.path, "/", model.char, "_top_fitting", plot.name,".pdf", sep=""), width = 10, height =7.5)
  par(mfrow=c(2,3))
  
  for(i in top.genes){
    # i="GID006843_3010436"  
    Top.fitting.list[[i]] <- data.frame(new.samps.tmp[, varialbs.plot], log(d.cpm[i,]), log(d.fit.cpm[i,]))
    colnames(Top.fitting.list[[i]]) <- c(varialbs.plot, "Log.Counts.in.CPM.RAW", "Log.Counts.in.CPM.FIT")
    rownames(Top.fitting.list[[i]]) <- row.names(fit$samples)
    coeffs <- round(fit$coefficients[i,], 2)
    
    AT <- AT.id[AT.id[,1] == i, 2]
    if.UP <- if.DOWN <- FALSE
    if(length(AT)!=0){
      if.UP <- AT %in% UP.genes[,1]
      if(if.UP)
        UP.coeffs[i,] <- coeffs
      if.DOWN <- AT %in% DOWN.genes[,1]
      if(if.DOWN)
        DOWN.coeffs[i,] <- coeffs
    }
    
    new.samps.tmp$tree_ID_nr <- new.samps.tmp$tree_ID
    levels(new.samps.tmp$tree_ID_nr) <- 1:length(levels(new.samps.tmp$tree_ID)) +1
    
    plot(Top.fitting.list[[i]][,1], Top.fitting.list[[i]][,2], main=paste(i, "/", AT, "\n", "Coefffs:",paste(coeffs, collapse=", "), "\n FDR:", top.tags[i, "FDR"], "\n UP:", if.UP, "DOWN:", if.DOWN), xlab = varialbs.plot, ylab="Log(Counts in CPM)", cex.main=0.7, pch=ifelse(new.samps.tmp$drough.control=="drought", 5, 1), col=3)
    points(Top.fitting.list[[i]][,1], Top.fitting.list[[i]][,3], col=1, pch="*", cex=2)
    
  }
  
  dev.off()
  
  
  pdf(paste(plots.path, "/", model.char, "_top_genes", plot.name,".pdf", sep=""), width = 10, height =5)
  
  for(j in 1:length(top.genes)){
    # j=1
    
    plot(0, type="n", main=paste0(top.genes[j]) ,xlim=c(min(new.samps.tmp$time_nr), max(new.samps.tmp$time_nr)), ylim=c(min(na.omit(d.cpm[top.genes[j], ])), max(na.omit(d.cpm[top.genes[j], ]))), xlab="Time", ylab="Gene Expression in cpm", xaxt = "n")
    axis(side=1, at=month.days[,2], labels=month.days[,1])
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[246])
    for(t in trees.order$legend){
      # t=trees.order$legend[1]
      lines(new.samps.tmp$time_nr[new.samps.tmp$tree_legend == t], d.cpm[top.genes[j], new.samps.tmp$tree_legend == t] , col=trees.order$color[trees.order$legend==t], type="b", pch=trees.order$pch[trees.order$legend==t], cex=trees.order$cex[trees.order$legend==t]) 
    }
    legend("topright", legend = trees.order$legend, col=trees.order$color, cex=0.5, text.col=trees.order$color)
    
  }

  dev.off()

  
}


# runs

library(edgeR)

model.formula = ~ FUL_gene
varialbs=c("FUL_gene")
varialbs.plot = "FUL_gene"

model.formula = ~ SVP_gene
varialbs=c("SVP_gene")
varialbs.plot = "SVP_gene"

model.formula = ~ FT_gene
varialbs=c("FT_gene")
varialbs.plot = "FT_gene"

elim.samps = c(November.samps)
plots.path="Plots_RUN_edgeR/"
plot.name=""
fit.clr = 1
LRTcoef=2
FDR=0.01


run.edgeR(x, new.samps, model.formula=model.formula, varialbs=varialbs, elim.samps=elim.samps, FDR=FDR, plots.path=plots.path, plot.name=plot.name, varialbs.plot=varialbs.plot, fit.clr = fit.clr, LRTcoef=LRTcoef, trees.order=trees.order, genes.full.description=genes.full.description, AT.id=AT.id)



model.formula = ~ Water.Potential
varialbs=c("Water.Potential")
varialbs.plot = "Water.Potential"

model.formula = ~ Soil.Moisture
varialbs=c("Soil.Moisture")
varialbs.plot = "Soil.Moisture"


elim.samps = c(flowered.samps, control.samps, November.samps, new.samps[new.samps$tree_ID=="8212", "sample_name"], new.samps[new.samps$tree_ID=="990", "sample_name"])
plots.path="Plots_RUN_edgeR"
plot.name=""
fit.clr = 1
LRTcoef=2
FDR=0.01


run.edgeR(x, new.samps, model.formula=model.formula, varialbs=varialbs, elim.samps=elim.samps, FDR=FDR, plots.path=plots.path, plot.name=plot.name, varialbs.plot=varialbs.plot, fit.clr = fit.clr, LRTcoef=LRTcoef, trees.order=trees.order, genes.full.description=genes.full.description, AT.id=AT.id)


####################################
### run.edgeR.robust()
####################################

source("/home/gosia/R/R_Shimizu_RNA_seq/Run_edgeR_robust.R")
library(edgeR)



model.formula = ~ FUL_gene
varialbs=c("FUL_gene")
varialbs.plot = "FUL_gene"

model.formula = ~ SVP_gene
varialbs=c("SVP_gene")
varialbs.plot = "SVP_gene"

model.formula = ~ FT_gene
varialbs=c("FT_gene")
varialbs.plot = "FT_gene"

elim.samps = c(November.samps)
plots.path="Plots_RUN_edgeR_robust/"
plot.name=""
fit.clr = 1
LRTcoef=2
FDR=0.01



run.edgeR.robust(x, new.samps, model.formula, varialbs, elim.samps, FDR=FDR, plots.path=plots.path, plot.name=plot.name, varialbs.plot, fit.clr = fit.clr, LRTcoef=LRTcoef, trees.order=trees.order, genes.full.description=genes.full.description, AT.id=AT.id)



model.formula = ~ Water.Potential
varialbs=c("Water.Potential")
varialbs.plot = "Water.Potential"

model.formula = ~ Soil.Moisture
varialbs=c("Soil.Moisture")
varialbs.plot = "Soil.Moisture"


elim.samps = c(flowered.samps, control.samps, November.samps, new.samps[new.samps$tree_ID=="8212", "sample_name"], new.samps[new.samps$tree_ID=="990", "sample_name"])
plots.path="Plots_RUN_edgeR_robust"
plot.name=""
fit.clr = 1
LRTcoef=2
FDR=0.8

run.edgeR.robust(x, new.samps, model.formula, varialbs, elim.samps, FDR=FDR, plots.path=plots.path, plot.name=plot.name, varialbs.plot, fit.clr = fit.clr, LRTcoef=LRTcoef, trees.order=trees.order, genes.full.description=genes.full.description, AT.id=AT.id)


#####################################################################################################

### clustering

#####################################################################################################

setwd("/home/Shared/data/seq/Shimizu_RNA_seq/")
load("Shimizu_workspace.Rdata")


############################################
### select genes for clustering
############################################

x.orig <- x
new.samps.orig <- new.samps

# elim.samps=c(flowered.samps, "K4_990_20081205", "K5_990_20090511")
elim.samps=NULL
x <- x.orig[,!names(x.orig) %in% elim.samps]
new.samps <- new.samps.orig[!new.samps.orig$sample_name %in% elim.samps, ]


library(limma)
library(edgeR)

d.org <- DGEList(x, group=new.samps$tree_ID)
d.org <- calcNormFactors(d.org)

### make sure a gene is expressed (CPM > 1) in more than 2 samples
d.cpm.org <- cpm(d.org, normalized.lib.sizes=TRUE)
dim(d.cpm.org)

# sum(d.cpm.org["GID031739_3527284",] > 1)
# sum( rowSums(d.cpm.org > 1) > 2 )

d.org <- d.org[ rowSums(d.cpm.org > 1) > 2, ]
d.cpm.org <- cpm(d.org, normalized.lib.sizes=TRUE)
d.cpm.org.l <- log(d.cpm.org  + min(d.cpm.org [d.cpm.org  != 0]))

count.data <- list()
count.data$d.cpm.org <- d.cpm.org
count.data$d.cpm.org.l <- d.cpm.org.l
selected.genes <- list()
selected.genes$all.genes <- rownames(d.cpm.org)
selected.genes$AT.genes <- selected.genes$all.genes[selected.genes$all.genes %in% AT.id[,1]]


############################################
### select genes for clustering 2 - with high variance
############################################


# save.image("Clustering_workspace_all_genes.Rdata")

### genes with high variance
dir.create("Clustering/Variance/", showWarnings=F, recursive=T)

pdf("Clustering/Variance/CPM_Mean-variance_trend.pdf")
mean <- apply(d.cpm.org.l, 1, mean)
sd <- apply(d.cpm.org.l, 1, sd)
plot(mean,sd, pch=16, main="Log(cpm)")

mean <- apply(d.cpm.org, 1, mean)
sd <- apply(d.cpm.org, 1, sd)
plot(mean,sd, pch=16, main="cpm", log="xy")

plotMeanVar(d.org, show.raw.vars=T)

dev.off()


### assings weights based on mean-variance trend - voom

pdf("Clustering/Variance/Voom_estimation_Mean-variance_trend.pdf")
d.voom <- voom(counts=d.org, plot = TRUE, span=0.6)
dev.off()

d.voom.weight <- (d.voom$E * d.voom$weights) # values in log2

# row.weights <- rowSums(d.voom$weights)

# pdf("Clustering/Voom_weighted_Mean-variance_w_trend.pdf")
# d.voom.mean <- unlist(lapply(1:nrow(d.voom.weight), function(v){ 
#   sum(d.voom.weight[v, ])/row.weights[v] 
# }))
# d.voom.sd <- unlist(lapply(1:nrow(d.voom.weight), function(v){ 
#   sqrt(sum(d.voom$weights[v,]*(d.voom$E[v, ]-d.voom.mean[v])^2)/row.weights[v])
# }))
# plot(d.voom.mean, d.voom.sd)
# dev.off()

pdf("Clustering/Variance/Voom_Mean-variance_trend.pdf")
d.voom.mean <- apply(d.voom.weight, 1, mean)
d.voom.sd <- apply(d.voom.weight, 1, sd)
#l <- lowess(d.voom.mean, d.voom.sd)
plot(d.voom.mean, d.voom.sd, pch=16, main="Voom weighted values")
#lines(l, col=2)
dev.off()


names(d.voom.sd) <- rownames(d.voom$E)
selected.genes$var.genes.voom <- names(sort(d.voom.sd, decreasing=TRUE)[1:5000])
count.data$d.voom.weight <- d.voom.weight


pdf("Clustering/Variance/Voom_sel_CPM_Mean-variance_trend.pdf")
mean <- apply(d.cpm.org.l, 1, mean)
sd <- apply(d.cpm.org.l, 1, sd)
plot(mean,sd, pch=16, main="Log(cpm)")
points(mean[selected.genes$var.genes.voom], sd[selected.genes$var.genes.voom], pch=16, col=2)
points(mean[genes.full.description[!is.na(genes.full.description$Flowering), "ID"]], sd[genes.full.description[!is.na(genes.full.description$Flowering), "ID"]], pch=16, col=3)

mean <- apply(d.cpm.org, 1, mean)
sd <- apply(d.cpm.org, 1, sd)
plot(mean,sd, pch=16, main="cpm", log="xy")
points(mean[selected.genes$var.genes.voom], sd[selected.genes$var.genes.voom], pch=16, col=2)
points(mean[genes.full.description[!is.na(genes.full.description$Flowering), "ID"]], sd[genes.full.description[!is.na(genes.full.description$Flowering), "ID"]], pch=16, col=3)

dev.off()


### variance stabilization - DESeq 

library(DESeq)
cds <- newCountDataSet(countData=d.org$counts, conditions=new.samps$sample_ID)
cds <- estimateSizeFactors( cds )
cds <- estimateDispersions( cds, method="blind" )
vsd <- getVarianceStabilizedData( cds )

pdf("Clustering/Variance/DESeq_Mean-variance_trend.pdf")
plot(rowMeans(vsd), sqrt(genefilter::rowVars(vsd)), main="DESeq VSD")
dev.off()


selected.genes$var.genes.vs <- names(sort(genefilter::rowVars(vsd), decreasing=TRUE)[1:5000])
count.data$vsd <- vsd

# vsd["GID063259_3539580",]

pdf("Clustering/Variance/DESeq_sel_CPM_Mean-variance_trend.pdf")

mean <- apply(d.cpm.org.l, 1, mean)
sd <- apply(d.cpm.org.l, 1, sd)
plot(mean,sd, pch=16, main="Log(cpm)")
points(mean[selected.genes$var.genes.vs], sd[selected.genes$var.genes.vs], pch=16, col=2)
points(mean[genes.full.description[!is.na(genes.full.description$Flowering), "ID"]], sd[genes.full.description[!is.na(genes.full.description$Flowering), "ID"]], pch=16, col=3)

mean <- apply(d.cpm.org, 1, mean)
sd <- apply(d.cpm.org, 1, sd)
plot(mean,sd, pch=16, main="cpm", log="xy")
points(mean[selected.genes$var.genes.vs], sd[selected.genes$var.genes.vs], pch=16, col=2)
points(mean[genes.full.description[!is.na(genes.full.description$Flowering), "ID"]], sd[genes.full.description[!is.na(genes.full.description$Flowering), "ID"]], pch=16, col=3)

dev.off()


"GID053198_3536068" %in% selected.genes$var.genes.vs



############################################
# load clustering workspace on server
############################################

dir.create("Clustering", showWarnings=F, recursive=T)

# install.packages("/home/gosia/R/packages/", "/home/gosia/R/libraries/3.0.2/")
library(cluster)
library(clusterSim)
library(lattice)
library(parallel)
library(gtools)
library(gplots) 
library(RColorBrewer)
source("/home/gosia/R/R_Shimizu_RNA_seq/heatmap2.r")


############################################
# PAM run
############################################

source("/home/gosia/R/R_Shimizu_RNA_seq/PAM_clustering.R")

pam.x <- count.data$d.cpm.org
# pam.x <- pam.x[1:200,]
dim(pam.x)


PAM.clustering(pam.x=pam.x, new.samps=new.samps, genes.full.description, dist.method= "euclidian", norm.method="norm", out.path="Clustering/PAM/", out.name="CPM_norm", prior.nr.cl=2:50, mc.cores=1, ylim=c(-2,6), clustering.samps=NA)



PAM.clustering(pam.x=pam.x, new.samps=new.samps, genes.full.description, dist.method= "euclidian", norm.method="norm", out.path="Clustering/PAM_per_tree/", out.name="8266", prior.nr.cl=2:20, mc.cores=1, ylim=c(-2,6), clustering.samps=new.samps$sample_name[new.samps$tree_ID=="8266"])



############################################
# Kmeans run
############################################

source("/home/gosia/R/R_Shimizu_RNA_seq/Kmeans_clustering.R")
# source("/home/gosia/R/R_Shimizu_RNA_seq/Kmeans_clustering_tmp.R")

pam.x <- count.data$d.cpm.org
#pam.x <- pam.x[1:200,]
dim(pam.x)



# all
kmeans.clustering(pam.x, new.samps, genes.full.description, colors.wp, norm.method="norm", out.path="Clustering/Kmeans_all/", out.name="CPM_norm", prior.nr.cl=10:40, mc.cores=5, ylim=c(-2, 5), clustering.samps="all", colors=c("yellow","blue"))


# all, no November
samps <- new.samps[!new.samps$sample_name %in% November.samps, "sample_name"]

kmeans.clustering(pam.x[,samps], new.samps[samps,], genes.full.description, colors.wp, norm.method="norm", out.path="Clustering/Kmeans_all_noNov/", out.name="CPM_norm_noNov", prior.nr.cl=10:40, mc.cores=5, ylim=c(-2, 5), clustering.samps="all", colors=c("yellow","blue"))



# drought trees (no flowering, no outlier F1_8266)
drought.samp <- new.samps[new.samps$tree_ID %in% c("970", "8266"), "sample_name"]
# no flowering sample, no outlier F1_8266
drought.samp <- drought.samp[!drought.samp %in% c("E7_8266_20090416", "F1_8266_20081202")]

kmeans.clustering(pam.x[,drought.samp], new.samps[drought.samp,], genes.full.description, colors.wp, norm.method="norm", out.path="Clustering/Kmeans_drought_trees/", out.name="CPM_norm_drought", prior.nr.cl=5:40, mc.cores=5, ylim=c(-2, 5), clustering.samps="all", colors=c("yellow","blue"))



# drought trees (no flowering, no November)
drought.samp <- new.samps[new.samps$tree_ID %in% c("970", "8266"), "sample_name"]
drought.samp <- drought.samp[!drought.samp %in% c("E7_8266_20090416")]
drought.samp <- drought.samp[!drought.samp %in% November.samps]

kmeans.clustering(pam.x[,drought.samp], new.samps[drought.samp,], genes.full.description, colors.wp, norm.method="norm", out.path="Clustering/Kmeans_drought_trees_noNov/", out.name="CPM_norm_drought_noNov", prior.nr.cl=5:40, mc.cores=5, ylim=c(-2, 5), clustering.samps="all", colors=c("yellow","blue"))



# per tree

kmeans.clustering(pam.x, new.samps, genes.full.description, colors.wp, norm.method="norm", out.path="Clustering/Kmeans_per_tree/", out.name="8266", prior.nr.cl=10:20, mc.cores=5, ylim=c(-2, 5), clustering.samps=new.samps$sample_name[new.samps$tree_ID=="8266"], colors=c("yellow","blue"))


kmeans.clustering(pam.x, new.samps, genes.full.description, colors.wp, norm.method="norm", out.path="Clustering/Kmeans_per_tree/", out.name="970", prior.nr.cl=5:20, mc.cores=5, ylim=c(-2, 5), clustering.samps=new.samps$sample_name[new.samps$tree_ID=="970"], colors=c("yellow","blue"))


kmeans.clustering(pam.x, new.samps, genes.full.description, colors.wp, norm.method="norm", out.path="Clustering/Kmeans_per_tree/", out.name="1099", prior.nr.cl=5:20, mc.cores=5, ylim=c(-2, 5), clustering.samps=new.samps$sample_name[new.samps$tree_ID=="1099"], colors=c("yellow","blue"))


kmeans.clustering(pam.x, new.samps, genes.full.description, colors.wp, norm.method="norm", out.path="Clustering/Kmeans_per_tree/", out.name="1377", prior.nr.cl=5:20, mc.cores=5, ylim=c(-2, 5), clustering.samps=new.samps$sample_name[new.samps$tree_ID=="1377"], colors=c("yellow","blue"))





#####################################################################################################

### GO analysis with topGO

#####################################################################################################


setwd("/home/Shared/data/seq/Shimizu_RNA_seq/")
load("Shimizu_workspace.Rdata")

dir.create("GO", showWarnings=F, recursive=T)

library(topGO)
library(org.At.tair.db)

# go.name <- "Kmeans_all"
# whole.clustering <- read.table("Clustering/Kmeans_all/CPM_norm_kmeans_table.xls", sep="\t", header=T)
# nr.clusters <- 24

go.name <- "Kmeans_drought_trees"
whole.clustering <- read.table("Clustering/Kmeans_drought_trees/CPM_norm_drought_kmeans_table.xls", sep="\t", header=T)
nr.clusters <- 20

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

# cluster <- 3
sel.cluster.genes <- whole.clustering$AT_ID[whole.clustering[, paste0("clustering", nr.clusters)]==cluster]

gene.vector=as.numeric(assayed.genes %in% sel.cluster.genes)
names(gene.vector) <- assayed.genes
names(assayed.genes) <- gene.vector


  for(go in c("BP","MF","CC")){
    # go = "CC"
    
    sampleGOdata <- new("topGOdata", description = "Simple session", ontology = go, allGenes = gene.vector, geneSel = fun.gene.sel , nodeSize = 10, annot = annFUN.org, mapping = "org.At.tair.db")
# Fisher's exact test which is based on gene counts, and a Kolmogorov-Smirnov like test which computes enrichment based on gene scores
# For the method classic each GO category is tested independently
    resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
    resultFisher.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher")
    

pValues.Fisher <- score(resultFisher)
topNodes.Fisher <- sum(pValues.Fisher < 0.05)
topNodes.Fisher
pValues.Fisher.elim <- score(resultFisher.elim)
topNodes.Fisher.elim <- sum(pValues.Fisher.elim < 0.0001)
topNodes.Fisher.elim

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

if(topNodes.Fisher > 0){
  
  allRes.Fisher <- GenTable(sampleGOdata, classicFisher = resultFisher, elimFisher = resultFisher.elim, orderBy = "classicFisher", ranksOf = "elimFisher", topNodes = topNodes.Fisher+1)
  
  allRes.Fisher$GO <- go
  allRes.Fisher$cluster <- cluster
  
  allRes.Fisher.merged <- rbind(allRes.Fisher.merged, allRes.Fisher[-c(topNodes.Fisher+1), , drop=F])
}


if(topNodes.Fisher > 0){
  
  allRes.Fisher.elim <- GenTable(sampleGOdata, classicFisher = resultFisher, elimFisher = resultFisher.elim, orderBy = "elimFisher", ranksOf = "classicFisher", topNodes = topNodes.Fisher.elim+1 )
  
  allRes.Fisher.elim$GO <- go
  allRes.Fisher.elim$cluster <- cluster
  
  allRes.Fisher.elim.merged <- rbind(allRes.Fisher.elim.merged, allRes.Fisher.elim[-c(topNodes.Fisher+1), , drop=F])
}


#     pdf(paste("GO/GO_", go, ".pdf", sep=""), width = 10, height = 10)
#     showSigOfNodes(sampleGOdata, score(resultKS), firstSigNodes = 5, useInfo = 'all')
#     dev.off()
    
    
  }
  

}

write.table(allRes.Fisher.merged, paste("GO/",go.name,"_GO_Fisher" , ".xls", sep=""), sep="\t", row.names=F)

write.table(allRes.Fisher.elim.merged, paste("GO/",go.name,"_GO_Fisher_elim" , ".xls", sep=""), sep="\t", row.names=F)



allRes.Fisher.elim.merged.s <- allRes.Fisher.elim.merged[allRes.Fisher.elim.merged$elimFisher < 0.0000001, ]

write.table(allRes.Fisher.elim.merged.s, paste("GO/",go.name,"_GO_Fisher_elim_1E-07" , ".xls", sep=""), sep="\t", row.names=F)



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


# go.name <- "Kmeans_all"
# whole.clustering <- read.table("Clustering/Kmeans_all/CPM_norm_kmeans_table.xls", sep="\t", header=T)
# nr.clusters <- 24

go.name <- "Kmeans_drought_trees"
whole.clustering.org <- read.table("Clustering/Kmeans_drought_trees/CPM_norm_drought_kmeans_table.xls", sep="\t", header=T)
nr.clusters <- 20


whole.clustering <- whole.clustering.org[!is.na(whole.clustering.org$AT_ID), c("ID", "AT_ID", paste0("clustering", nr.clusters))]
colnames(whole.clustering) <- c("ID", "AT_ID","clustering")


gene.set <- Athaliana.flowering.genes[,1]
name.out <- "Athaliana_flowering_genes"


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
  
  if(any(fisher.x < 10))
    fisher.table[cluster, "pvalues"] <- NA
  
  fisher.table[cluster, "X1"] <- fisher.x[1,1]
  fisher.table[cluster, "X2"] <- fisher.x[1,2]
  fisher.table[cluster, "X3"] <- fisher.x[2,1]
  fisher.table[cluster, "X4"] <- fisher.x[2,2]
  fisher.table[cluster,  "cluster.size"] <- fisher.x[1,1] + fisher.x[2,1]
 
  
}

fisher.table <- data.frame(cluster=1:nr.clusters, fisher.table)

write.table(fisher.table, paste0(dir.out,"/", name.out, "_", go.name, "_Fisher_pvalues.xls"), sep="\t", row.names=F)


##############################
### gene sets
##############################

setwd("/home/Shared/data/seq/Shimizu_RNA_seq/")
load("Shimizu_workspace.Rdata")

dir.out <- "Fisher_test/Gene_sets"

dir.create(dir.out, showWarnings=F, recursive=T)


# go.name <- "Kmeans_all"
# whole.clustering <- read.table("Clustering/Kmeans_all/CPM_norm_kmeans_table.xls", sep="\t", header=T)
# nr.clusters <- 24

go.name <- "Kmeans_drought_trees"
whole.clustering.org <- read.table("Clustering/Kmeans_drought_trees/CPM_norm_drought_kmeans_table.xls", sep="\t", header=T)
nr.clusters <- 20


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


go.name <- "Kmeans_all"
whole.clustering <- read.table("Clustering/Kmeans_all/CPM_norm_kmeans_table.xls", sep="\t", header=T)
nr.clusters <- 24

# go.name <- "Kmeans_drought_trees"
# whole.clustering.org <- read.table("Clustering/Kmeans_drought_trees/CPM_norm_drought_kmeans_table.xls", sep="\t", header=T)
# nr.clusters <- 20


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






