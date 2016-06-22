#####################################################################################################

### GLM // gene expression vs Soil Moisture & Water Potential & Temperature
### MODEL SELECTION; fit all combinations of models; 
### GO analysis

#####################################################################################################

setwd("/Users/gosia/DATA/Shimizu_RNA_seq/")


####################################################
### Correlation between factors
####################################################

setwd("/Users/gosia/DATA/Shimizu_RNA_seq/")
load(file="Data/new_samps_interpolation.RData")

names(new.samps)

# library(Hmisc)
# rcorr(new.samps$Soil.Moisture, new.samps$Water.Potential, type="pearson")

cor(new.samps$Soil.Moisture, new.samps$Water.Potential, method="pearson", use="complete.obs")
cor(new.samps$Soil.Moisture, new.samps$Water.Potential, method="spearman", use="complete.obs")

plot(new.samps$Soil.Moisture, new.samps$Water.Potential)


cor(new.samps$Temperature[new.samps$drough.control=="control"], new.samps$Water.Potential[new.samps$drough.control=="control"], method="pearson", use="complete.obs")
plot(new.samps$Temperature[new.samps$drough.control=="control"], new.samps$Water.Potential[new.samps$drough.control=="control"])

cor(new.samps$Temperature[new.samps$drough.control=="control"], new.samps$Soil.Moisture[new.samps$drough.control=="control"], method="pearson", use="complete.obs")
plot(new.samps$Temperature[new.samps$drough.control=="control"], new.samps$Soil.Moisture[new.samps$drough.control=="control"])


library(graphics)

new.samps$tree_ID_nr <- new.samps$tree_ID
levels(new.samps$tree_ID_nr) <- 1:length(levels(new.samps$tree_ID)) +1

pdf("Plots_Raw/Correlation.pdf", width = 10, height = 10)
pairs(new.samps[,c("Soil.Moisture", "Water.Potential", "Temperature")], col=new.samps$tree_ID_nr, pch=19)
dev.off()

pdf("Plots_Raw/Correlation28.pdf", width = 10, height = 10)
pairs(new.samps[,c("Soil.Moisture28", "Water.Potential28", "Temperature28")], col=new.samps$tree_ID_nr, pch=19)
dev.off()

pdf("Plots_Raw/Correlation14.pdf", width = 10, height = 10)
pairs(new.samps[,c("Soil.Moisture14", "Water.Potential14", "Temperature14")], col=new.samps$tree_ID_nr, pch=19)
dev.off()


### Partial correlation coefficient
### From formulas in Sheskin, 3e
### a,b=variables to be correlated, c=variable to be partialled out of both
pcor = function(a,b,c)
{
  (cor(a,b)-cor(a,c)*cor(b,c))/sqrt((1-cor(a,c)^2)*(1-cor(b,c)^2))
}
### end of script

#####################################################################################################

### load data

#####################################################################################################

setwd("/Users/gosia/DATA/Shimizu_RNA_seq/")

### samples & factors
load(file="Data/new_samps_interpolation.RData")
head(new.samps)
rownames(new.samps) <- new.samps$sample_name
### counts
x <- read.table("Data/raw_data-clean.csv", sep=",", header=T, row.names=1)
x <- x[, new.samps$sample_name]


flowered.samps <- c("E7_8266_20090416") # flowered sample
control.samps <- rownames(new.samps[new.samps$drough.control=="control",])



### Control genes

AT.id <- read.table("Data/genes_descr_control/best_hit_blast_result_Sl_predicted_exons.txt", sep=",", stringsAsFactors=FALSE)
head(AT.id)

UP.genes <- read.table("Data/genes_descr_control/mDr_Day10_drought_up_regulated_genes_Harb_etal.csv", sep=";", header=T, stringsAsFactors=FALSE)
DOWN.genes <- read.table("Data/genes_descr_control/mDr_Day10_drought_down_regulated_genes_Harb_etal.csv", sep=";", header=T, stringsAsFactors=FALSE)

Drought.genes <- read.table("Data/genes_descr_control/mDr_Day10_drought_regulated_genes_Harb_etal.csv", sep=";", header=T, stringsAsFactors=FALSE)

# MolEcol.genes.c1 <- read.table("Data/genes_descr_control/gene_list_in_clusterI_Mol_Ecol.csv", sep=";", header=T, stringsAsFactors=FALSE)
# MolEcol.genes.c2 <- read.table("Data/genes_descr_control/gene_list_in_clusterII_Mol_Ecol.csv", sep=";", header=T, stringsAsFactors=FALSE)

### DE unigenes for S.beccariana overlap with BLAST A.thaliana
MolEcol.DE.genes <- read.table("Data/genes_descr_control/gene_list_in_all_clusters_DEgenes_Mol_Ecol.csv", sep=";", header=T, stringsAsFactors=FALSE)


# Athaliana.flowering.genes <- read.table("Data/genes_descr_control/gene_list_flowering_related_MolEcol_Athaliana_S4.csv", sep=";", header=T, stringsAsFactors=FALSE)
# head(Athaliana.flowering.genes )

####################################################
###  clustering
####################################################

library(edgeR)
d <- DGEList(x, group=new.samps$tree_ID)
d <- calcNormFactors(d)

# make sure a gene is expressed (CPM > 1) in more than 2 samples
cps <- cpm(d, normalized.lib.sizes=TRUE)
d <- d[ rowSums(cps>1) > 2, ]

sample.dist <- dist(t(d$counts))
sample.clust <- hclust(sample.dist, method = "complete")
plot(sample.clust)

heatmap(d$counts, Rowv = NA)

heatmap(cor(d$counts, method = "pearson"), symm = TRUE)

library(gplots)
heatmap(diffExpMatrix,Colv=NA, col = redgreen(75),labRow=NA)


####################################################
###  Fitting model automaticly
####################################################

### run.edgeR()

model.formula = ~ Soil.Moisture + Water.Potential + Temperature
varialbs=c("Soil.Moisture", "Water.Potential", "Temperature")
elim.samps=flowered.samps
n.top=200
plots.path="PlotsRUNedgeR"
plot.name=""
control.genes=MolEcol.DE.genes
varialbs.plot = "Soil.Moisture"
fit.clr = 2
LRTcoef=2


run.edgeR <- function(x, new.samps, model.formula, varialbs, elim.samps, n.top=200, plots.path="PlotsTest", plot.name="", AT.id, UP.genes, DOWN.genes, varialbs.plot, fit.clr = 3, LRTcoef=2){ 
  
  x <- x[, new.samps$sample_name] + 1  
  for(i in varialbs){
    x <- x[,!is.na(new.samps[, i])]
    new.samps <- new.samps[!is.na(new.samps[, i]), ]
  } 
  x <- x[,!names(x) %in% elim.samps]
  new.samps <- new.samps[!new.samps$sample_name %in% elim.samps, ]
  
  library(edgeR)
  d <- DGEList(x, group=new.samps$tree_ID)
  d <- calcNormFactors(d)
  
  ### make sure a gene is expressed (CPM > 1) in more than 2 samples
  d.cpm <- cpm(d, normalized.lib.sizes=TRUE)
  d <- d[ rowSums(d.cpm>1) > 2, ]
  # cat(paste("*", nrow(d), "genes expressed in more than two samples \n"))
  
  
  # design model matrix
  design <- model.matrix(model.formula, data=new.samps)
  
  # estimate dispersion  
  d <- estimateGLMCommonDisp(d,design)
  d <- estimateGLMTrendedDisp(d,design)
  d <- estimateGLMTagwiseDisp(d,design)
  
  # glmFit
  fit <- glmFit(d,design)

  # make contrast to test any difference from control
#   mc <- makeContrasts(contrasts=c("Intercept", "Soil.Moisture"),levels=c("Intercept", "Soil.Moisture"))
#   colnames(mc) <- c("Intercept", "Soil.Moisture")
#   lrt <- glmLRT(fit, contrast=mc)
  
  lrt <- glmLRT(fit, coef=LRTcoef)

  ## plot TOP genes // expression vs soil moisture
  top.tags <- topTags(lrt, n=n.top)
  top.genes <- rownames(top.tags$table)
  
  d.fit <- DGEList(fit$fitted.values, group=new.samps$tree_ID[match(names(x), new.samps$sample_name)])
  d.fit <- calcNormFactors(d.fit)
  
  d.cpm <- cpm(d, normalized.lib.sizes=TRUE)
  d.fit.cpm <- cpm(d.fit, normalized.lib.sizes=TRUE)
  
  library(stringr)
  model.char <- str_replace_all(as.character(model.formula), " ", "") # can be gsub()
  model.char <- str_replace_all(model.char[2], "\\.", "_")
  
  #output.name <- str_replace_all(paste("Top_fitting_", model.char, "_cpm", plot.name, sep=""),"[[:punct:]]", "." )
  #assign("Top.fitting.list", list())
  Top.fitting.list <- list()
  UP.coeffs <- matrix(NA, length(top.genes), ncol(design))
  rownames(UP.coeffs) <- top.genes
  
  DOWN.coeffs <- matrix(NA, length(top.genes), ncol(design))
  rownames(DOWN.coeffs) <- top.genes
  
  dir.create(plots.path, recursive=T, showWarnings = FALSE)
  
  pdf( paste(plots.path, "/Top_fitting_", model.char, "_cpm", plot.name,".pdf", sep=""), width = 9, height =6)
  par(mfrow=c(2,3))
  
  for(i in top.genes){
    # i="GID006843_3010436"  
    Top.fitting.list[[i]] <- data.frame(new.samps[, varialbs.plot], log(d.cpm[i,]), log(d.fit.cpm[i,]))
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
    
    new.samps$tree_ID_nr <- new.samps$tree_ID
    levels(new.samps$tree_ID_nr) <- 1:length(levels(new.samps$tree_ID)) +1
    
    plot(Top.fitting.list[[i]][,1], Top.fitting.list[[i]][,2], main=paste(i, "/", AT, "\n", "Coefffs:",paste(coeffs, collapse=", "), "\n FDR:", top.tags$table[i, "FDR"], "\n UP:", if.UP, "DOWN:", if.DOWN), xlab = varialbs.plot, ylab="Log(Counts in CPM)", cex.main=0.7, pch=ifelse(new.samps$drough.control=="drought", 5, 1), col=new.samps$tree_ID_nr)
    points(Top.fitting.list[[i]][,1], Top.fitting.list[[i]][,3], col=1, pch="*", cex=2)
    
  }
  
  dev.off()
  
  #save(Top.fitting.list, file=paste(plots.path, "/Top_fitting_", model.char, "_cpm", plot.name,".RData", sep=""))
  
  UP.coeffs <- UP.coeffs[complete.cases(UP.coeffs), ]  
  DOWN.coeffs <- DOWN.coeffs[complete.cases(DOWN.coeffs), ]
  
#   invisible(list(Top.genes=top.tags$table, Top.fitting.list=Top.fitting.list, UP.coeffs=UP.coeffs, DOWN.coeffs=DOWN.coeffs))
  
  invisible(list(Coeff=lrt$coefficients, Dev=lrt$deviance, Des=lrt$design, Tab=lrt$table , UP.coeffs=UP.coeffs, DOWN.coeffs=DOWN.coeffs))
  
}



run.edgeR.models.fit <- function(x, new.samps, varialbs, elim.samps, out.path="Models", out.name="Models_fitting"){ 
  
  x <- x[, new.samps$sample_name] + 1  
  for(i in varialbs){
    x <- x[,!is.na(new.samps[, i])]
    new.samps <- new.samps[!is.na(new.samps[, i]), ]
  } 
  x <- x[,!names(x) %in% elim.samps]
  new.samps <- new.samps[!new.samps$sample_name %in% elim.samps, ]
  
  library(edgeR)
  d <- DGEList(x, group=new.samps$tree_ID)
  d <- calcNormFactors(d)
  
  ### make sure a gene is expressed (CPM > 1) in more than 2 samples
  d.cpm <- cpm(d, normalized.lib.sizes=TRUE)
  d <- d[ rowSums(d.cpm>1) > 2, ]
  
  models.results <- list()
  
  for(c in length(varialbs):1){
    # c=2
    models.comb <- combn(varialbs, c)
    
    for(m in 1:ncol(models.comb)){
      # m=1
      
      model.formula <- as.formula(paste0("~", paste(models.comb[,m], collapse="+")))
      cat("Fitting: ", paste(models.comb[,m], collapse="+"), "\n")
      design <- model.matrix(model.formula, data=new.samps)
      
      # estimate dispersion  
      d <- estimateGLMCommonDisp(d,design)
      d <- estimateGLMTrendedDisp(d,design)
      d <- estimateGLMTagwiseDisp(d,design)
      
      # glmFit
      fit <- glmFit(d,design)
      
      # make contrast to test any difference from control
      #mc <- makeContrasts(contrasts=c("Intercept", "Soil.Moisture"),levels=c("Intercept", "Soil.Moisture"))
      
      mc <- diag(rep(1, length(models.comb[,m])+1))    
      rownames(mc) <- colnames(mc) <- c("Intercept", models.comb[,m])
      
      LRT <- list()
      for(v in 1:length(models.comb[,m])){
        # v=1
        cat(" Testing: ", colnames(mc)[1+v], "\n")
        lrt <- glmLRT(fit, contrast=mc[,1+v]) 
        LRT[[colnames(mc)[1+v]]] <- lrt[[c("table")]]  
      }
      
      models.results[[paste(models.comb[,m], collapse="+")]] <- list(COEFFS = fit$coefficients, DEV = fit$deviance, DES = fit$design, LRT = LRT)
      
      #str(models.results)
      
    }   
  }
  
  dir.create(out.path, recursive=T, showWarnings = FALSE)
  save(models.results, file=paste(out.path, "/", out.name,".RData", sep=""))
  invisible(models.results)
  
}


### RUN run.edgeR.models.fit()

varialbs=c("Soil.Moisture", "Water.Potential", "Temperature")
elim.samps=c(control.samps, flowered.samps)
out.path="Models"
out.name="Models_fitting_NoControl"

Mf <- run.edgeR.models.fit(x, new.samps, varialbs=c("Soil.Moisture", "Water.Potential", "Temperature"), elim.samps=c(control.samps, flowered.samps), out.path="Models", out.name="Models_fitting_NoControl")

str(Mf)
summary(Mf)



Mf28 <- run.edgeR.models.fit(x, new.samps, varialbs=c("Soil.Moisture28", "Water.Potential28", "Temperature28"), elim.samps=c(control.samps, flowered.samps), out.path="Models", out.name="Models_fitting_28_NoControl")



# head(Mf[[1]]$LRT[[1]])
# genes <- row.names(Mf[[1]]$LRT[[1]])
# g <- which(genes=="GID059315_3538405") # WP
# dev <- vector()
# dev[1] <- Mf[["Soil.Moisture+Water.Potential+Temperature"]]$DEV[g]
# names(Mf[["Soil.Moisture+Water.Potential+Temperature"]]$LRT)
# which.max(c(Mf[["Soil.Moisture+Water.Potential+Temperature"]]$LRT[["Soil.Moisture"]][g,"PValue"], 
#           Mf[["Soil.Moisture+Water.Potential+Temperature"]]$LRT[["Water.Potential"]][g,"PValue"], 
#           Mf[["Soil.Moisture+Water.Potential+Temperature"]]$LRT[["Temperature"]][g,"PValue"]))
# 
# dev[2] <- Mf[["Soil.Moisture+Water.Potential"]]$DEV[g]
# ChiS <- diff(dev)
# pchisq(q=ChiS, df=1, lower.tail = FALSE)
# names(Mf[["Soil.Moisture+Water.Potential"]]$LRT)
# which.max(c(Mf[["Soil.Moisture+Water.Potential"]]$LRT[["Soil.Moisture"]][g,"PValue"], 
#             Mf[["Soil.Moisture+Water.Potential"]]$LRT[["Water.Potential"]][g,"PValue"]))
# dev[3] <- Mf[["Water.Potential"]]$DEV[g]
# ChiS <- diff(dev)
# pchisq(q=ChiS, df=1, lower.tail = FALSE)
# Mf[["Water.Potential"]]$LRT[["Water.Potential"]][g,"PValue"]
# Mf[["Soil.Moisture"]]$LRT[["Soil.Moisture"]][g,"PValue"]
# Mf[["Temperature"]]$LRT[["Temperature"]][g,"PValue"]
# Mf[["Water.Potential"]]$DEV[g]
# Mf[["Soil.Moisture"]]$DEV[g]
# Mf[["Temperature"]]$DEV[g]


### FUN run.edgeR.model.select

load("Models/Models_fitting_NoControl.RData")
Mf <- models.results
full.model <- "Soil.Moisture+Water.Potential+Temperature"

Mf <- Mf28
full.model <- "Soil.Moisture28+Water.Potential28+Temperature28"


genes <- row.names(Mf[[1]]$LRT[[1]])
tr <- 0.01
gene.models <- list()


for(g in 1:length(genes)){
  # g <- which(genes=="GID062074_3539215") 
  
  models.tree <- vector()
  dev <- vector()
  Chi.S <- vector()
  Chi.pv <- vector()
  
  m <- 1
  models.tree[1] <- full.model
  model.var <- unlist(strsplit(models.tree[m], split="+", fixed=TRUE))
  lrt.pvs <- matrix(NA, length(model.var), length(model.var))
  colnames(lrt.pvs) <- model.var
  
  dev[m] <- Mf[[models.tree[m]]]$DEV[g]  
  lrt.pvs.temp <- sapply(model.var, function(v) Mf[[models.tree[1]]]$LRT[[v]][g,"PValue"])
  lrt.pvs[m, names(lrt.pvs.temp)] <- lrt.pvs.temp
  
  while(max(na.omit(lrt.pvs[m,])) >= tr && m < ncol(lrt.pvs)){
    
    var.elim <- model.var[ which.max(na.omit(lrt.pvs[m,])) ]
    
    model.var <- setdiff(model.var, var.elim)
    models.tree[m+1] <- paste(model.var, collapse="+", sep="")
    
    dev[m+1] <- Mf[[models.tree[m+1]]]$DEV[g]
    
    Chi.S[m] <- dev[m+1]-dev[m]
    Chi.pv[m] <- pchisq(q=Chi.S[m], df=1, lower.tail = FALSE)
    
    lrt.pvs.temp <- sapply(model.var, function(v) Mf[[models.tree[m+1]]]$LRT[[v]][g,"PValue"])
    lrt.pvs[m+1, names(lrt.pvs.temp)] <- lrt.pvs.temp
    
    m <- m+1
  }
  
  if(m < ncol(lrt.pvs)){
    gene.models[[genes[g]]] <- list(lrt.pvs=lrt.pvs, Chi.pv=Chi.pv, final.model=models.tree[m])
  }else{
    if(na.omit(lrt.pvs[m,]) <= tr){
      gene.models[[genes[g]]] <- list(lrt.pvs=lrt.pvs, Chi.pv=Chi.pv, final.model=models.tree[m])
    }else{
      gene.models[[genes[g]]] <- list(lrt.pvs=lrt.pvs, Chi.pv=Chi.pv, final.model="Intercept")
    }
    
  }
  
} # end for


save(gene.models, file="Models/Models_selecting_NoControl.RData")

save(gene.models, file="Models/Models_selecting_28_NoControl.RData")


### version that do not select model 
for(g in 1:length(genes)){
  # g <- which(genes=="GID059315_3538405") # WP
  
  models.tree <- vector()
  dev <- vector()
  Chi.S <- vector()
  Chi.pv <- vector()
  
  models.tree[1] <- "Soil.Moisture+Water.Potential+Temperature"
  model.var <- unlist(strsplit(models.tree[1], split="+", fixed=TRUE))
  
  lrt.pvs <- matrix(NA, length(model.var), length(model.var))
  colnames(lrt.pvs) <- model.var
  
  dev[1] <- Mf[[models.tree[1]]]$DEV[g]
  
  lrt.pvs.temp <- sapply(model.var, function(v) Mf[[models.tree[1]]]$LRT[[v]][g,"PValue"])
  lrt.pvs[1, names(lrt.pvs.temp)] <- lrt.pvs.temp
  var.elim <- model.var[ which.max(na.omit(lrt.pvs[1,])) ]
  
  model.var <- setdiff(model.var, var.elim)
  models.tree[2] <- paste(model.var, collapse="+", sep="")
  
  dev[2] <- Mf[[models.tree[2]]]$DEV[g]
  
  Chi.S[1] <- dev[2]-dev[1]
  Chi.pv[1] <- pchisq(q=Chi.S[1], df=1, lower.tail = FALSE)
  
  lrt.pvs.temp <- sapply(model.var, function(v) Mf[[models.tree[2]]]$LRT[[v]][g,"PValue"])
  lrt.pvs[2, names(lrt.pvs.temp)] <- lrt.pvs.temp
  lrt.max.pv <- max(lrt.pvs)
  var.elim <- model.var[ which.max(na.omit(lrt.pvs[2,])) ]
  
  model.var <- setdiff(model.var, var.elim)
  models.tree[3] <- paste(model.var, collapse="+", sep="")
  
  dev[3] <- Mf[[models.tree[3]]]$DEV[g]
  
  Chi.S[2] <- dev[3]-dev[2]
  Chi.pv[2] <- pchisq(q=Chi.S[2], df=1, lower.tail = FALSE)
  
  lrt.pvs.temp <- sapply(model.var, function(v) Mf[[models.tree[3]]]$LRT[[v]][g,"PValue"])
  lrt.pvs[3, names(lrt.pvs.temp)] <- lrt.pvs.temp
  lrt.max.pv[3] <- max(lrt.pvs)
  
  
  gene.models[[genes[g]]] <- list(lrt.pvs=lrt.pvs, Chi.pv=Chi.pv)
  
  
}


### RUN run.egdeR()


model.var <- c("Soil.Moisture", "Water.Potential", "Temperature")


for (i in model.var){
  # i="Soil.Moisture"
  run.edgeR(x, new.samps, model.formula = as.formula(paste("~", i)), varialbs=c(i), elim.samps=flowered.samps, n.top=200, plots.path="PlotsRUNedgeR_30Aug", plot.name="", AT.id, UP.genes, DOWN.genes,  varialbs.plot = i, fit.clr = 4, LRTcoef=2)
  
  run.edgeR(x, new.samps, model.formula = as.formula(paste("~", i)), varialbs=c(i), elim.samps=c(control.samps, flowered.samps), n.top=200, plots.path="PlotsRUNedgeR_30Aug", plot.name="NoContrl", AT.id, UP.genes, DOWN.genes,  varialbs.plot = i, fit.clr = 4, LRTcoef=2)
  
}

### gene clusters by model  

load("Models/Models_fitting_NoControl.RData")
load("Models/Models_selecting_NoControl.RData")

load("Models/Models_fitting_28_NoControl.RData")
load("Models/Models_selecting_28_NoControl.RData")

genes <- names(gene.models)
models <- c(names(models.results), "Intercept")

clusters <- vector("list", length(models))
names(clusters) <- models

clusters.AT <- vector("list", length(models))
names(clusters.AT) <- models

for(g in genes){
  # g="GID000003_777"
  clusters[[gene.models[[g]]$final.model]] <- c(clusters[[gene.models[[g]]$final.model]], g)
  
  if(g %in% AT.id[,1])
  clusters.AT[[gene.models[[g]]$final.model]] <- c(clusters.AT[[gene.models[[g]]$final.model]], AT.id[g == AT.id[,1],2])
  
}

save(clusters, clusters.AT, file="Models/Models_clusters_NoControl.RData")
save(clusters, clusters.AT, file="Models/Models_clusters_28_NoControl.RData")

#####################################################################################################

### GO analysis with topGO

#####################################################################################################

load("Models/Models_fitting_NoControl.RData")
load("Models/Models_selecting_NoControl.RData")
load("Models/Models_clusters_NoControl.RData")

genes <- names(gene.models)
models <- c(names(models.results), "Intercept")

length(clusters[["Intercept"]])
length(genes)

length(clusters[["Water.Potential"]])
length(clusters[["Soil.Moisture"]])
length(clusters[["Temperature"]])
length(unique(clusters.AT[["Water.Potential"]]))
length(unique(clusters.AT[["Soil.Moisture"]]))
length(unique(clusters.AT[["Temperature"]]))


length(clusters[["Water.Potential28"]])
length(clusters[["Soil.Moisture28"]])
length(clusters[["Temperature28"]])
length(unique(clusters.AT[["Water.Potential28"]]))
length(unique(clusters.AT[["Soil.Moisture28"]]))
length(unique(clusters.AT[["Temperature28"]]))


head(MolEcol.DE.genes)

pdf("Plot_Hist/Hist_clust_WP.pdf", width = 10, height = 10)
hist(MolEcol.DE.genes[ MolEcol.DE.genes[,1] %in% clusters.AT[["Water.Potential"]], 2], breaks=7, main="Water.Potential")
dev.off()

pdf("Plot_Hist/Hist_clust_SM.pdf", width = 10, height = 10)
hist(MolEcol.DE.genes[ MolEcol.DE.genes[,1] %in% clusters.AT[["Soil.Moisture"]], 2], breaks=7, main="Soil.Moisture")
dev.off()

pdf("Plot_Hist/Hist_clust_T.pdf", width = 10, height = 10)
hist(MolEcol.DE.genes[ MolEcol.DE.genes[,1] %in% clusters.AT[["Temperature"]], 2], breaks=7, main="Temperature")
dev.off()


pdf("Plot_Hist/Hist_clust_WP28.pdf", width = 10, height = 10)
hist(MolEcol.DE.genes[ MolEcol.DE.genes[,1] %in% clusters.AT[["Water.Potential28"]], 2], breaks=7, main="Water.Potential28")
dev.off()

pdf("Plot_Hist/Hist_clust_SM28.pdf", width = 10, height = 10)
hist(MolEcol.DE.genes[ MolEcol.DE.genes[,1] %in% clusters.AT[["Soil.Moisture28"]], 2], breaks=7, main="Soil.Moisture28")
dev.off()

pdf("Plot_Hist/Hist_clust_T28.pdf", width = 10, height = 10)
hist(MolEcol.DE.genes[ MolEcol.DE.genes[,1] %in% clusters.AT[["Temperature28"]], 2], breaks=7, main="Temperature28")
dev.off()

### GO

library(topGO)
library(org.At.tair.db)

fun.gene.sel <- function(gene.vector) {
  return(gene.vector <- ifelse(gene.vector==0, FALSE, TRUE))
}

assayed.genes <- unique(unlist(clusters.AT, use.names = FALSE))

for(m in models[c(5,6,7)]){

de.genes <- unique(clusters.AT[[m]])

gene.vector=as.numeric(assayed.genes %in% de.genes)
names(gene.vector) <- assayed.genes
names(assayed.genes) <- gene.vector

for(go in c("BP","MF","CC")){

sampleGOdata <- new("topGOdata", description = "Simple session", ontology = go, allGenes = gene.vector, geneSel = fun.gene.sel , nodeSize = 10, annot = annFUN.org, mapping = "org.At.tair.db")

resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")

allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,classicKS = resultKS, elimKS = resultKS.elim, orderBy = "elimKS", ranksOf = "elimKS", topNodes = 10)

write.table(allRes, paste("GO/GO_", m , go, ".csv", sep=""), sep=";")

pdf(paste("Plots_GO/GO_", m , go, ".pdf", sep=""), width = 10, height = 10)
showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')
dev.off()

}

}








