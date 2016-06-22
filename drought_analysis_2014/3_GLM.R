
#####################################################################################################

### GLM - models with one covariate 

#####################################################################################################

####################################
### run.edgeR()
####################################

# function

setwd("/home/Shared/data/seq/Shimizu_RNA_seq/")
load("Shimizu_workspace.Rdata")

library(edgeR)
source("/home/gosia/R/R_Shimizu_RNA_seq/Run_edgeR.R")

#### TO RUN:


model.formula = ~ Water.Potential
varialbs=c("Water.Potential")
varialbs.plot = "Water.Potential"

plots.path="Plots_RUN_edgeR"
fit.clr = 1
LRTcoef="Water.Potential"
FDR=0.05

#elim.samps = c(flowered.samps, control.samps)
elim.samps = c(flowered.samps)
plot.name="_ALL"

ee <- run.edgeR(x, new.samps, model.formula=model.formula, varialbs=varialbs, elim.samps=elim.samps, FDR=FDR, plots.path=plots.path, plot.name=plot.name, varialbs.plot=varialbs.plot, fit.clr = fit.clr, LRTcoef=LRTcoef, trees.order=trees.order, genes.full.description=genes.full.description, AT.id=AT.id)



elim.samps = c(flowered.samps, control.samps)
plot.name="_DROUGHT"

ee <- run.edgeR(x, new.samps, model.formula=model.formula, varialbs=varialbs, elim.samps=elim.samps, FDR=FDR, plots.path=plots.path, plot.name=plot.name, varialbs.plot=varialbs.plot, fit.clr = fit.clr, LRTcoef=LRTcoef, trees.order=trees.order, genes.full.description=genes.full.description, AT.id=AT.id)


elim.samps = c(flowered.samps, new.samps$sample_name[new.samps$tree_ID != 8266] )
plot.name="_8266"

ee <- run.edgeR(x, new.samps, model.formula=model.formula, varialbs=varialbs, elim.samps=elim.samps, FDR=FDR, plots.path=plots.path, plot.name=plot.name, varialbs.plot=varialbs.plot, fit.clr = fit.clr, LRTcoef=LRTcoef, trees.order=trees.order, genes.full.description=genes.full.description, AT.id=AT.id)



elim.samps = c(flowered.samps, new.samps$sample_name[new.samps$tree_ID != 970] )
plot.name="_970"

ee <- run.edgeR(x, new.samps, model.formula=model.formula, varialbs=varialbs, elim.samps=elim.samps, FDR=FDR, plots.path=plots.path, plot.name=plot.name, varialbs.plot=varialbs.plot, fit.clr = fit.clr, LRTcoef=LRTcoef, trees.order=trees.order, genes.full.description=genes.full.description, AT.id=AT.id)




new.samps$tree_ID <- as.factor(new.samps$tree_ID)

model.formula = ~ Water.Potential + tree_ID - 1
elim.samps = c(flowered.samps)
plot.name="_ALL"


ee <- run.edgeR(x, new.samps, model.formula=model.formula, varialbs=varialbs, elim.samps=elim.samps, FDR=FDR, plots.path=plots.path, plot.name=plot.name, varialbs.plot=varialbs.plot, fit.clr = fit.clr, LRTcoef=LRTcoef, trees.order=trees.order, genes.full.description=genes.full.description, AT.id=AT.id)









# old runs: with genes as covariates 

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



####################################
### run.voom()
####################################

# function

setwd("/home/Shared/data/seq/Shimizu_RNA_seq/")
load("Shimizu_workspace.Rdata")


source("/home/gosia/R/R_Shimizu_RNA_seq/Run_voom.R")


model.formula = ~ Water.Potential
varialbs=c("Water.Potential")
varialbs.plot = "Water.Potential"

plots.path="Plots_RUN_voom"
fit.clr = 1
LRTcoef="Water.Potential"
FDR=0.05
# voom.method = "robust"
voom.method="ls" 

#elim.samps = c(flowered.samps, control.samps)
elim.samps = c(flowered.samps)
plot.name="_ALL_ls"


tt <- run.voom(x, new.samps, model.formula=model.formula, varialbs=varialbs, elim.samps=elim.samps, FDR=FDR, plots.path=plots.path, plot.name=plot.name, varialbs.plot=varialbs.plot, fit.clr = fit.clr, LRTcoef=LRTcoef, trees.order=trees.order, genes.full.description=genes.full.description, AT.id=AT.id, voom.method=voom.method)




model.formula = ~ Water.Potential
varialbs=c("Water.Potential")
varialbs.plot = "Water.Potential"

plots.path="Plots_RUN_voom"
fit.clr = 1
LRTcoef="Water.Potential"
FDR=0.05
voom.method = "robust"
# voom.method="ls" 

elim.samps = c(flowered.samps, control.samps)
plot.name="_DROUGHT_robust"

tt <- run.voom(x, new.samps, model.formula=model.formula, varialbs=varialbs, elim.samps=elim.samps, FDR=FDR, plots.path=plots.path, plot.name=plot.name, varialbs.plot=varialbs.plot, fit.clr = fit.clr, LRTcoef=LRTcoef, trees.order=trees.order, genes.full.description=genes.full.description, AT.id=AT.id, voom.method=voom.method)




model.formula = ~ Water.Potential
varialbs=c("Water.Potential")
varialbs.plot = "Water.Potential"

plots.path="Plots_RUN_voom"
fit.clr = 1
LRTcoef="Water.Potential"
FDR=0.05
# voom.method = "robust"
voom.method="ls" 

elim.samps = c(flowered.samps, new.samps$sample_name[new.samps$tree_ID != 8266] )
plot.name="_8266_ls"

tt <- run.voom(x, new.samps, model.formula=model.formula, varialbs=varialbs, elim.samps=elim.samps, FDR=FDR, plots.path=plots.path, plot.name=plot.name, varialbs.plot=varialbs.plot, fit.clr = fit.clr, LRTcoef=LRTcoef, trees.order=trees.order, genes.full.description=genes.full.description, AT.id=AT.id, voom.method=voom.method)




model.formula = ~ Water.Potential
varialbs=c("Water.Potential")
varialbs.plot = "Water.Potential"

plots.path="Plots_RUN_voom"
fit.clr = 1
LRTcoef="Water.Potential"
FDR=0.05
# voom.method = "robust"
voom.method="ls" 

elim.samps = c(flowered.samps, new.samps$sample_name[new.samps$tree_ID != 970] )
plot.name="_970_ls"

tt <- run.voom(x, new.samps, model.formula=model.formula, varialbs=varialbs, elim.samps=elim.samps, FDR=FDR, plots.path=plots.path, plot.name=plot.name, varialbs.plot=varialbs.plot, fit.clr = fit.clr, LRTcoef=LRTcoef, trees.order=trees.order, genes.full.description=genes.full.description, AT.id=AT.id, voom.method=voom.method)





new.samps$tree_ID <- as.factor(new.samps$tree_ID)

model.formula = ~ Water.Potential + tree_ID - 1
varialbs=c("Water.Potential")
varialbs.plot = "Water.Potential"

plots.path="Plots_RUN_voom"
fit.clr = 1
LRTcoef="Water.Potential"
FDR=0.05
voom.method = "robust"
# voom.method="ls" 

#elim.samps = c(flowered.samps, control.samps)
elim.samps = c(flowered.samps)
plot.name="_ALL_robust"

tt <- run.voom(x, new.samps, model.formula=model.formula, varialbs=varialbs, elim.samps=elim.samps, FDR=FDR, plots.path=plots.path, plot.name=plot.name, varialbs.plot=varialbs.plot, fit.clr = fit.clr, LRTcoef=LRTcoef, trees.order=trees.order, genes.full.description=genes.full.description, AT.id=AT.id, voom.method=voom.method)




####################################
### run.edgeR.robust()
####################################
# BioC 14

setwd("/home/Shared/data/seq/Shimizu_RNA_seq/")
load("Shimizu_workspace.Rdata")

library(edgeR)

source("/home/gosia/R/R_Shimizu_RNA_seq/Run_edgeR_robust.R")


model.formula = ~ Water.Potential
varialbs=c("Water.Potential")
varialbs.plot = "Water.Potential"

# model.formula = ~ Water.Potential14
# varialbs=c("Water.Potential14")
# varialbs.plot = "Water.Potential14"
# 
# model.formula = ~ Soil.Moisture
# varialbs=c("Soil.Moisture")
# varialbs.plot = "Soil.Moisture"


# per drought tree
plots.path="Plots_RUN_edgeR_robust"
fit.clr = 1
LRTcoef = "Water.Potential"
FDR = 0.1


elim.samps = c(flowered.samps)
plot.name="_ALL"

try(run.edgeR.robust(x, new.samps, model.formula, varialbs, elim.samps, FDR=FDR, plots.path=plots.path, plot.name=plot.name, varialbs.plot, fit.clr = fit.clr, LRTcoef=LRTcoef, trees.order=trees.order, genes.full.description=genes.full.description, AT.id=AT.id))



elim.samps = c(flowered.samps, control.samps)
plot.name="_DROUGHT"


try(run.edgeR.robust(x, new.samps, model.formula, varialbs, elim.samps, FDR=FDR, plots.path=plots.path, plot.name=plot.name, varialbs.plot, fit.clr = fit.clr, LRTcoef=LRTcoef, trees.order=trees.order, genes.full.description=genes.full.description, AT.id=AT.id))



elim.samps = c(flowered.samps, new.samps$sample_name[new.samps$tree_ID != 8266] )
plot.name="_8266"

try(run.edgeR.robust(x, new.samps, model.formula, varialbs, elim.samps, FDR=FDR, plots.path=plots.path, plot.name=plot.name, varialbs.plot, fit.clr = fit.clr, LRTcoef=LRTcoef, trees.order=trees.order, genes.full.description=genes.full.description, AT.id=AT.id))



elim.samps = c(flowered.samps, new.samps$sample_name[new.samps$tree_ID != 970] )
plot.name="_970"


try(run.edgeR.robust(x, new.samps, model.formula, varialbs, elim.samps, FDR=FDR, plots.path=plots.path, plot.name=plot.name, varialbs.plot, fit.clr = fit.clr, LRTcoef=LRTcoef, trees.order=trees.order, genes.full.description=genes.full.description, AT.id=AT.id))






####################################
### run.RUV.edgeR()
####################################


# install.packages("/home/gosia/R/packages/RUVSeq/EDASeq_1.9.3.tar.gz", dependencies=TRUE, lib="/home/gosia/R/libraries/3.0.2")
# install.packages("/home/gosia/R/packages/RUVSeq/RUVSeq_0.99.1.tar.gz", dependencies=TRUE, lib="/home/gosia/R/libraries/3.0.2")
# install.packages("/home/gosia/R/packages/RUVSeq/zebrafishRNASeq_0.99.1.tar.gz", dependencies=TRUE, lib="/home/gosia/R/libraries/3.0.2")

setwd("/home/Shared/data/seq/Shimizu_RNA_seq/")
load("Shimizu_workspace.Rdata")

library(EDASeq)
library(RUVSeq)
library(edgeR)

source("/home/gosia/R/R_Shimizu_RNA_seq/Run_RUVSeq.R")

## negative control genes = housekeeping genes

neg.contrl.AT <- read.table("Data/genes_descr_control/Housekeeping_Genes/Table4_Arabidopsis_Housekeepieng_Gene_ID.txt", header=TRUE, sep="\t")

neg.contrl.AT$Arabidopsis.Housekeeping.Gene.ID <- toupper(neg.contrl.AT$Arabidopsis.Housekeeping.Gene.ID)
neg.contrl.AT$Housekeeping <- TRUE
head(neg.contrl.AT)


neg.contrl <- merge(genes.full.description, neg.contrl.AT, by.x="AT_ID", by.y="Arabidopsis.Housekeeping.Gene.ID", all.x=TRUE, sort=FALSE)
head(neg.contrl)

dim(neg.contrl)
sum(!is.na(neg.contrl$Housekeeping))

neg.contrl.ID <- neg.contrl[!is.na(neg.contrl$Housekeeping), "ID"]


# RUN RUVSeq normalisation 

elim.samps = c(flowered.samps)
plots.path="Plots_RUN_RUVSeq_normalization"
plot.name="RUVg_norm"


set1 <- run.RUVSeq(x, new.samps, neg.contrl.ID, elim.samps, plots.path=plots.path, plot.name=plot.name)


# RUN DE

model.formula = ~ Water.Potential + W_1
varialbs=c("Water.Potential")
varialbs.plot = "Water.Potential"


plots.path="Plots_RUN_edgeR_RUVSeq"
fit.clr = 1
LRTcoef = "Water.Potential"
FDR=0.05


elim.samps = c(flowered.samps)
plot.name="_ALL"

rr <- run.edgeR.RUVSeq(set1, new.samps, model.formula=model.formula, varialbs=varialbs, elim.samps=elim.samps, FDR=FDR, plots.path=plots.path, plot.name=plot.name, varialbs.plot=varialbs.plot, fit.clr = fit.clr, LRTcoef=LRTcoef, trees.order=trees.order, genes.full.description=genes.full.description, AT.id=AT.id)



elim.samps = c(flowered.samps, control.samps)
plot.name="_DROUGHT"


rr <- run.edgeR.RUVSeq(set1, new.samps, model.formula=model.formula, varialbs=varialbs, elim.samps=elim.samps, FDR=FDR, plots.path=plots.path, plot.name=plot.name, varialbs.plot=varialbs.plot, fit.clr = fit.clr, LRTcoef=LRTcoef, trees.order=trees.order, genes.full.description=genes.full.description, AT.id=AT.id)


elim.samps = c(flowered.samps, new.samps$sample_name[new.samps$tree_ID != 8266] )
plot.name="_8266"


rr <- run.edgeR.RUVSeq(set1, new.samps, model.formula=model.formula, varialbs=varialbs, elim.samps=elim.samps, FDR=FDR, plots.path=plots.path, plot.name=plot.name, varialbs.plot=varialbs.plot, fit.clr = fit.clr, LRTcoef=LRTcoef, trees.order=trees.order, genes.full.description=genes.full.description, AT.id=AT.id)



elim.samps = c(flowered.samps, new.samps$sample_name[new.samps$tree_ID != 970] )
plot.name="_970"


rr <- run.edgeR.RUVSeq(set1, new.samps, model.formula=model.formula, varialbs=varialbs, elim.samps=elim.samps, FDR=FDR, plots.path=plots.path, plot.name=plot.name, varialbs.plot=varialbs.plot, fit.clr = fit.clr, LRTcoef=LRTcoef, trees.order=trees.order, genes.full.description=genes.full.description, AT.id=AT.id)



