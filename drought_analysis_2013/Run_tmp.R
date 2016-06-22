
setwd("/home/Shared/data/seq/Shimizu_RNA_seq/")
load("Shimizu_workspace.Rdata")

library(EDASeq)
library(RUVSeq)
library(edgeR)

source("/home/gosia/R/R_Shimizu_RNA_seq/Run_RUVSeq.R")

## negative control genes = housekeeping genes
neg.contrl.AT <- read.table("_Data/genes_descr_control/Housekeeping_Genes/Table4_Arabidopsis_Housekeepieng_Gene_ID.txt", header=TRUE, sep="\t")
neg.contrl.AT$Arabidopsis.Housekeeping.Gene.ID <- toupper(neg.contrl.AT$Arabidopsis.Housekeeping.Gene.ID)
neg.contrl.AT$Housekeeping <- TRUE
colnames(neg.contrl.AT) <- c("AT_id", "Housekeeping")
head(neg.contrl.AT)


neg.contrl <- merge(genes.full.description, neg.contrl.AT, by="AT_id", all.x=TRUE, sort=FALSE)
head(neg.contrl)

dim(neg.contrl)
sum(!is.na(neg.contrl$Housekeeping))

neg.contrl.ID <- neg.contrl[!is.na(neg.contrl$Housekeeping), "SL_id"]


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
FDR=0.1


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



