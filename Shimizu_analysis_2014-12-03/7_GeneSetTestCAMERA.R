
#####################################################################################################

### CAMERA Gene set testing

#####################################################################################################

setwd("/home/Shared/data/seq/Shimizu_RNA_seq/")
load("Shimizu_workspace.Rdata")

library(edgeR)

####

model.formula = ~ Water.Potential28
varialbs=c("Water.Potential28")
varialbs.plot = "Water.Potential28"

plots.path="Plots_RUN_camera"

LRTcoef="Water.Potential28"


elim.samps = c(flowered.samps)
plot.name="_ALL"

####



####




dir.create(plots.path, recursive=T, showWarnings = FALSE)

library(stringr)
model.char <- str_replace_all(as.character(model.formula), " ", "") # can be gsub()
model.char <- str_replace_all(model.char[2], "\\.", "_")


x.tmp <- x[, new.samps$sample_name]  
x.tmp <- x.tmp[,!names(x.tmp) %in% elim.samps]
new.samps.tmp <- new.samps[!new.samps$sample_name %in% elim.samps, ]


  for(i in varialbs){
    x.tmp <- x.tmp[,!is.na(new.samps.tmp[, i])]
    new.samps.tmp <- new.samps.tmp[!is.na(new.samps.tmp[, i]), ]
  } 

d <- DGEList(x.tmp, group=new.samps.tmp$tree_ID)
d <- calcNormFactors(d)

### make sure a gene is expressed (CPM > 1) in more than 2 samples
d.cpm <- cpm(d, normalized.lib.sizes=TRUE)


d <- d[ rowSums(d.cpm>1) > ncol(d.cpm)/2, ]
# cat(paste("*", nrow(d), "genes expressed in more than two samples \n"))


# design model matrix
design <- model.matrix(model.formula, data=new.samps.tmp)
design

# voom transformation  
v <- voom(d,design,plot=FALSE)


indexFlowering <- rownames(v) %in% rownames( genes.full.description[!is.na(genes.full.description[, "Flowering"]), ] )

genesDrounght_regulation <- genes.full.description[!is.na(genes.full.description[, "Drought_regulation"]) , ] 

indexDrounght_Up <- rownames(v) %in% rownames( genesDrounght_regulation[ genesDrounght_regulation[, "Drought_regulation"] == "Up", ] )

indexDrounght_Down <- rownames(v) %in% rownames( genesDrounght_regulation[ genesDrounght_regulation[, "Drought_regulation"] == "Down", ] )

index <- list(indexFlowering, indexDrounght_Up, indexDrounght_Down)



camera(v,index,design)

