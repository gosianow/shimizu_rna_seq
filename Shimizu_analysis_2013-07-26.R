#####################################################################################################

# data load

#####################################################################################################

setwd("/Users/gosia/DATA/Shimizu_RNA_seq/")

### get table of counts

# ch <- as.character(read.table("data/raw_data.csv", sep=",",nrow=1,stringsAsFactors=FALSE)[1,-c(1,2,57)])
# x <- read.table("data/raw_data.csv", sep=",",skip=3, header=F)
# x <- x[,-c(2,57)]
# colnames(x) <- c("genes", ch)
# write.table(x, "data/raw_data-clean.csv", sep=",", row.names = F, quote=F)

x <- read.table("data/raw_data-clean.csv", sep=",", header=T, row.names=1)
head(x)
names(x)

### parse sample information
samps <- read.table("data/sample_list.csv", sep=";", header=TRUE)
head(samps)
names(samps)
#write.table(samps, "sample_list.xls", sep="\t", row.names = F)

### time format change 

samps$time_ch <- strptime(samps$year.month.day, "%d.%m.%y")
samps$time_nr <- as.numeric(samps$time_ch)

samps$flowered <- samps$tree_ID=="8266"
samps$tree_ID <- as.factor(samps$tree_ID)

names(samps)

#####################################################################################################

### analysis by MARK

#####################################################################################################

### edgeR

library(edgeR)
d <- DGEList(x, group=samps$tree_ID)
d <- calcNormFactors(d)


### make sure a gene is expressed (CPM > 1) in more than 2 samples
cps <- cpm(d, normalized.lib.sizes=TRUE)
d <- d[ rowSums(cps>1) > 2, ]

### reorder by tree and time
o <- order(d$samples$group, samps$tm)
samps <- samps[o,]
d <- d[,o]
tmO <- tmO[o]
samps$tm <- samps$tm[o]

### create d$genes -> TAIR annotation / AT ID / functional description

blast <- read.table("data/best_hit_blast_result_Sl_predicted_exons.txt",sep=",")
m <- match( rownames(d$counts), blast$V1)
d$genes <- data.frame(assembl_id=rownames(d$counts), at_id=as.character(blast$V2[m]), at_symbol=NA, 
                      stringsAsFactors=FALSE)


library(org.At.tair.db)
nna <- !is.na(d$genes$at_id)
syms <- mget(d$genes$at_id[nna], org.At.tairSYMBOL)
syms <- sapply(syms, function(u) if(is.na(u[1])) NA else paste(u,collapse=";"))
d$genes$at_symbol[nna] <- syms

tfd <- read.table("data/TAIR10_functional_descriptions",sep="\t", header=TRUE, quote="", 
                  comment.char="", stringsAsFactors=FALSE)
t10id <- gsub("\\.[1-9]","",tfd$Model_name)
m <- match(d$genes$at_id[nna],t10id)
d$genes$description[nna] <- as.character(tfd$Short_description[m])

# write.table(d$genes, "data/genes_description.xls", sep="\t", row.names = F, quote=F)

cps <- cpm(d, normalized.lib.sizes=TRUE)
head(cps)
d$genes <- cbind(d$genes, round(cps,1))

head(d$genes)


### plots: MDS

library(ggplot2)
library(ggdendro)

cols <- c("black","blue","orange","darkgreen","red","salmon")

pdf("out/mds3.pdf",w=10,h=10)
mds <- plotMDS(d, col=cols[as.numeric(d$samples$group)],main="method=logFC, prior.count=2")
mds10 <- plotMDS(d,col=cols[as.numeric(d$samples$group)],main="method=logFC, prior.count=10",
                 prior.count=10)
mdsb <- plotMDS(d,col=cols[as.numeric(d$samples$group)],method="bcv")
hc <- hclust(as.dist(mdsb$distance.matrix))
#plot(hc, col=as.numeric(d$samples$group), hang=-1)
ggdendrogram(hc, rotate=TRUE, size=4, theme_dendro=FALSE, 
             color=(cols[as.numeric(d$samples$group)])[hc$order])
dev.off()


### plots: smoothScatter

p <- function(...) smoothScatter(..., nrpoints=0, nbin=64, add=TRUE, transformation = function(x) x^.1)
g <- unique(d$samples$group)
pdf("out/all.pdf",w=10,h=10)
for(i in g)
pairs( log2(1+d$counts[,as.character(d$samples$group)==i]), pch=".", panel=p, lower.panel=NULL)
dev.off()


### take subset
dbak <- d

keep <- samps$tree_ID != "990" & samps$developmental_stage != "flower_bud"
d <- d[,keep]

### design model matrix

design <- model.matrix(~-1 + tree_ID, data=samps[keep,])

#design <- model.matrix(~0+tree_ID, data=samps[keep,])

design <- design[, colSums(design)>0]

# estimate dispersion

d <- estimateGLMCommonDisp(d,design)
d <- estimateGLMTrendedDisp(d,design)
d <- estimateGLMTagwiseDisp(d,design)

### glmFit fits genewise negative binomial glms

fit <- glmFit(d,design)

# make contrast to test any difference from control
mc <- makeContrasts(contrasts=c("(tree_ID970+tree_ID8212)/2-(tree_ID1099+tree_ID1377)/2", 
                                "tree_ID8266-(tree_ID1099+tree_ID1377)/2"),levels=colnames(design))
colnames(mc) <- c("sheet_vs_control","flower_vs_control")

# glmLRT conducts likelihood ratio tests for one or more coefficients in the linear model

lrt <- glmLRT(fit, contrast=mc)

tt <- topTags(lrt,n=nrow(d))$table
write.table(tt,"out/crude_global_analysis_12june2013.xls", row.names=FALSE, sep="\t")

ttt <- topTags(lrt,n=1000)$table
write.table(ttt[,c("assembl_id", "logFC.sheet_vs_control", "logFC.flower_vs_control" ,"logCPM",
                   "LR", "PValue", "FDR" )],"out/TOP1000mc.xls", row.names=FALSE, sep="\t")


### plot expression levels over time for TOP DE genes

plotGene <- function(d,id="GID000013_4400", tm, cols=c("black","blue","blue","blue","orange","black"), 
                     FUN=sqrt,main=NULL,...) 
{
  cps <- cpm(d, normalized.lib.sizes=TRUE)
  x <- cps[id,]
  g <- d$samples$group
  ug <- unique(as.character(g))
  # main title for plot
  if (is.null(main)) {
    main <- id
    m <- match(id, d$genes$assembl_id)
    if( !is.na(d$genes$at_id[m]) ) main <- paste0(main," // ",d$genes$at_id[m])
    if( nchar(d$genes$at_symbol[m]) > 0 ) main <- paste0(main," // ",d$genes$at_symbol[m])
  }
  for(tr in 1:length(ug)) {
    k <- g==ug[tr]
    if(tr==1) {
      plot(tm[k], FUN(x[k]), pch=19, lwd=3, type="b", ylim=range(sqrt(x)),col=cols[tr], 
           xlab="Date", ylab="expression (counts per million)",main=main,...)
    } else {
      points(tm[k], FUN(x[k]), pch=19, lwd=3, type="b",col=cols[tr])
    }
  }
  
  dt <- as.numeric(as.POSIXlt("2009-01-01"))
  abline(v=dt,lwd=4,col="grey")
  abline(v=as.numeric(as.POSIXlt(paste0("2009-",1:12,"-01"))),lwd=1,col="grey")
  abline(v=as.numeric(as.POSIXlt(paste0("2008-",11:12,"-01"))),lwd=1,col="grey")
  
  legend("topright",as.character(ug),col=cols,pch=19,lwd=3,bg="white")
}


rn <- topTags(lrt,n=120)$table$assembl_id # top 120 genes

pdf("out/differential_drought_vs_control_13june2013.pdf",w=18,h=12)
par(mfrow=c(2,3))
for(ii in rn)
  plotGene(d,ii,tm[keep],cols=c("grey30","darkblue","blue","grey50","orange"))
dev.off()



### additional contrast to check

mc1 <- makeContrasts(contrasts=c("(tree_ID970+tree_ID8212)/2-(tree_ID1099+tree_ID1377)/2"),
                     levels=colnames(design))
colnames(mc1) <- c("sheet_vs_control")

lrt1 <- glmLRT(fit, contrast=mc1)

ttt1 <- topTags(lrt1,n=1000)$table
names(ttt1)
write.table(ttt1[,c("assembl_id", "logFC" ,"logCPM",
                   "LR", "PValue", "FDR" )],"out/TOP1000mc1.xls", row.names=FALSE, sep="\t")


mc2 <- makeContrasts(contrasts=c("tree_ID8266-(tree_ID1099+tree_ID1377)/2"),levels=colnames(design))
colnames(mc2) <- c("flower_vs_control")

lrt2 <- glmLRT(fit, contrast=mc2)

ttt2 <- topTags(lrt2,n=1000)$table
write.table(ttt2[,c("assembl_id", "logFC","logCPM",
                    "LR", "PValue", "FDR" )],"out/TOP1000mc2.xls", row.names=FALSE, sep="\t")


#####################################################################################################

### preparation of variables describing samples - Soild Moisture & Water Potential

#####################################################################################################

setwd("/Users/gosia/DATA/Shimizu_RNA_seq/")

####################################################
### samps
####################################################

### parse sample information
samps <- read.table("data/sample_list.csv", sep=";", stringsAsFactors=FALSE, header=TRUE)
head(samps)
names(samps)

# no genes variables &  no samps from November 2009
new.samps <- samps[!samps$sample_name%in% c("N9_970_20090114","M1_1099_20090114","I5_8212_20090119","F8_8266_20090119"), c(1:16)]

### time format change 
new.samps$time_ch <- strptime(new.samps$year.month.day, "%d.%m.%y")
new.samps$time_nr <- as.numeric(new.samps$time_ch)

new.samps$flowered <- new.samps$tree_ID=="8266"

names(new.samps)

### list of all unique days in 2008 and 2009
ad <- read.table("Data/Unique_days_short.csv", sep=";")
all.days <- as.numeric(strptime(ad[1:731,], "%d.%m.%y"))


# columns 8-16 come from Lambir_meteorological_data_...xls


# FUN that matches two tables 
match.func <- function(DF1, col.match1, DF2, col.match2, col.values2){
  
  DF1[, col.values2] <- NA
  
  for(i in unique(DF1[,col.match1])){
    # i=unique(DF1[,col.match1])[1]
    DF1[DF1[,col.match1]==i, col.values2] <- ifelse(length(DF2[DF2[,col.match2]==i, col.values2])==0, NA, DF2[DF2[,col.match2]==i, col.values2])
    
  }
  invisible(return(DF1))
}

### FUN that fits loess and calculates rolled means 
loess.and.roll <- function(intrp.x, intrp.y, intrp.points, table1, col.values="Soil.Moisture", plots.path="PlotsLoess", plot.name="", span=0.2){
  #table1 has to have "time_nr" column
  
  # FUN that matches two tables 
  match.func <- function(DF1, col.match1, DF2, col.match2, col.values2){
    
    DF1[, col.values2] <- NA
    
    for(i in unique(DF1[,col.match1])){
      # i=unique(DF1[,col.match1])[1]
      DF1[DF1[,col.match1]==i, col.values2] <- ifelse(length(DF2[DF2[,col.match2]==i, col.values2])==0, NA, DF2[DF2[,col.match2]==i, col.values2])
      
    }
    invisible(DF1)
  }
  
  intrp.points <- sort(intrp.points)
  
  data <- data.frame(intrp.x, intrp.y)
  data <- data[order(data[,1]), ]
  colnames(data) <- c("time_nr", col.values)
  
  ### loess
  table2.loess <- loess(as.formula(paste(col.values, "~", "time_nr")), data, span = span, control = loess.control(surface = "direct"))
  table2.pred <- predict(table2.loess , intrp.points, se = TRUE)

  table2 <- data.frame(intrp.points, table2.pred$fit)
  colnames(table2) <- c("time_nr", col.values)
  
  ### rollmean over 14 and 28 days
  library(zoo)
  col.values14 <- paste(col.values, "14", sep="")
  table2[,col.values14] <- NA
  table2[14:nrow(table2), col.values14] <- rollmean(table2[,col.values], 14)
  
  col.values28 <- paste(col.values, "28", sep="")
  table2[,col.values28] <- NA
  table2[28:nrow(table2), col.values28] <- rollmean(table2[,col.values], 28)
  
  table1 <- match.func(DF1=table1, col.match1="time_nr", DF2=table2, col.match2="time_nr", col.values2=col.values)
  table1 <- match.func(DF1=table1, col.match1="time_nr", DF2=table2, col.match2="time_nr", col.values2=col.values14)
  table1 <- match.func(DF1=table1, col.match1="time_nr", DF2=table2, col.match2="time_nr", col.values2=col.values28)
  
  dir.create(plots.path, recursive=T, showWarnings=FALSE)
  
  library(stringr)
  
  pdf( paste(plots.path, "/Loess_", str_replace_all(paste(col.values, plot.name, sep=""),"[[:punct:]]", "_"),".pdf", sep=""), width = 10, height = 5)
  plot(intrp.x, intrp.y, pch=20, col=1, main=plot.name ,ylab=col.values, xlab="Time", xlim=c(min(table1[,"time_nr"], intrp.x), max(table1[,"time_nr"], intrp.x)), ylim=c(min(table2.pred$fit, intrp.y), max(table2.pred$fit, intrp.y)))
  abline(v=table1[,"time_nr"], col="grey")
  lines(intrp.points, table2.pred$fit, col=5)
  points(table1[,"time_nr"], table1[,col.values], col=4)
  dev.off()
  
  invisible(table1)
  
}


### FUN that interpolates and calculates rolled means 
interpolate.and.roll <- function(intrp.x, intrp.y, intrp.points, table1, col.values, plots.path="PlotsInterpolate", plot.name=""){
  #table1 has to have "time_nr" column
  
  # FUN that matches two tables 
  match.func <- function(DF1, col.match1, DF2, col.match2, col.values2){
    
    DF1[, col.values2] <- NA
    
    for(i in unique(DF1[,col.match1])){
      # i=unique(DF1[,col.match1])[1]
      DF1[DF1[,col.match1]==i, col.values2] <- ifelse(length(DF2[DF2[,col.match2]==i, col.values2])==0, NA, DF2[DF2[,col.match2]==i, col.values2])
      
    }
    invisible(DF1)
  }
  
  library(zoo)
  
  interp.intrp.points <- approx(intrp.x, intrp.y, xout=intrp.points)
  
  table2<-  data.frame(time_nr=interp.intrp.points$x[!is.na(interp.intrp.points$y)], interp.intrp.points$y[!is.na(interp.intrp.points$y)])
  
  colnames(table2) <- c("time_nr", col.values)
  
  ### rollmean over 14 and 28 days
  col.values14 <- paste(col.values, "14", sep="")
  table2[,col.values14] <- NA
  table2[14:nrow(table2), col.values14] <- rollmean(table2[,col.values], 14)
  
  col.values28 <- paste(col.values, "28", sep="")
  table2[,col.values28] <- NA
  table2[28:nrow(table2), col.values28] <- rollmean(table2[,col.values], 28)
  
  table1 <- match.func(DF1=table1, col.match1="time_nr", DF2=table2, col.match2="time_nr", col.values2=col.values)
  table1 <- match.func(DF1=table1, col.match1="time_nr", DF2=table2, col.match2="time_nr", col.values2=col.values14)
  table1 <- match.func(DF1=table1, col.match1="time_nr", DF2=table2, col.match2="time_nr", col.values2=col.values28)
  
  dir.create(plots.path, recursive=T, showWarnings=FALSE)
  
  library(stringr)
  
  pdf( paste(plots.path, "/Interp_", str_replace_all(paste(col.values, plot.name, sep=""),"[[:punct:]]", "_"),".pdf", sep=""), width = 10, height = 5)
  plot(intrp.x, intrp.y, type="l", col=6, ylab=col.values, xlab="Time")
  abline(v=table1[,"time_nr"], col="grey")
  points(interp.intrp.points, col=1, pch=20)
  points(table1[,"time_nr"], table1[,col.values])
  dev.off()
  
  
  invisible(table1)
  
}


####################################################
### Soild Moisture
####################################################

### get soil moisture data
sm <- read.table("data/meteor/Soil_Moisture.csv", sep=";", stringsAsFactors=FALSE, header=TRUE)
sm <- sm[1:271,]

### time format change
sm$Year <- ifelse(sm$DayUniq <= 366, 2008, 2009)
sm$year.month.day <- paste(sm$Day,sm$Month, sm$Year, sep=".")
sm$time_ch <- strptime(sm$year.month.day, "%d.%m.%Y")
sm$time_nr <- as.numeric(sm$time_ch)

### plots of raw data

pdf("Plots/Soil_moisture.pdf", width = 10, height = 5)
plot(sm$time_nr, sm$DE970, type="l", col=1, ylim=c(-0.2,1.1), main="Soil Moisture", ylab="Soil moisture", xlab="Time")
lines(sm$time_nr, sm$C970, col=1, lty=2)
lines(sm$time_nr, sm$DE8266, col=3)
lines(sm$time_nr, sm$C8266, col=3, lty=2)
lines(sm$time_nr, sm$DE8212, col=4)
lines(sm$time_nr, sm$C8212a, col=4, lty=2)
lines(sm$time_nr, sm$C8212b, col=4, lty=3)
legend("bottomleft", c("DE970", "C970","DE8266","C8266" ,"DE8212","C8212a" ,"C8212b"), col=c(1, 1, 3, 3, 4, 4, 4), lty=c(1,2,1,2,1,2,3), cex=0.5)
dev.off()


pdf("Plots/Soil_moisture_Control.pdf", width = 10, height = 5)
plot(sm$time_nr, sm$C970, type="l", col=2, ylim=c(0,1.1), xlim=c(1.220e+09, 1257289200), main="Control trees", ylab="Soil moisture", xlab="Time")
lines(sm$time_nr, sm$C8266, col=3)
lines(sm$time_nr, sm$C8212a, col=4)
lines(sm$time_nr, sm$C8212b, col=5)
abline(v=new.samps$time_nr[new.samps$tree_ID %in% c(990, 1099, 1377)], col="grey")
dev.off()

pdf("Plots/Soil_moisture_DE970.pdf", width = 10, height = 5)
plot(sm$time_nr, sm$DE970, type="l", col=2, ylim=c(0,1.1), xlim=c(1.220e+09, 1257289200), main="DE970", ylab="Soil moisture", xlab="Time")
abline(v=samps$time_nr[samps$tree_ID==970], col="grey")
dev.off()

pdf("Plots/Soil_moisture_DE8266.pdf", width = 10, height = 5)
plot(sm$time_nr, sm$DE8266, type="l", col=2, ylim=c(0,1.1), xlim=c(1.220e+09, 1257289200), main="DE8266", ylab="Soil moisture", xlab="Time")
abline(v=samps$time_nr[samps$tree_ID==8266], col="grey")
dev.off()

pdf("Plots/Soil_moisture_DE8212.pdf", width = 10, height = 5)
plot(sm$time_nr, sm$DE8212, type="l", col=2, ylim=c(0,1.1), xlim=c(1.220e+09, 1257289200), main="DE8212", ylab="Soil moisture", xlab="Time")
abline(v=samps$time_nr[samps$tree_ID==8212], col="grey")
dev.off()

# plots with intepolation

pdf("Plots/Soil_moisture_Control_lowess_interpolation.pdf", width = 10, height = 5)
plot(control.l, col=2, ylim=c(0,1.1), type="l", main="Control trees", ylab="Soil moisture", xlab="Time")
abline(v=new.samps$time_nr[new.samps$tree_ID %in% c(990, 1099, 1377)], col="grey")
Control.interp <- approx(control.l$x, control.l$y, xout=new.samps$time_nr[new.samps$tree_ID %in% c(990, 1099, 1377)])
points(Control.interp, col=4)
dev.off()

pdf("Plots/Soil_moisture_DE970_interpolation.pdf", width = 10, height = 5)
plot(sm$time_nr, sm$DE970, col=2, ylim=c(0,1.1), type="l", main="970 tree", ylab="Soil moisture", xlab="Time")
abline(v=new.samps$time_nr[new.samps$tree_ID == 970], col="grey")
DE970.interp <- approx(sm$time_nr, sm$DE970, xout=new.samps$time_nr[new.samps$tree_ID == 970])
points(DE970.interp, col=4)
dev.off()

pdf("Plots/Soil_moisture_DE8212_interpolation.pdf", width = 10, height = 5)
plot(sm$time_nr, sm$DE8212, col=2, ylim=c(0,1.1), type="l", main="8212 tree", ylab="Soil moisture", xlab="Time")
abline(v=new.samps$time_nr[new.samps$tree_ID == 8212], col="grey")
DE8212.interp <- approx(sm$time_nr, sm$DE8212, xout=new.samps$time_nr[new.samps$tree_ID == 8212])
points(DE8212.interp, col=4)
dev.off()

pdf("Plots/Soil_moisture_DE8266_interpolation.pdf", width = 10, height = 5)
plot(sm$time_nr, sm$DE8266, col=2, ylim=c(0,1.1), type="l", main="8266 tree", ylab="Soil moisture", xlab="Time")
abline(v=new.samps$time_nr[new.samps$tree_ID == 8266], col="grey")
DE8266.interp <- approx(sm$time_nr, sm$DE8266, xout=new.samps$time_nr[new.samps$tree_ID == 8266])
points(DE8266.interp, col=4)
dev.off()


### lowess + interpolation // loess

control <- rbind(as.matrix(sm[!is.na(sm$C970), c("time_nr", "C970")]),
                 as.matrix(sm[!is.na(sm$C8266), c("time_nr", "C8266")]),
                 as.matrix(sm[!is.na(sm$C8212a), c("time_nr", "C8212a")]),
                 as.matrix(sm[!is.na(sm$C8212b), c("time_nr", "C8212b")]))
control <- as.data.frame(control)
colnames(control) <- c("time_nr", "Soil.Moisture")
control <- control[order(control$time_nr), ]
head(control)



### lowess
control.l <- lowess(control[,1], control[,2], f=0.01)

pdf("Plots/Soil_moisture_Control_lowess.pdf", width = 10, height = 5)
plot(control, col=2, ylim=c(0,1.1), pch=20, main="Control trees", ylab="Soil moisture", xlab="Time")
lines(control.l)
dev.off()


### interpolate.and.roll

new.samps.control <- interpolate.and.roll(intrp.x=control.l$x, intrp.y=control.l$y, intrp.points=all.days, table1=new.samps[new.samps$drough.control=="control",], col.values="Soil.Moisture", plot.name="Control")

new.samps.DE970 <- interpolate.and.roll(intrp.x=sm$time_nr, intrp.y=sm$DE970, intrp.points=all.days, table1=new.samps[new.samps$tree_ID == 970,], col.values="Soil.Moisture", plot.name="DE970")

new.samps.DE8212 <- interpolate.and.roll(intrp.x=sm$time_nr, intrp.y=sm$DE8212, intrp.points=all.days, table1=new.samps[new.samps$tree_ID == 8212,], col.values="Soil.Moisture", plot.name="DE8212")

new.samps.DE8266 <- interpolate.and.roll(intrp.x=sm$time_nr, intrp.y=sm$DE8266, intrp.points=all.days, table1=new.samps[new.samps$tree_ID == 8266,], col.values="Soil.Moisture", plot.name="DE8266")

new.samps <- rbind(new.samps.control, new.samps.DE970, new.samps.DE8212, new.samps.DE8266)

names(new.samps)

### loess.and.roll

new.samps.control <- loess.and.roll(intrp.x=control[,1], intrp.y=control[,2], intrp.points=all.days, table1=new.samps[new.samps$drough.control=="control",], col.values="Soil.Moisture", plots.path="PlotsLoess", plot.name="Control", span=3/nrow(control)*6)

new.samps.DE970 <- loess.and.roll(intrp.x=sm$time_nr, intrp.y=sm$DE970, intrp.points=all.days, table1=new.samps[new.samps$tree_ID == 970,], col.values="Soil.Moisture", plots.path="PlotsLoess", plot.name="DE970", span=3/nrow(sm)*8)

new.samps.DE8212 <- loess.and.roll(intrp.x=sm$time_nr, intrp.y=sm$DE8212, intrp.points=all.days, table1=new.samps[new.samps$tree_ID == 8212,], col.values="Soil.Moisture", plots.path="PlotsLoess", plot.name="DE8212", span=3/nrow(sm)*8)

# negative values for sample35
new.samps.DE8212[new.samps.DE8212$sample_num=="sample35", c("Soil.Moisture", "Soil.Moisture14", "Soil.Moisture28")] <- NA

new.samps.DE8266 <- loess.and.roll(intrp.x=sm$time_nr, intrp.y=sm$DE8266 ,intrp.points=all.days, table1=new.samps[new.samps$tree_ID == 8266,], col.values="Soil.Moisture", plots.path="PlotsLoess", plot.name="DE8266", span=3/nrow(sm)*8)

new.samps <- rbind(new.samps.control, new.samps.DE970, new.samps.DE8212, new.samps.DE8266)

names(new.samps)


####################################################
### Water Potential
####################################################

### prepare Mean.Water.Potential var 

wp <- read.table("Data/meteor/Water_potential_from_Inoue.csv", sep=";", header=T)
head(wp)
names(wp)

### make wp values in (0,1) except SD
library(stringr)
wp[,-c(1, which(str_detect(names(wp), "SD")==T))] <- wp[,-c(1, which(str_detect(names(wp), "SD")==T))] + 1

### change time format
wp$time_ch <- strptime(wp$Date, "%y.%m.%d")
wp$time_nr <- as.numeric(wp$time_ch)

#tree.ids <- unique(as.character(new.samps[,"tree_ID"]))
tree.ids <- c("970" , "1099", "8266", "8212", "1377", "990" )
exp.tree.ids <-  c("DE970" , "C1099", "DE8266", "DE8212", "C1377", "C990" )
interp <- list()

### PLOTS of raw WP & sampling points & interpolation in sampling points
for(i in 1:length(tree.ids)){
  # i=5 
  tree.sampl <- paste(exp.tree.ids[i], "sMean", sep="")
  
  pdf(paste("Plots/Water_potential_", exp.tree.ids[i], ".pdf", sep=""), width = 10, height = 5)
  plot(wp$time_nr[!is.na(wp[,tree.sampl])], wp[!is.na(wp[,tree.sampl]), tree.sampl], type="l", col=3, ylim=c(0,1), xlim=c(1.220e+09, 1257289200), main=exp.tree.ids[i], ylab="Mean Water Potential", xlab="Time")
  abline(v=new.samps$time_nr[new.samps$tree_ID==tree.ids[i]], col="grey")
  dev.off()
  
  pdf(paste("Plots/Water_potential_", exp.tree.ids[i], "_interpolation.pdf", sep=""), width = 10, height = 5)
  plot(wp$time_nr[!is.na(wp[,tree.sampl])], wp[!is.na(wp[,tree.sampl]), tree.sampl], type="l", col=3, ylim=c(0,1), main=exp.tree.ids[i], ylab="Mean Water Potential", xlab="Time")
  abline(v=new.samps$time_nr[new.samps$tree_ID==tree.ids[i]], col="grey")
  interp[[tree.ids[i]]] <- approx(wp$time_nr[!is.na(wp[,tree.sampl])], wp[!is.na(wp[,tree.sampl]), tree.sampl], xout=new.samps$time_nr[new.samps$tree_ID==tree.ids[i]])
  points(interp[[tree.ids[i]]], col=4)
  dev.off()
  
}

# for(i in 1:length(tree.ids)){
# new.samps$Mean.Water.Potential[new.samps$tree_ID == tree.ids[i]] <- interp[[tree.ids[i]]]$y
# }

pdf(paste("Plots/Water_potential.pdf", sep=""), width = 10, height = 5)
i=1
tree.sampl <- paste(exp.tree.ids[i], "sMean", sep="")
plot(wp$time_nr[!is.na(wp[,tree.sampl])], wp[!is.na(wp[,tree.sampl]), tree.sampl], type="l", col=i, ylim=c(0,1), xlim=c(min(wp$time_nr), max(wp$time_nr)), main="Water Potential", ylab="Mean Water Potential", xlab="Time")
for(i in 2:length(tree.ids)){
  # i=5
  tree.sampl <- paste(exp.tree.ids[i], "sMean", sep="")
  lines(wp$time_nr[!is.na(wp[,tree.sampl])], wp[!is.na(wp[,tree.sampl]), tree.sampl], col=i)
}
legend("topleft", exp.tree.ids, col=1:length(tree.ids), lty=rep(1, 6) , cex=0.4)
dev.off()


### interpolation and rolled means 

tree.ids <- c("970" , "1099", "8266", "8212", "1377", "990" )
exp.tree.ids <-  c("DE970" , "C1099", "DE8266", "DE8212", "C1377", "C990" )

new.samps.wp <- NULL

### intepolation
for(i in 1:length(tree.ids)){
  # i=1
  tree.sampl <- paste(exp.tree.ids[i], "sMean", sep="")
  new.samps.wp <- rbind(new.samps.wp, new.samps.wp.tmp <- interpolate.and.roll(intrp.x=wp$time_nr[!is.na(wp[,tree.sampl])], intrp.y=wp[!is.na(wp[,tree.sampl]),tree.sampl], intrp.points=all.days, table1=new.samps[new.samps$tree_ID == tree.ids[i],] , col.values="Water.Potential", plot.name=tree.sampl))
  
}

names(new.samps.wp)
new.samps <- new.samps.wp

save(new.samps, file="Data/new_samps_interpolation.RData")


tree.ids <- c("970" , "1099", "8266", "8212", "1377", "990" )
exp.tree.ids <-  c("DE970" , "C1099", "DE8266", "DE8212", "C1377", "C990" )

new.samps.wp <- NULL

### loess
for(i in 1:length(tree.ids)){
  # i=1
  tree.sampl <- paste(exp.tree.ids[i], "sMean", sep="")
  new.samps.wp <- rbind(new.samps.wp, new.samps.wp.tmp <- loess.and.roll(intrp.x=wp$time_nr[!is.na(wp[,tree.sampl])], intrp.y=wp[!is.na(wp[,tree.sampl]),tree.sampl], intrp.points=all.days, table1=new.samps[new.samps$tree_ID == tree.ids[i],], col.values="Water.Potential", plots.path="PlotsLoess", plot.name=tree.sampl, span=0.5))
  
}

names(new.samps.wp)
new.samps <- new.samps.wp

save(new.samps, file="Data/new_samps_loess.RData")


####################################################
### Temperature
####################################################

tempr <- read.table("Data/meteor/Lambir_temperature_data.csv", head=T, sep=";")
head(tempr)
nrow(tempr)

tempr$time_ch <- strptime(tempr$Time, "%d.%m.%y %H:%M")
tempr$time_nr <- as.numeric(tempr$time_ch)
library(stringr)
tempr$time_day <- substring(tempr$Time, 1, 8)
tempr$time_day_ch <- strptime(tempr$time_day, "%d.%m.%y")
tempr$time_day_nr <- as.numeric(tempr$time_day_ch)

tempr$time_hour <- substring(tempr$Time, 10, 11)
tempr$time_hour_nr <- as.numeric(tempr$time_hour)

### calculate average temp per day 
days <- unique(tempr$time_day_ch)
tempr.day <- data.frame(time_day_ch=days, time_day_nr=as.numeric(days) ,temp_avg_day=0)

for(i in tempr.day$time_day_ch){
  #i=days[1]
  tempr.day[tempr.day$time_day_ch==i , "temp_avg_day"] <- mean(tempr$Temp[tempr$time_day_ch==i])
  
}

### some plots of raw temperature

pdf("Plots/Temperature.pdf", width = 10, height = 5)
plot(tempr$time_nr, tempr$Temp, pch=".", main="Temperature", xlab="Time", ylab="Temperature")
abline(v=new.samps$time_nr, col="grey")
dev.off()


pdf("Plots/Temperature_h12.pdf", width = 10, height = 5)
plot(tempr$time_nr[tempr$time_hour_nr==12], tempr$Temp[tempr$time_hour_nr==12], pch=".", main="Temperature at 12:00", xlab="Time", ylab="Temperature")
abline(v=new.samps$time_nr, col="grey")
dev.off()

pdf("Plots/Temperature_avg_day.pdf", width = 10, height = 5)
plot(tempr.day$time_day_nr, tempr.day$temp_avg_day, pch=20, main="Average Temperature per day", xlab="Time", ylab="Temperature")
abline(v=new.samps$time_nr, col="grey")
dev.off()

### plot some rolled means and sums of temperature

library(zoo)
tempr.day$temp_sum14 <- NA
tempr.day$temp_sum14[14:length(tempr.day$temp_sum14)] <- rollsum(tempr.day$temp_avg_day, 14)

pdf("Plots/Temperature_sum14.pdf", width = 10, height = 5)
plot(tempr.day$time_day_nr, tempr.day$temp_sum14, pch=20, main="Temperature Sum over 14 days", xlab="Time", ylab="Temperature")
abline(v=new.samps$time_nr, col="grey")
dev.off()

tempr.day$temp_sum28 <- NA
tempr.day$temp_sum28[28:length(tempr.day$temp_sum28)] <- rollsum(tempr.day$temp_avg_day, 28)

pdf("Plots/Temperature_sum28.pdf", width = 10, height = 5)
plot(tempr.day$time_day_nr, tempr.day$temp_sum28, pch=20, main="Temperature Sum over 28 days", xlab="Time", ylab="Temperature")
abline(v=new.samps$time_nr, col="grey")
dev.off()


tempr.day$temp_avg14 <- NA
tempr.day$temp_avg14[14:length(tempr.day$temp_avg14)] <- rollmean(tempr.day$temp_avg_day, 14)

pdf("Plots/Temperature_avg14.pdf", width = 10, height = 5)
plot(tempr.day$time_day_nr, tempr.day$temp_avg14, pch=20, main="Temperature Average over 14 days", xlab="Time", ylab="Temperature")
abline(v=new.samps$time_nr, col="grey")
dev.off()

tempr.day$temp_avg28 <- NA
tempr.day$temp_avg28[28:length(tempr.day$temp_avg28)] <- rollmean(tempr.day$temp_avg_day, 28)

pdf("Plots/Temperature_avg28.pdf", width = 10, height = 5)
plot(tempr.day$time_day_nr, tempr.day$temp_avg28, pch=20, main="Temperature Average over 28 days", xlab="Time", ylab="Temperature")
abline(v=new.samps$time_nr, col="grey")
dev.off()

head(new.samps)

# new.samps$Temperature_avg14 <- NA
# new.samps$Temperature_avg28 <- NA
# new.samps$Temperature_avg_day <- NA
# 
# for(i in unique(new.samps$time_nr)){
#   
#   new.samps[new.samps$time_nr == i, "Temperature_avg14"] <- tempr.day$temp_avg14[tempr.day$time_day_nr==i]
#   new.samps[new.samps$time_nr == i, "Temperature_avg28"] <- tempr.day$temp_avg28[tempr.day$time_day_nr==i]
#   new.samps[new.samps$time_nr == i, "Temperature_avg_day"] <- tempr.day$temp_avg_day[tempr.day$time_day_nr==i]
# }

### intepolate & roll

new.samps.t <- interpolate.and.roll(intrp.x=tempr.day$time_day_nr, intrp.y=tempr.day$temp_avg_day, intrp.points=all.days, table1=new.samps, col.values="Temperature")

names(new.samps.t)
new.samps <- new.samps.t

save(new.samps, file="Data/new_samps_intepolation.RData")
write.table(new.samps, "Data/new_samps_intepolation.csv", quote=F, sep=";", row.names=F)

### loess

new.samps.t <- loess.and.roll(intrp.x=tempr.day$time_day_nr, intrp.y=tempr.day$temp_avg_day, intrp.points=all.days, table1=new.samps, col.values="Temperature", plots.path="PlotsLoess", span=0.2)

names(new.samps.t)
new.samps <- new.samps.t

save(new.samps, file="Data/new_samps_loess.RData")
write.table(new.samps, "Data/new_samps_loess.csv", quote=F, sep=";", row.names=F)


#####################################################################################################

### GLM // gene expression vs Soil Moisture & Water Potential

#####################################################################################################

setwd("/Users/gosia/DATA/Shimizu_RNA_seq/")

### loess samples
load(file="Data/new_samps_loess.RData")
head(new.samps)
rownames(new.samps) <- new.samps$sample_name

x <- read.table("Data/raw_data-clean.csv", sep=",", header=T, row.names=1)
head(x)
x <- x[, new.samps$sample_name]


AT.id <- read.table("Data/genes_descr_control/best_hit_blast_result_Sl_predicted_exons.txt", sep=",", stringsAsFactors=FALSE)

UP.genes <- read.table("Data/genes_descr_control/mDr_Day10_drought_up_regulated_genes_Harb_etal.csv", sep=";", header=T, stringsAsFactors=FALSE)

DOWN.genes <- read.table("Data/genes_descr_control/mDr_Day10_drought_down_regulated_genes_Harb_etal.csv", sep=";", header=T, stringsAsFactors=FALSE)


MolEcol.genes.c1 <- read.table("Data/genes_descr_control/gene_list_in_clusterI_Mol_Ecol.csv", sep=";", header=T, stringsAsFactors=FALSE)

MolEcol.genes.c2 <- read.table("Data/genes_descr_control/gene_list_in_clusterII_Mol_Ecol.csv", sep=";", header=T, stringsAsFactors=FALSE)


####################################################
###  Soil Moisture
####################################################


### except samples from November and 8266flowered (id: E7_8266_20090416)

x <- x[,!is.na(new.samps$Soil.moisture)]
new.samps <- new.samps[!is.na(new.samps$Soil.moisture), ]

x <- x[,!names(x)=="E7_8266_20090416"]
new.samps <- new.samps[!new.samps$sample_name=="E7_8266_20090416", ]



# ### nls
# param <- nls(expr ~ a + b * sm, data=data.frame(sm=soil.moist , expr=log(expr.back[1,]+1)), start = list(a = 0, b = 0))
# s.param <- summary(param)
# s.param$parameters
# lines(seq(0,1,0.1), s.param$parameters["a","Estimate"] + s.param$parameters["b","Estimate"] * seq(0,1,0.1), col=2)
# ### glm
# param.glm <- glm(expr ~ sm, data=data.frame(sm=soil.moist , expr=log(expr.back[1,]+1)), family = gaussian)
# summary(param.glm)


### edgeR

library(edgeR)
d <- DGEList(x, group=new.samps$tree_ID[match(names(x), new.samps$sample_name)])
d <- calcNormFactors(d)

### make sure a gene is expressed (CPM > 1) in more than 2 samples
d.cpm <- cpm(d, normalized.lib.sizes=TRUE)
d <- d[ rowSums(d.cpm>1) > 2, ]

# design model matrix
design <- model.matrix(~ Soil.moisture, data=new.samps)

# estimate dispersion

d <- estimateGLMCommonDisp(d,design)
d <- estimateGLMTrendedDisp(d,design)
d <- estimateGLMTagwiseDisp(d,design)


# glmFit fits genewise negative binomial glms, all with the same design matrix but possibly different dispersions, offsets and weights

fit <- glmFit(d,design)
fit

lrt <- glmLRT(fit)
lrt


genes <- rownames(d$counts)

## plot TOP genes // expression vs soil moisture
top.sm <- topTags(lrt, n=200)
head(top.sm$table)
genes.top.sm <- rownames(top.sm$table)

AT.id <- read.table("Data/genes_descr_control/best_hit_blast_result_Sl_predicted_exons.txt", sep=",", stringsAsFactors=FALSE)
head(AT.id)

UP.genes <- read.table("Data/genes_descr_control/mDr_Day10_drought_up_regulated_genes_Harb_etal.csv", sep=";", header=T, stringsAsFactors=FALSE)
head(UP.genes)

DOWN.genes <- read.table("Data/genes_descr_control/mDr_Day10_drought_down_regulated_genes_Harb_etal.csv", sep=";", header=T, stringsAsFactors=FALSE)


pdf("Plots/Top_fitting_Soil_Moisture_counts.pdf", width = 9, height =6)
par(mfrow=c(2,3))
for(i in genes.top.sm){
  # i=genes.top.sm[1]
  
  AT <- AT.id[AT.id[,1] == i, 2]
  
  if.UP <- AT %in% UP.genes$Gene
  if.DOWN <- AT %in% DOWN.genes$Gene
  
  coeffs <- round(fit$coefficients[i,], 2)
  
  plot(fit$design[,2], fit$counts[i,], main=paste(i, "/", AT, "\n", "Coefffs:",coeffs[1], ",",coeffs[2], "\n FDR:", top.sm$table[i, "FDR"], "\n UP:", if.UP, "DOWN:", if.DOWN), xlab="Soil moisture", ylab="Counts", cex.main=0.7)
  points(fit$design[,2], fit$fitted.values[i,], col=2)
  
}
dev.off()


d.fit <- DGEList(fit$fitted.values, group=new.samps$tree_ID[match(names(x), new.samps$sample_name)])
d.fit <- calcNormFactors(d.fit)

d.cpm <- cpm(d, normalized.lib.sizes=TRUE)
d.fit.cpm <- cpm(d.fit, normalized.lib.sizes=TRUE)




pdf("Plots/Top_fitting_Soil_Moisture_cpm.pdf", width = 9, height =6)
par(mfrow=c(2,3))
for(i in genes.top.sm){
  # i=genes.top.sm[1]
  
  AT <- AT.id[AT.id[,1] == i, 2]
  
  if.UP <- AT %in% UP.genes$Gene
  if.DOWN <- AT %in% DOWN.genes$Gene
  
  coeffs <- round(fit$coefficients[i,], 2)
  
  plot(fit$design[,2], d.cpm[i,], main=paste(i, "/", AT, "\n", "Coefffs:",coeffs[1], ",",coeffs[2], "\n FDR:", top.sm$table[i, "FDR"], "\n UP:", if.UP, "DOWN:", if.DOWN), xlab="Soil moisture", ylab="Counts", cex.main=0.7, pch=ifelse(new.samps$drough.control=="drought", 5, 1))
  points(fit$design[,2], d.fit.cpm[i,], col=2)
  
}
dev.off()


### plot expression vs time // drought experiment & control


pdf("Plots/Top_Time_Soil_Moisture_cpm.pdf", width = 9, height =6)
par(mfrow=c(2,1))
plot(new.samps$time_nr, new.samps$Soil.moisture, pch=ifelse(new.samps$drough.control=="drought", 5, 1), col=ifelse(new.samps$drough.control=="drought", 2, 4))

plot(new.samps$time_nr, d.cpm[i,], pch=ifelse(new.samps$drough.control=="drought", 5, 1), col=ifelse(new.samps$drough.control=="drought", 2, 4))
dev.off()


plot(new.samps$time_nr[new.samps$drough.control=="drought"], d.cpm[i,new.samps$drough.control=="drought"], pch=5)

plot(new.samps$time_nr[new.samps$drough.control=="control"], d.cpm[i,new.samps$drough.control=="control"], pch=1)



# ### find out the AT id of RD29A, RD20B, COR15A drought response genes
# library(org.At.tair.db)
# at.x <- as.list(org.At.tair.db::org.At.tairSYMBOL)
# lookup <- data.frame( id=rep(names(at.x), sapply(at.x, length)), alias=unlist(at.x,use.names=FALSE) )
# 
# lookup[which(lookup$alias=="RD29A"),]
# which(UP.genes$Gene=="AT5G52310")
# which(DOWN.genes$Gene=="AT5G52310")
# which(AT.id$V2=="AT5G52310")
# 
# lookup[which(lookup$alias=="COR15A"),]
# which(UP.genes$Gene=="AT2G42540")
# which(DOWN.genes$Gene=="AT2G42540")
# which(AT.id$V2=="AT2G42540")

####################################################
###  Mean Water Potential
####################################################


### except samples for which interpolation was not possible and 8266flowered (id: E7_8266_20090416)

x <- x[,!is.na(new.samps$Mean.Water.Potential)]
new.samps <- new.samps[!is.na(new.samps$Mean.Water.Potential), ]

x <- x[,!names(x)=="E7_8266_20090416"]
new.samps <- new.samps[!new.samps$sample_name=="E7_8266_20090416", ]

### edgeR

library(edgeR)
d <- DGEList(x, group=new.samps$tree_ID[match(names(x), new.samps$sample_name)])
d <- calcNormFactors(d)

### make sure a gene is expressed (CPM > 1) in more than 2 samples
d.cpm <- cpm(d, normalized.lib.sizes=TRUE)
d <- d[ rowSums(d.cpm>1) > 2, ]

# design model matrix
design <- model.matrix(~ Mean.Water.Potential, data=new.samps)

# estimate dispersion

d <- estimateGLMCommonDisp(d,design)
d <- estimateGLMTrendedDisp(d,design)
d <- estimateGLMTagwiseDisp(d,design)


# glmFit fits genewise negative binomial glms, all with the same design matrix but possibly different dispersions, offsets and weights

fit <- glmFit(d,design)
fit$deviance

lrt <- glmLRT(fit)
lrt


genes <- rownames(d$counts)

## plot TOP genes // expression vs soil moisture
top.wp <- topTags(lrt, n=200)
head(top.wp$table)
genes.top.wp <- rownames(top.wp$table)

d.fit <- DGEList(fit$fitted.values, group=new.samps$tree_ID[match(names(x), new.samps$sample_name)])
d.fit <- calcNormFactors(d.fit)

d.cpm <- cpm(d, normalized.lib.sizes=TRUE)
d.fit.cpm <- cpm(d.fit, normalized.lib.sizes=TRUE)


pdf("Plots/Top_fitting_Mean_Water_Potential_cpm.pdf", width = 9, height =6)
par(mfrow=c(2,3))
for(i in genes.top.wp){
  # i=genes.top.wp[1]  
  AT <- AT.id[AT.id[,1] == i, 2]  
  if.UP <- AT %in% UP.genes$Gene
  if.DOWN <- AT %in% DOWN.genes$Gene  
  coeffs <- round(fit$coefficients[i,], 2)
  
  plot(fit$design[,2], d.cpm[i,], main=paste(i, "/", AT, "\n", "Coefffs:",coeffs[1], ",",coeffs[2], "\n FDR:", top.sm$table[i, "FDR"], "\n UP:", if.UP, "DOWN:", if.DOWN), xlab="Mean Water Potential", ylab="Counts", cex.main=0.7, pch=ifelse(new.samps$drough.control=="drought", 5, 1))
  points(fit$design[,2], d.fit.cpm[i,], col=3) 
}
dev.off()

####


intersect(genes.top.wp, genes.top.sm)


####################################################
###  Fitting model automaticly
####################################################

### run.edgeR()

run.edgeR <- function(x, new.samps, model.formula, varialbs, elim.samps, n.top=200, plots.path="PlotsTest", plot.name="", AT.id, UP.genes, DOWN.genes, design.col.plot = 2, fit.clr = 3){ 
  
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
  
  # design model matrix
  design <- model.matrix(model.formula, data=new.samps)
  
  # estimate dispersion  
  d <- estimateGLMCommonDisp(d,design)
  d <- estimateGLMTrendedDisp(d,design)
  d <- estimateGLMTagwiseDisp(d,design)
  
  # glmFit
  fit <- glmFit(d,design)
  lrt <- glmLRT(fit)
  
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
  
  dir.create(plots.path, recursive=T)
  
  pdf( paste(plots.path, "/Top_fitting_", model.char, "_cpm", plot.name,".pdf", sep=""), width = 9, height =6)
  par(mfrow=c(2,3))
  
  for(i in top.genes){
    # i="GID006843_3010436"  
    Top.fitting.list[[i]] <- data.frame(fit$design[, design.col.plot], log(d.cpm[i,]), log(d.fit.cpm[i,]))
    colnames(Top.fitting.list[[i]]) <- c(colnames(fit$design)[design.col.plot], "Log.Counts.in.CPM.RAW", "Log.Counts.in.CPM.FIT")
    rownames(Top.fitting.list[[i]]) <- row.names(fit$samples)
    coeffs <- round(fit$coefficients[i,], 2)
    
    AT <- AT.id[AT.id[,1] == i, 2]  
    if(length(AT)!=0){
      if.UP <- AT %in% UP.genes$Gene
      if(if.UP)
        UP.coeffs[i,] <- coeffs
      if.DOWN <- AT %in% DOWN.genes$Gene  
      if(if.DOWN)
        DOWN.coeffs[i,] <- coeffs
    }
    
    plot(fit$design[, design.col.plot], log(d.cpm[i,]), main=paste(i, "/", AT, "\n", "Coefffs:",paste(coeffs, collapse=", "), "\n FDR:", top.tags$table[i, "FDR"], "\n UP:", if.UP, "DOWN:", if.DOWN), xlab = colnames(fit$design)[design.col.plot], ylab="Log(Counts in CPM)", cex.main=0.7, pch=ifelse(new.samps$drough.control=="drought", 5, 1))
    points(fit$design[, design.col.plot], log(d.fit.cpm[i,]), col=fit.clr)
    
  }
  
  dev.off()
  
  #save(Top.fitting.list, file=paste(plots.path, "/Top_fitting_", model.char, "_cpm", plot.name,".RData", sep=""))
  
  UP.coeffs <- UP.coeffs[complete.cases(UP.coeffs), ]  
  DOWN.coeffs <- DOWN.coeffs[complete.cases(DOWN.coeffs), ]
  
  invisible(return(list(Top.fitting.list=Top.fitting.list, UP.coeffs=UP.coeffs, DOWN.coeffs=DOWN.coeffs)))
  
}


### RUN run.egdeR()

setwd("/Users/gosia/DATA/Shimizu_RNA_seq/")

load(file="Data/new_samps.RData")
names(new.samps)

x <- read.table("Data/raw_data-clean.csv", sep=",", header=T, row.names=1)

AT.id <- read.table("Data/genes_descr_control/best_hit_blast_result_Sl_predicted_exons.txt", sep=",", stringsAsFactors=FALSE)

UP.genes <- read.table("Data/genes_descr_control/mDr_Day10_drought_up_regulated_genes_Harb_etal.csv", sep=";", header=T, stringsAsFactors=FALSE)

DOWN.genes <- read.table("Data/genes_descr_control/mDr_Day10_drought_down_regulated_genes_Harb_etal.csv", sep=";", header=T, stringsAsFactors=FALSE)

### 

elim.samps <- c("E7_8266_20090416") # flowered sample

### Soil.moisture

model.formula <- ~ Soil.moisture
varialbs <- c("Soil.moisture")

Out.SM <- run.edgeR(x, new.samps, model.formula, varialbs, elim.samps, n.top=200, plots.path="PlotsRUNedgeR", plot.name="", AT.id, UP.genes, DOWN.genes, design.col.plot = 2, fit.clr = 2)

### Water potential

class(model.formula <- ~ Mean.Water.Potential + drough.control)
#class(as.formula("~ Mean.Water.Potential + drough.control"))
varialbs <- c("Mean.Water.Potential", "drough.control")

Out.MWP <- run.edgeR(x, new.samps, model.formula, varialbs, elim.samps, n.top=200, plots.path="PlotsTest", plot.name="", AT.id, UP.genes, DOWN.genes, design.col.plot = 2, fit.clr = 3)

Out.MWP$UP.coeffs


#intersect

### Temperature

model.formula <- ~ Temperature_avg_day
varialbs <- c("Temperature_avg_day")

Out.T <- run.edgeR(x, new.samps, model.formula, varialbs, elim.samps, n.top=200, plots.path="PlotsRUNedgeR", plot.name="", AT.id, UP.genes, DOWN.genes, design.col.plot = 2, fit.clr = 4)



model.formula <- ~ Temperature_avg14
varialbs <- c("Temperature_avg14")

Out.T <- run.edgeR(x, new.samps, model.formula, varialbs, elim.samps, n.top=200, plots.path="PlotsRUNedgeR", plot.name="", AT.id, UP.genes, DOWN.genes, design.col.plot = 2, fit.clr = 4)


model.formula <- ~ Temperature_avg28
varialbs <- c("Temperature_avg28")

Out.T <- run.edgeR(x, new.samps, model.formula, varialbs, elim.samps, n.top=200, plots.path="PlotsRUNedgeR", plot.name="", AT.id, UP.genes, DOWN.genes, design.col.plot = 2, fit.clr = 5)

Out.T$Top.fitting.list
### check outliers

load("PlotsTest/Top_fitting_Mean_Water_Potential+drough_control_cpm.RData")


rownames(Top.fitting.list$GID013606_3470013)[which(Top.fitting.list$GID013606_3470013["Log.Counts.in.CPM"] > 1)]

rownames(Top.fitting.list$GID049681_3534549)[which(Top.fitting.list$GID049681_3534549["Log.Counts.in.CPM"] > 2)]

rownames(Top.fitting.list$GID064624_3540005)[which(Top.fitting.list$GID064624_3540005["Log.Counts.in.CPM"] > 4)]

rownames(Top.fitting.list$GID047654_3533666)[which(Top.fitting.list$GID047654_3533666["Log.Counts.in.CPM"] > 1)]














