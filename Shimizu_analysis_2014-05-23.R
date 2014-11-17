#####################################################################################################
# + preparation of variables describing samples
# + crude analysis - MDS plots
# + plots of expression for flowering genes
# + PLOTS OF SELECTED GENES
# + GLM - models with one covariate
# + clustering
# + GO analysis with topGO
# + Fisher's exact tests



# + RUV analysis 
# + Fisher exact test for MolEcol clusters: Sh.l <-> Sh.b 

# + CAMERA Gene set testing

#####################################################################################################


#####################################################################################################

### preparation of variables describing samples 

#####################################################################################################

setwd("/home/Shared/data/seq/Shimizu_RNA_seq/")

####################################################
### samps
####################################################

dir.create("Samples_out")

### parse sample information
samps <- read.table("Data/sample_list.csv", sep=";", stringsAsFactors=FALSE, header=TRUE)
head(samps)
names(samps)

### CORRECT THE DATE!!!
samps$year.month.day.old <- samps$year.month.day
samps$year.month.day <- paste(substr(samps$sample_name, nchar(samps$sample_name)-1, nchar(samps$sample_name)), substr(samps$sample_name, nchar(samps$sample_name)-3, nchar(samps$sample_name)-2), substr(samps$sample_name, nchar(samps$sample_name)-5, nchar(samps$sample_name)-4), sep=".")

samps[,c("sample_name","year.month.day", "year.month.day.old")]


# columns 8-16 come from Lambir_meteorological_data_...xls
# no genes variables 
new.samps <- samps[,c("sample_num", "sample_name","sample_ID","tree_ID","drough.control","developmental_stage","year.month.day")]

# define param for ploting: colors, legend
rownames(new.samps) <- new.samps$sample_name
new.samps$tree_ID <- as.factor(new.samps$tree_ID)
new.samps$tree_col <- new.samps$tree_ID
levels(new.samps$tree_col) <- c("orange", "cyan3", "green3", "blue", "red", "magenta3")
new.samps$tree_col <- as.character(new.samps$tree_col)
new.samps$tree_col[new.samps$sample_name=="E7_8266_20090416"] <- "darkmagenta"
new.samps$tree_legend <- paste(new.samps$tree_ID, new.samps$drough.control, new.samps$developmental_stage, sep="-")
new.samps$short.name <- paste(new.samps$tree_ID, paste(substr(new.samps$year.month.day, 7,8), substr(new.samps$year.month.day, 4,5), substr(new.samps$year.month.day, 1,2), sep="."), sep="_")
head(new.samps)
### time format change 
new.samps$time_ch <- strptime(new.samps$year.month.day, "%d.%m.%y")
new.samps$time_nr <- as.numeric(new.samps$time_ch)
new.samps$time_ch[which.min(new.samps$time_nr)]
new.samps$time_ch[which.max(new.samps$time_nr)]
new.samps$flowered <- new.samps$tree_ID=="8266"
write.table(new.samps,"Samples_out/samps.xls", quote=FALSE, sep="\t", row.names=FALSE)

# object to plot legends
trees.order <- data.frame(legend=c("990-control-leaf_bud", "1099-control-leaf_bud", "1377-control-leaf_bud", "970-drought-leaf_bud", "8212-drought-leaf_bud", "8266-drought-leaf_bud", "8266-drought-flower_bud"), color=c("cyan3", "green3", "blue", "orange", "red", "magenta3", "darkmagenta"), pch=c(16, 16, 16, 18, 18, 18, 18), cex=c(1, 1, 1, 1, 1, 1, 2) , tree_ID = c("990", "1099", "1377", "970", "8212", "8266", ""), condition=c("C", "C","C", "DE", "DE", "DE", "DE"), stringsAsFactors=F) 


### list of all unique days in 2008 and 2009
ad <- read.table("Data/Unique_days.csv", sep=";")
all.days <- data.frame(days.ch = ad[,], days.nr = as.numeric(strptime(ad[,], "%d.%m.%y")))
head(all.days)
library(stringr)
month.days <- all.days[str_sub(all.days[,1], 1, 2)=="01",]

ad.short <- read.table("Data/Unique_days_short.csv", sep=";")
all.days.short <- data.frame(days.ch = ad.short[,], days.nr = as.numeric(strptime(ad.short[,], "%d.%m.%y")))
head(all.days.short)
library(stringr)
month.days.short <- all.days.short[str_sub(all.days.short[,1], 1, 2)=="01",]


####################################################
### Soild Moisture
####################################################

### get soil moisture data
sm <- read.table("Data/meteor/Soil_Moisture.csv", sep=";", stringsAsFactors=FALSE, header=TRUE)
head(sm)
sm <- sm[1:271,]

### time format change
sm$Year <- ifelse(sm$DayUniq <= 366, 2008, 2009)
sm$year.month.day <- paste(sm$Day,sm$Month, sm$Year, sep=".")
sm$time_ch <- strptime(sm$year.month.day, "%d.%m.%Y")
sm$time_nr <- as.numeric(sm$time_ch)


### plots of raw data
dir.create(path="Plots_Raw", showWarnings=FALSE, recursive=TRUE)

pdf("Plots_Raw/Soil_moisture.pdf", width = 10, height = 5)
plot(sm$time_nr, sm$DE970, type="l", col=trees.order$color[trees.order$tree_ID=="970"], ylim=c(-0.1,1.1), main="Soil Moisture", ylab="Soil moisture", xlab="Time", xaxt = "n", lwd=1.5)
axis(side=1, at=month.days[,2], labels=month.days[,1])
lines(sm$time_nr, sm$C970, col=trees.order$color[trees.order$tree_ID=="970"], lty=2, lwd=1.5)
lines(sm$time_nr, sm$DE8266, col=trees.order$color[trees.order$tree_ID=="8266"], lwd=1.5)
lines(sm$time_nr, sm$C8266, col=trees.order$color[trees.order$tree_ID=="8266"], lty=2, lwd=1.5)
lines(sm$time_nr, sm$DE8212, col=trees.order$color[trees.order$tree_ID=="8212"], lwd=1.5)
lines(sm$time_nr, sm$C8212a, col=trees.order$color[trees.order$tree_ID=="8212"], lty=2, lwd=1.5)
lines(sm$time_nr, sm$C8212b, col=trees.order$color[trees.order$tree_ID=="8212"], lty=3, lwd=1.5)
legend("bottomleft", c("DE970", "C970","DE8266","C8266" ,"DE8212","C8212a" ,"C8212b"), col=c(trees.order$color[trees.order$tree_ID=="970"], trees.order$color[trees.order$tree_ID=="970"], trees.order$color[trees.order$tree_ID=="8266"], trees.order$color[trees.order$tree_ID=="8266"], trees.order$color[trees.order$tree_ID=="8212"], trees.order$color[trees.order$tree_ID=="8212"], trees.order$color[trees.order$tree_ID=="8212"]), lty=c(1,2,1,2,1,2,3), cex=0.5, lwd=1)
dev.off()


pdf("Plots_Raw/Soil_moisture_Control.pdf", width = 10, height = 5)
plot(sm$time_nr, sm$C970, type="l", col=trees.order$color[trees.order$tree_ID=="970"], ylim=c(0,1.1), xlim=c(min(all.days.short[,2]), max(all.days.short[,2])), main="Control trees", ylab="Soil moisture", xlab="Time",  xaxt = "n")
axis(side=1, at=month.days[,2], labels=month.days[,1])
lines(sm$time_nr, sm$C8266, col=trees.order$color[trees.order$tree_ID=="8266"])
lines(sm$time_nr, sm$C8212a, col=trees.order$color[trees.order$tree_ID=="8212"])
lines(sm$time_nr, sm$C8212b, col=trees.order$color[trees.order$tree_ID=="8212"])
abline(v=new.samps$time_nr[new.samps$tree_ID %in% c(990, 1099, 1377)], col="grey")
dev.off()

pdf("Plots_Raw/Soil_moisture_DE970.pdf", width = 10, height = 5)
plot(sm$time_nr, sm$DE970, type="l", col=trees.order$color[trees.order$tree_ID=="970"], ylim=c(0,1.1), xlim=c(min(all.days.short[,2]), max(all.days.short[,2])), main="DE970", ylab="Soil moisture", xlab="Time",  xaxt = "n")
axis(side=1, at=month.days[,2], labels=month.days[,1])
abline(v=new.samps$time_nr[new.samps$tree_ID==970], col="grey")
dev.off()

pdf("Plots_Raw/Soil_moisture_DE8266.pdf", width = 10, height = 5)
plot(sm$time_nr, sm$DE8266, type="l", col=trees.order$color[trees.order$tree_ID=="8266"], ylim=c(0,1.1), xlim=c(min(all.days.short[,2]), max(all.days.short[,2])), main="DE8266", ylab="Soil moisture", xlab="Time",  xaxt = "n")
axis(side=1, at=month.days[,2], labels=month.days[,1])
abline(v=new.samps$time_nr[new.samps$tree_ID==8266], col="grey")
dev.off()

pdf("Plots_Raw/Soil_moisture_DE8212.pdf", width = 10, height = 5)
plot(sm$time_nr, sm$DE8212, type="l", col=trees.order$color[trees.order$tree_ID=="8212"], ylim=c(0,1.1), xlim=c(min(all.days.short[,2]), max(all.days.short[,2])), main="DE8212", ylab="Soil moisture", xlab="Time",  xaxt = "n")
axis(side=1, at=month.days[,2], labels=month.days[,1])
abline(v=new.samps$time_nr[new.samps$tree_ID==8212], col="grey")
dev.off()


# merge control samples
control <- rbind(as.matrix(sm[!is.na(sm$C970), c("time_nr", "C970")]),
                 as.matrix(sm[!is.na(sm$C8266), c("time_nr", "C8266")]),
                 as.matrix(sm[!is.na(sm$C8212a), c("time_nr", "C8212a")]),
                 as.matrix(sm[!is.na(sm$C8212b), c("time_nr", "C8212b")]))
control <- as.data.frame(control)
colnames(control) <- c("time_nr", "Soil.Moisture")
control <- control[order(control$time_nr), ]



### lowess for control
dir.create("Plots_Interpolate")
control.l <- lowess(control[,1], control[,2], f=0.02)

pdf("Plots_Interpolate/Soil_moisture_Control_lowess.pdf", width = 10, height = 5)
plot(control, col="grey", ylim=c(0,1.1), pch=20, main="Control trees", ylab="Soil moisture", xlab="Time",  xaxt = "n")
axis(side=1, at=month.days[,2], labels=month.days[,1])
lines(control.l, lwd=2)
dev.off()


# Interpolate and roll 

source("/home/gosia/R/R_Shimizu_RNA_seq/Interpolate_and_roll.R")

new.samps.control <- interpolate.and.roll(intrp.x=control.l$x, intrp.y=control.l$y, intrp.points=all.days[,2], table1=new.samps[new.samps$drough.control=="control",], roll.days=c(14, 28), col.values="Soil.Moisture", plot.name="Control", color=trees.order$color[trees.order$tree_ID=="1377"])

new.samps.DE970 <- interpolate.and.roll(intrp.x=sm$time_nr, intrp.y=sm$DE970, intrp.points=all.days[,2], table1=new.samps[new.samps$tree_ID == 970,], roll.days=c(14, 28), col.values="Soil.Moisture", plot.name="DE970", color=trees.order$color[trees.order$tree_ID=="970"])

new.samps.DE8212 <- interpolate.and.roll(intrp.x=sm$time_nr, intrp.y=sm$DE8212, intrp.points=all.days[,2], table1=new.samps[new.samps$tree_ID == 8212,], roll.days=c(14, 28), col.values="Soil.Moisture", plot.name="DE8212", color=trees.order$color[trees.order$tree_ID=="8212"])

new.samps.DE8266 <- interpolate.and.roll(intrp.x=sm$time_nr, intrp.y=sm$DE8266, intrp.points=all.days[,2], table1=new.samps[new.samps$tree_ID == 8266,], roll.days=c(14, 28), col.values="Soil.Moisture", plot.name="DE8266",  color=trees.order$color[trees.order$tree_ID=="8266"])

new.samps.sm <- rbind(new.samps.control, new.samps.DE970, new.samps.DE8212, new.samps.DE8266)

new.samps <- new.samps.sm

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


tree.ids <- trees.order$tree_ID[1:6]
exp.tree.ids <- paste0(trees.order$condition[1:6], trees.order$tree_ID[1:6])

### PLOTS of raw WP & sampling points 
for(i in 1:length(tree.ids)){
  # i=5 
  tree.sampl <- paste(exp.tree.ids[i], "sMean", sep="")
  
  pdf(paste("Plots_Raw/Water_potential_", exp.tree.ids[i], ".pdf", sep=""), width = 10, height = 5)
  plot(wp$time_nr[!is.na(wp[,tree.sampl])], wp[!is.na(wp[,tree.sampl]), tree.sampl], type="l", col=trees.order$color[trees.order$tree_ID==tree.ids[i]], ylim=c(0.2,0.8), xlim=c(min(all.days.short[,2]), max(all.days.short[,2])), main=exp.tree.ids[i], ylab="Mean Water Potential", xlab="Time", xaxt = "n")
  axis(side=1, at=month.days[,2], labels=month.days[,1])
  abline(v=new.samps$time_nr[new.samps$tree_ID==tree.ids[i]], col="grey")
  dev.off()
  
  
}


pdf(paste("Plots_Raw/Water_potential.pdf", sep=""), width = 10, height = 5)
i=1
tree.sampl <- paste(exp.tree.ids[i], "sMean", sep="")
plot(wp$time_nr[!is.na(wp[,tree.sampl])], wp[!is.na(wp[,tree.sampl]), tree.sampl], type="l", col=trees.order$color[trees.order$tree_ID==tree.ids[i]], ylim=c(0.2,0.8), xlim=c(min(wp$time_nr), max(wp$time_nr)), main="Water Potential", ylab="Mean Water Potential", xlab="Time", xaxt = "n")
axis(side=1, at=month.days[,2], labels=month.days[,1])
for(i in 2:length(tree.ids)){
  # i=5
  tree.sampl <- paste(exp.tree.ids[i], "sMean", sep="")
  lines(wp$time_nr[!is.na(wp[,tree.sampl])], wp[!is.na(wp[,tree.sampl]), tree.sampl], col=trees.order$color[trees.order$tree_ID==tree.ids[i]])
}
legend("bottomleft", trees.order$legend[1:6], col=trees.order$color[1:6], lty=rep(1, 6) , cex=0.6)
dev.off()




### add avg values to the May samples

wp <- wp[order(wp$time_nr), ]
new.samps <- new.samps[order(new.samps$time_nr), ]

tail(wp)

wp.temp <- as.data.frame(matrix(NA, 2, ncol(wp)))
colnames(wp.temp) <- colnames(wp)

wp.temp$time_nr <- c(1241992800, 1242079200)
wp.temp$time_ch <- strptime(c("2009-05-11", "2009-05-12"), "%Y-%m-%d")

May11.trees <-  c("DE970" , "C1099", "C1377", "C990")
May12.trees <-  c("DE8266")

for(i in 1:length(May11.trees)){
  # i=1
  tree.sampl <- paste(May11.trees[i], "sMean", sep="")
  wp.temp[wp.temp$time_nr == 1241992800, tree.sampl] <- mean(tail(  na.omit(wp[, tree.sampl]), 2))  
}

for(i in 1:length(May12.trees)){
  # i=1
  tree.sampl <- paste(May12.trees[i], "sMean", sep="")
  wp.temp[wp.temp$time_nr == 1242079200, tree.sampl] <- mean(tail(  na.omit(wp[, tree.sampl]), 2))  
}

wp <- rbind(wp, wp.temp)


### interpolation and rolled means 
new.samps.wp <- NULL

### intepolation
for(i in 1:length(tree.ids)){
  # i=1
  tree.sampl <- paste(exp.tree.ids[i], "sMean", sep="")
  new.samps.wp <- rbind(new.samps.wp, new.samps.wp.tmp <- interpolate.and.roll(intrp.x=wp$time_nr[!is.na(wp[,tree.sampl])], intrp.y=wp[!is.na(wp[,tree.sampl]),tree.sampl], intrp.points=all.days[,2], table1=new.samps[new.samps$tree_ID == tree.ids[i],] , roll.days=c(14, 28), col.values="Water.Potential", plot.name=tree.sampl))
  
}

names(new.samps.wp)

new.samps <- new.samps.wp


####################################################
### Temperature
####################################################

tempr <- read.table("Data/meteor/Lambir_temperature_data.csv", head=T, sep=";")
head(tempr)

library(stringr)

tempr$time_ch <- strptime(substring(tempr$Time, 1, 8), "%d.%m.%y")
tempr$time_nr <- as.numeric(tempr$time_ch)

tempr$time_hour <- substring(tempr$Time, 10, 11)
tempr$time_hour_nr <- as.numeric(tempr$time_hour)

pdf("Plots_Raw/Temperature_raw.pdf", width = 10, height = 5)
plot(tempr$time_nr, tempr$Temp, pch=".", main="Temperature", xlab="Time", ylab="Temperature")
#abline(v=new.samps$time_nr, col="grey")
dev.off()

pdf("Plots_Raw/Temperature_over_day.pdf", width = 10, height = 5)
boxplot(Temp ~ time_hour, data=tempr, xlab="Hour", ylab="Temperature")
dev.off()


### calculate average temp per day 
days <- unique(tempr$time_ch)
tempr.day <- data.frame(time_ch=days, time_nr=as.numeric(days) , TempAvg=0, TempMax=0, TempMin=0)

for(i in 1:nrow(tempr.day)){
  # i=1
  tempr.day[i , "TempAvg"] <- mean(tempr$Temp[tempr$time_ch==days[i]])
  tempr.day[i , "TempMax"] <- max(tempr$Temp[tempr$time_ch==days[i] & tempr$time_hour_nr %in% 10:18])
  tempr.day[i , "TempMin"] <- min(tempr$Temp[tempr$time_ch==days[i] & tempr$time_hour_nr %in% 10:18])
}



### lowess for control
dir.create("Plots_Interpolate")
tempr.l <- lowess(tempr.day[,"time_nr"], tempr.day[,"TempAvg"], f=0.0.5)

pdf("Plots_Interpolate/Temperature_lowess.pdf", width = 10, height = 5)
plot(tempr.day[,"time_nr"], tempr.day[,"TempAvg"], col="grey", pch=20, main="Control trees", ylab="Soil moisture", xlab="Time",  xaxt = "n")
axis(side=1, at=month.days[,2], labels=month.days[,1])
lines(tempr.l, lwd=2)
dev.off()




### some plots of avg, max, min temperature
pdf("Plots_Raw/Temperature_Avg.pdf", width = 10, height = 5)
plot(tempr.day$time_nr, tempr.day$TempAvg, pch=20, main="Average Temperature", xlab="Time", ylab="Temperature",xlim=c(min(all.days.short[,2]), max(all.days.short[,2])),  xaxt = "n")
axis(side=1, at=month.days[,2], labels=month.days[,1])
abline(v=new.samps$time_nr, col="grey")
dev.off()

pdf("Plots_Raw/Temperature_Max.pdf", width = 10, height = 5)
plot(tempr.day$time_nr, tempr.day$TempMax, pch=20, main="Max Temperature between 10 and 18", xlab="Time", ylab="Temperature",xlim=c(min(all.days.short[,2]), max(all.days.short[,2])),  xaxt = "n")
axis(side=1, at=month.days[,2], labels=month.days[,1])
abline(v=new.samps$time_nr, col="grey")
dev.off()

pdf("Plots_Raw/Temperature_Min.pdf", width = 10, height = 5)
plot(tempr.day$time_nr, tempr.day$TempMin, pch=20, main="Min Temperature between 10 and 18", xlab="Time", ylab="Temperature",xlim=c(min(all.days.short[,2]), max(all.days.short[,2])),  xaxt = "n")
axis(side=1, at=month.days[,2], labels=month.days[,1])
abline(v=new.samps$time_nr, col="grey")
dev.off()


### rolled means of temperature

library(zoo)

roll.days = c(14, 28)
variables = c("TempAvg", "TempMax", "TempMin")

for(v in variables){
  # v="TempAvg"
  new.samps <- match.func(DF1=new.samps, col.match1="time_nr", DF2=tempr.day, col.match2="time_nr", col.values2=v)
  
  for(r in roll.days){
    # r=14
    name <- paste(v, r, sep="")
    tempr.day[, name] <- NA
    tempr.day[r:length(tempr.day[, name]), name] <- rollmean(tempr.day[,v], r)
    
    new.samps <- match.func(DF1=new.samps, col.match1="time_nr", DF2=tempr.day, col.match2="time_nr", col.values2=name)
  }
}

head(new.samps)

dir.create(path="Samples_out", showWarnings=F, recursive=T)

save(new.samps, file="Samples_out/new_samps_interpolation.RData")

write.table(new.samps, "Samples_out/new_samps_interpolation.xls", sep="\t",  row.names = F)





#####################################################################################################

### Load data 

#####################################################################################################

setwd("/home/Shared/data/seq/Shimizu_RNA_seq/")

### samples & factors
new.samps <- read.table("Samples_out/new_samps_interpolation.xls", header=TRUE, sep="\t", stringsAsFactors=FALSE)

rownames(new.samps) <- new.samps$sample_name

# colors reprezening WP level 

colors.wp <- new.samps[, c("sample_name", "Water.Potential")]
colors.wp$Water.Potential <- 100*round(colors.wp$Water.Potential, digits=2)
colors.palette <- as.character( colorRampPalette(c("red","blue"))(diff(range(na.omit(colors.wp$Water.Potential)))+1))
colors.wp$colors.wp <- colors.palette[colors.wp$Water.Potential - min(na.omit(colors.wp$Water.Potential))+1]
colors.wp$colors.wp[is.na(colors.wp$colors.wp)] <- "#838B8B" 
rownames(colors.wp) <- colors.wp$sample_name
colors.wp <- colors.wp[,"colors.wp", drop=F]

new.samps <- merge(new.samps, colors.wp, by=0, all.x=T)
rownames(new.samps) <- new.samps$sample_name


### counts
x <- read.table("Data/raw_data-clean.csv", sep=",", header=T, row.names=1)
x <- x[, new.samps$sample_name]


### IMPORTANT - ORDER
samps.order <- row.names(new.samps[order(new.samps$drough.control, new.samps$tree_ID, new.samps$time_nr), ])
x <- x[, samps.order]
new.samps <- new.samps[samps.order,]


# normalise between samples
library(edgeR)
d.org <- DGEList(x, group=new.samps$tree_ID)
d.org <- calcNormFactors(d.org)

d.cpm <- cpm(d.org, normalized.lib.sizes=TRUE)
d.cpm.l <- log(d.cpm + min(d.cpm[d.cpm != 0]))



flowered.samps <- c("E7_8266_20090416") # flowered sample
control.samps <- rownames(new.samps[new.samps$drough.control=="control",])


# object to plot legends
trees.order <- data.frame(legend=c("990-control-leaf_bud", "1099-control-leaf_bud", "1377-control-leaf_bud", "970-drought-leaf_bud", "8212-drought-leaf_bud", "8266-drought-leaf_bud", "8266-drought-flower_bud"), color=c("cyan3", "green3", "blue", "orange", "red", "magenta3", "darkmagenta"), pch=c(16, 16, 16, 18, 18, 18, 18), cex=c(1, 1, 1, 1, 1, 1, 2) , tree_ID = c("990", "1099", "1377", "970", "8212", "8266", ""), condition=c("C", "C","C", "DE", "DE", "DE", "DE"), stringsAsFactors=F) 


### list of all unique days in 2008 and 2009
ad <- read.table("Data/Unique_days.csv", sep=";")
all.days <- data.frame(days.ch = ad[,], days.nr = as.numeric(strptime(ad[,], "%d.%m.%y")))
head(all.days)
library(stringr)
month.days <- all.days[str_sub(all.days[,1], 1, 2)=="01",]

ad.short <- read.table("Data/Unique_days_short.csv", sep=";")
all.days.short <- data.frame(days.ch = ad.short[,], days.nr = as.numeric(strptime(ad.short[,], "%d.%m.%y")))
head(all.days.short)
library(stringr)
month.days.short <- all.days.short[str_sub(all.days.short[,1], 1, 2)=="01",]

###### files with control genes

AT.id <- read.table("Data/genes_descr_control/best_hit_blast_result_Sl_predicted_exons.txt", sep=",", stringsAsFactors=FALSE)
head(AT.id)

genes.description <- read.table("Data/genes_descr_control/genes_description.xls", sep="\t", header=T)
head(genes.description)

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

### crude analysis - MDS plots

#####################################################################################################

setwd("/home/Shared/data/seq/Shimizu_RNA_seq/")
load("Shimizu_workspace.Rdata")

all(colnames(x)==rownames(new.samps))

library(edgeR)
d <- DGEList(x, group=new.samps$tree_ID)
d <- calcNormFactors(d)

# make sure a gene is expressed (CPM > 1) in more than 2 samples
cps <- cpm(d, normalized.lib.sizes=TRUE)
d <- d[ rowSums(cps>1) > 10, ]


# plots: MDS, 
library(ggplot2)
library(EDASeq)
dir.create("Plots_MDS")


pdf("Plots_MDS/mds_500.pdf",w=10,h=10)
mds <- plotMDS(d, col=new.samps$tree_col, top=500, main="method=logFC, prior.count=2", labels=new.samps$short.name)
mds10 <- plotMDS(d,col=new.samps$tree_col, main="method=logFC, prior.count=10", prior.count=10, top=500, labels=new.samps$short.name)
mdsb <- plotMDS(d, col=new.samps$tree_col, method="bcv", top=500, labels=new.samps$short.name)
dev.off()


pdf("Plots_MDS/mds_500_WP.pdf",w=10,h=10)
mds <- plotMDS(d, col=new.samps$colors.wp, top=500, main="method=logFC, prior.count=2", labels=new.samps$short.name)
mds10 <- plotMDS(d,col=new.samps$colors.wp,main="method=logFC, prior.count=10", prior.count=10, top=500, labels=new.samps$short.name)
mdsb <- plotMDS(d, col=new.samps$colors.wp, method="bcv", top=500, labels=new.samps$short.name)
dev.off()


pdf("Plots_MDS/mds_1000.pdf",w=10,h=10)
mds <- plotMDS(d, col=new.samps$tree_col, top=1000, main="method=logFC, prior.count=2", labels=new.samps$short.name)
mds10 <- plotMDS(d,col=new.samps$tree_col,main="method=logFC, prior.count=10", prior.count=10, top=1000, labels=new.samps$short.name)
mdsb <- plotMDS(d, col=new.samps$tree_col, method="bcv", top=1000, labels=new.samps$short.name)
dev.off()


pdf("Plots_MDS/pca.pdf",w=10,h=10)
plotPCA(d$counts, col = new.samps$tree_col)
dev.off()


pdf("Plots_MDS/rle.pdf",w=15,h=7)
plotRLE(d$counts, col = new.samps$tree_col, outline = FALSE)
plotRLE(cpm(d, normalized.lib.sizes=TRUE), col = new.samps$tree_col, outline = FALSE)
dev.off()


#####################################################################################################

### plots of expression for flowering genes

#####################################################################################################


#AT.genes <- Athaliana.flowering.genes.FT.SVP
AT.genes <- Athaliana.flowering.genes

#elim.samps=c(flowered.samps)
elim.samps=NULL
x <- x[,!names(x) %in% elim.samps]
new.samps <- new.samps[!new.samps$sample_name %in% elim.samps, ]


colnames(x) == rownames(new.samps)


library(edgeR)
d.org <- DGEList(x, group=new.samps$tree_ID)
d.org <- calcNormFactors(d.org)

d.cpm <- cpm(d.org, normalized.lib.sizes=TRUE)
d.cpm.l <- log(d.cpm + min(d.cpm[d.cpm != 0]))


dir.create("Plots_of_flowering_genes/", showWarnings=F, recursive=T)


pdf(paste0("Plots_of_flowering_genes/" , "Flowering_genes_from_S4_table_Nov_fl" ,".pdf"), h=5, w=10)
#info.table <- NULL

for(g in 1:nrow(AT.genes)){
  # g=1
  cat(paste(g, ", "))
  genes <- AT.id[AT.id[,2] == AT.genes[g, 1], 1]
  
  if(length(genes!=0)){
    for(j in 1:length(genes)){
      # j=1
      
      #info.table <- rbind(info.table, data.frame(AT.ID=AT.genes[g, 1], AT.description=AT.genes[g, 2] , ID=genes[j]))
      
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

plot.vars <- c("Water.Potential", "Soil.Moisture")

for(v in plot.vars){
  plot(0, type="n", xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(new.samps[,v])), max(na.omit(new.samps[,v]))), xlab="Time", ylab=v, xaxt = "n")
  for(t in trees.order$legend){
    lines(new.samps$time_nr[new.samps$tree_legend==t], new.samps[new.samps$tree_legend==t, v], col=trees.order$color[trees.order$legend==t], type="b", pch=trees.order$pch[trees.order$legend==t], cex=trees.order$cex[trees.order$legend==t]) 
  }
  axis(side=1, at=month.days[,2], labels=month.days[,1])
  legend("bottomright", legend = trees.order$legend, col=trees.order$color, cex=0.7, text.col=trees.order$color)
  
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


dev.off()


#####################################################################################################

### PLOTS OF SELECTED GENES

#####################################################################################################

# have to load the workspace!
setwd("/home/Shared/data/seq/Shimizu_RNA_seq/")
load("Shimizu_workspace.Rdata")

source("/home/gosia/R/R_Shimizu_RNA_seq/Plot_function.R")



plot.data <- d.cpm

### select genes for plotting
clustering.results <- read.table("Clustering/Kmeans_all/CPM_norm_all_kmeans_table.xls", sep="\t", header=T, stringsAsFactors=F)
tot.nr.cluster <- 32 
nr.cluster <- 24 #12
gene.cluster.ID <- clustering.results[clustering.results[,paste0("clustering", tot.nr.cluster)] == nr.cluster, c("ID","AT_ID")]
gene.cluster.AT_ID <- unique(na.omit(clustering.results[clustering.results[,paste0("clustering", tot.nr.cluster)] == nr.cluster, "AT_ID"]))
MolEcol.clustered.genes <- read.table("Data/genes_descr_control/gene_list_in_all_clusters_DEgenes_Mol_Ecol.csv", sep=";", header=T)
MolEcol.nr.cluster <- 1
gene.set.AT_ID <- unique(MolEcol.clustered.genes[MolEcol.clustered.genes[,2] == MolEcol.nr.cluster,1])
genesAT_ID <- intersect(gene.cluster.AT_ID, gene.set.AT_ID)

genesID <- gene.cluster.ID[gene.cluster.ID$AT_ID %in% genesAT_ID, "ID"]


samps.to.plot <- new.samps$sample_name
out.path <- "Plots_of_selected_genes"
out.name <- "Overlap_with_MolEcol_cluster1"
plot.WP=F
ylab="Gene expression in CPM"



plot.genes(plot.data=plot.data, genesID=genesID, new.samps=new.samps, samps.to.plot=samps.to.plot, genes.full.description=genes.full.description, trees.order=trees.order, month.days=month.days, out.path=out.path, out.name=out.name, plot.WP=plot.WP, ylab=ylab)



# plot negative control genes


setwd("/home/Shared/data/seq/Shimizu_RNA_seq/")
load("Shimizu_workspace.Rdata")

source("/home/gosia/R/R_Shimizu_RNA_seq/Plot_function.R")

all(colnames(x)==rownames(new.samps))

library(edgeR)
d <- DGEList(x, group=new.samps$tree_ID)
d <- calcNormFactors(d)

# make sure a gene is expressed (CPM > 1) in more than 2 samples
cps <- cpm(d, normalized.lib.sizes=TRUE)



plot.data <- cps
genesID <- neg.contrl.ID

samps.to.plot <- new.samps$sample_name
out.path <- "Plots_RUN_RUVSeq"
out.name <- "Negative_control_genes_Housekeeping"
plot.WP=F
ylab="Gene expression in CPM"



plot.genes(plot.data=plot.data, genesID=genesID, new.samps=new.samps, samps.to.plot=samps.to.plot, genes.full.description=genes.full.description, trees.order=trees.order, month.days=month.days, out.path=out.path, out.name=out.name, plot.WP=plot.WP, ylab=ylab)




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
x.orig <- x.orig[,new.samps.orig$sample_name]
x <- x.orig[,!colnames(x.orig) %in% elim.samps]
new.samps <- new.samps.orig[!rownames(new.samps.orig) %in% elim.samps, ]

all(colnames(x)==rownames(new.samps))

library(limma)
library(edgeR)

d.org <- DGEList(x, group=new.samps$tree_ID)

# normalization between samples, default = "TMM"
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


write.table(data.frame(contigs=rownames(d.org$counts) , d.org$counts), "Clustering/Data_original.xls", quote=FALSE, sep="\t", row.names=FALSE)


write.table(data.frame(contigs=rownames(d.cpm.org) , d.cpm.org), "Clustering/Data_normalized_cpm.xls", quote=FALSE, sep="\t", row.names=FALSE)



############################################
# Kmeans run
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



source("/home/gosia/R/R_Shimizu_RNA_seq/Kmeans_clustering.R")
# source("/home/gosia/R/R_Shimizu_RNA_seq/Kmeans_clustering_tmp.R")

pam.x <- count.data$d.cpm.org
#pam.x <- pam.x[1:200,]
dim(pam.x)



# all
kmeans.clustering(pam.x, new.samps, genes.full.description, colors.wp, norm.method="norm", out.path="Clustering/Kmeans_all2/", out.name="CPM_norm_all", prior.nr.cl=10:50, mc.cores=15, ylim=c(-2, 5), clustering.samps="all", colors=c("yellow","blue"))



# drought trees (no flowering)
drought.samp <- new.samps[new.samps$tree_ID %in% c("970", "8266"), "sample_name"]
# no flowering sample E7_8266_20090416
drought.samp <- drought.samp[!drought.samp %in% c("E7_8266_20090416")]

kmeans.clustering(pam.x[,drought.samp], new.samps[drought.samp,], genes.full.description, colors.wp, norm.method="norm", out.path="Clustering/Kmeans_drought_trees/", out.name="CPM_norm_drought", prior.nr.cl=5:40, mc.cores=15, ylim=c(-2, 5), clustering.samps="all", colors=c("yellow","blue"))




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


go.name <- "Kmeans_all"
whole.clustering <- read.table("Clustering/Kmeans_all/CPM_norm_all_kmeans_table.xls", sep="\t", header=T)
nr.clusters <- 32


# go.name <- "Kmeans_drought_trees"
# whole.clustering <- read.table("Clustering/Kmeans_drought_trees/CPM_norm_drought_kmeans_table.xls", sep="\t", header=T)
# nr.clusters <- 22



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
  
  # cluster <- 1
  sel.cluster.genes <- whole.clustering$AT_ID[whole.clustering[, paste0("clustering", nr.clusters)]==cluster]
  
  gene.vector=as.numeric(assayed.genes %in% sel.cluster.genes)
  names(gene.vector) <- assayed.genes
  names(assayed.genes) <- gene.vector
  
  
  for(go in c("BP","MF","CC")){
    # go = "BP"
    
    sampleGOdata <- new("topGOdata", description = "Simple session", ontology = go, allGenes = gene.vector, geneSel = fun.gene.sel , nodeSize = 10, annot = annFUN.org, mapping = "org.At.tair.db")
    # Fisher's exact test which is based on gene counts, and a Kolmogorov-Smirnov like test which computes enrichment based on gene scores
    # For the method classic each GO category is tested independently
    resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
    resultFisher.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher")
    
    
    pValues.Fisher <- score(resultFisher)
    #topNodes.Fisher <- sum(pValues.Fisher < 0.05)
    topNodes.Fisher <- length(pValues.Fisher)
    pValues.Fisher.elim <- score(resultFisher.elim)
    #topNodes.Fisher.elim <- sum(pValues.Fisher.elim < 0.0001)
    topNodes.Fisher.elim <- length(pValues.Fisher.elim)
    
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
    
   
      allRes.Fisher <- GenTable(sampleGOdata, classicFisher = resultFisher, elimFisher = resultFisher.elim, orderBy = "classicFisher", ranksOf = "elimFisher", topNodes = topNodes.Fisher)      
      allRes.Fisher$GO <- go
      allRes.Fisher$cluster <- cluster      
      allRes.Fisher.merged <- rbind(allRes.Fisher.merged, allRes.Fisher)
      

      allRes.Fisher.elim <- GenTable(sampleGOdata, classicFisher = resultFisher, elimFisher = resultFisher.elim, orderBy = "elimFisher", ranksOf = "classicFisher", topNodes = topNodes.Fisher.elim)      
      allRes.Fisher.elim$GO <- go
      allRes.Fisher.elim$cluster <- cluster   
      allRes.Fisher.elim.merged <- rbind(allRes.Fisher.elim.merged, allRes.Fisher.elim)

    
    
    #     pdf(paste("GO/GO_", go, ".pdf", sep=""), width = 10, height = 10)
    #     showSigOfNodes(sampleGOdata, score(resultKS), firstSigNodes = 5, useInfo = 'all')
    #     dev.off()
    
    
  }
  
  
}



allRes.Fisher.merged$classicFisher.adj <- p.adjust(p=allRes.Fisher.merged$classicFisher, method = "BH")


write.table(allRes.Fisher.merged, paste("GO/",go.name,"_GO_Fisher" , ".xls", sep=""), sep="\t", row.names=F)

#try(write.table(allRes.Fisher.merged[allRes.Fisher.merged$classicFisher.adj <= 0.05, ], paste("GO/",go.name,"_GO_Fisher_FDR0.05" , ".xls", sep=""), sep="\t", row.names=F))

allRes.Fisher.elim.merged$elimFisher.adj <- p.adjust(p=allRes.Fisher.elim.merged$elimFisher, method = "BH")

write.table(allRes.Fisher.elim.merged, paste("GO/",go.name,"_GO_Fisher_elim" , ".xls", sep=""), sep="\t", row.names=F)

#try(write.table(allRes.Fisher.elim.merged[allRes.Fisher.elim.merged$elimFisher.adj <= 0.05, ], paste("GO/",go.name,"_GO_Fisher_elim_FDR0.05" , ".xls", sep=""), sep="\t", row.names=F))


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


go.name <- "Kmeans_all"
whole.clustering.org <- read.table("Clustering/Kmeans_all/CPM_norm_all_kmeans_table.xls", sep="\t", header=T)
nr.clusters <- 32


# go.name <- "Kmeans_drought_trees"
# whole.clustering.org <- read.table("Clustering/Kmeans_drought_trees/CPM_norm_drought_kmeans_table.xls", sep="\t", header=T)
# nr.clusters <- 22



whole.clustering <- whole.clustering.org[!is.na(whole.clustering.org$AT_ID), c("ID", "AT_ID", paste0("clustering", nr.clusters))]
colnames(whole.clustering) <- c("ID", "AT_ID","clustering")


gene.set <- Athaliana.flowering.genes[,1]
name.out <- "Athaliana_flowering_genes"


fisher.pvalues <- matrix(0, 1 ,nr.clusters)
colnames(fisher.pvalues) <- paste0( "CL", 1:nr.clusters)


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
  fisher.pvalues[1, cluster] <- fisher.out$p.value
  
  if(any(fisher.x < 10)){    
    fisher.table[cluster, "pvalues"] <- NA
    fisher.pvalues[1, cluster] <- NA
  }

  
  fisher.table[cluster, "X1"] <- fisher.x[1,1]
  fisher.table[cluster, "X2"] <- fisher.x[1,2]
  fisher.table[cluster, "X3"] <- fisher.x[2,1]
  fisher.table[cluster, "X4"] <- fisher.x[2,2]
  fisher.table[cluster,  "cluster.size"] <- fisher.x[1,1] + fisher.x[2,1]
  
  
}

fisher.table <- data.frame(cluster=1:nr.clusters, fisher.table)

write.table(fisher.table, paste0(dir.out,"/", name.out, "_", go.name, "_Fisher_pvalues.xls"), sep="\t", row.names=F)


fisher.pvalues.org <- fisher.pvalues

fisher.pvalues <- fisher.pvalues.org


fisher.pvalues.adj <- p.adjust(p=fisher.pvalues, method = "BH")
fisher.pvalues.adj <- matrix(fisher.pvalues.adj, dim(fisher.pvalues)[1], dim(fisher.pvalues)[2])
colnames(fisher.pvalues.adj) <- paste0( "CL", 1:nr.clusters)

fisher.pvalues.adj <- data.frame(gene.set="Athaliana.flowering.genes", fisher.pvalues.adj)

fisher.pvalues <- data.frame(gene.set="Athaliana.flowering.genes", fisher.pvalues)

write.table(fisher.pvalues, paste0(dir.out,"/_ALL_PVALUES_", go.name, "_Fisher_pvalues.xls"), sep="\t", row.names=F)

write.table(fisher.pvalues.adj, paste0(dir.out,"/_ALL_PVALUES_", go.name, "_Fisher_pvalues_ADJ.xls"), sep="\t", row.names=F)



##############################
### gene sets
##############################

setwd("/home/Shared/data/seq/Shimizu_RNA_seq/")
load("Shimizu_workspace.Rdata")

dir.out <- "Fisher_test/Gene_sets"

dir.create(dir.out, showWarnings=F, recursive=T)


# go.name <- "Kmeans_all"
# whole.clustering <- read.table("Clustering/Kmeans_all/CPM_norm_all_kmeans_table.xls", sep="\t", header=T)
# nr.clusters <- 32


go.name <- "Kmeans_drought_trees"
whole.clustering <- read.table("Clustering/Kmeans_drought_trees/CPM_norm_drought_kmeans_table.xls", sep="\t", header=T)
nr.clusters <- 22


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


# go.name <- "Kmeans_all"
# whole.clustering <- read.table("Clustering/Kmeans_all/CPM_norm_all_kmeans_table.xls", sep="\t", header=T)
# nr.clusters <- 32


go.name <- "Kmeans_drought_trees"
whole.clustering <- read.table("Clustering/Kmeans_drought_trees/CPM_norm_drought_kmeans_table.xls", sep="\t", header=T)
nr.clusters <- 22


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



##############################
### MolEcol clusters: Sh.l <-> Sh.b
##############################


setwd("/home/Shared/data/seq/Shimizu_RNA_seq/")
load("Shimizu_workspace.Rdata")

dir.out <- "Fisher_test/MolEcol_clusters_ShL_ShB"

dir.create(dir.out, showWarnings=F, recursive=T)


# go.name <- "Kmeans_all_unigenesONLY"
# whole.clustering.org <- read.table("Clustering/Kmeans_all/CPM_norm_all_kmeans_table.xls", sep="\t", header=T)
# nr.clusters <- 32


go.name <- "Kmeans_drought_trees_unigenesONLY"
whole.clustering.org <- read.table("Clustering/Kmeans_drought_trees/CPM_norm_drought_kmeans_table.xls", sep="\t", header=T)
nr.clusters <- 22


whole.clustering <- whole.clustering.org[!is.na(whole.clustering.org$ID), c("ID", "AT_ID", paste0("clustering", nr.clusters))]
colnames(whole.clustering) <- c("ID", "AT_ID","clustering")
head(whole.clustering)


MolEcol.clustered.genes <- read.table("Data/genes_descr_control/gene_list_in_all_clusters_unigenes_DEgenes_Mol_Ecol.csv", header=T, sep=",")
head(MolEcol.clustered.genes)


unigene.match <- read.table("Data/genes_descr_control/reciprocal_best_hit_Sb_vs_Sl.txt", sep="\t")
head(unigene.match)


sum(whole.clustering[, 1] %in% unigene.match[,2])


whole.clustering2 <- merge(whole.clustering, unigene.match, by.x=1, by.y=2, all.x=TRUE)
names(whole.clustering2)[4] <- "unigene_ID"

sum(!is.na(whole.clustering2[,4]))


# by AT_ID
sum(whole.clustering2[, 2] %in% MolEcol.clustered.genes[, 2])

length(unique(whole.clustering2[whole.clustering2[, 2] %in% MolEcol.clustered.genes[, 2], 2]))

# by unigene
sum(whole.clustering2[, 4] %in% MolEcol.clustered.genes[, 1])

length(unique(whole.clustering2[whole.clustering2[, 4] %in% MolEcol.clustered.genes[, 1], 4]))

###


whole.clustering2 <- whole.clustering2[!is.na(whole.clustering2$unigene_ID),]


fisher.pvalues <- matrix(0, max(MolEcol.clustered.genes[,3]), nr.clusters)
colnames(fisher.pvalues) <- paste0( "CL", 1:nr.clusters)


for(i in 1:max(MolEcol.clustered.genes[,3])){
  # i=1
  
  gene.set <- MolEcol.clustered.genes[MolEcol.clustered.genes[,3]==i,1]
  name.out <- paste0("MolEcol_cluster", i)
  
  
  fisher.table <- matrix(0,nr.clusters,6)
  colnames(fisher.table) <- c("X1", "X2", "X3", "X4", "pvalues", "cluster.size")
  
  for(cluster in 1:nr.clusters){
    
    # cluster <- 1
    
    clustering.tmp <- whole.clustering2
    
    clustering.tmp$gene.set <- ifelse(whole.clustering2$unigene_ID %in%  gene.set, 1, 2)
    clustering.tmp$clustering <- ifelse(whole.clustering2$clustering==cluster, 1, 2)
    
    
    fisher.x <- as.matrix(with(clustering.tmp, table(gene.set, clustering)))
    
    ### !!! very inportant to set up alternative hipothesis
    fisher.out <- fisher.test(fisher.x, alternative="greater")
    
    fisher.table[cluster, "pvalues"] <- fisher.out$p.value
    fisher.pvalues[i, cluster] <- fisher.out$p.value
    
    if(any(fisher.x <= 10)){
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

fisher.pvalues.adj <- data.frame(gene.set=paste0("MolEcol_cluster", 1:max(MolEcol.clustered.genes[,3])), fisher.pvalues.adj)

fisher.pvalues <- data.frame(gene.set=paste0("MolEcol_cluster", 1:max(MolEcol.clustered.genes[,3])), fisher.pvalues)

write.table(fisher.pvalues, paste0(dir.out,"/_ALL_PVALUES_", go.name, "_Fisher_pvalues.xls"), sep="\t", row.names=F)

write.table(fisher.pvalues.adj, paste0(dir.out,"/_ALL_PVALUES_", go.name, "_Fisher_pvalues_ADJ.xls"), sep="\t", row.names=F)



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





































