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



