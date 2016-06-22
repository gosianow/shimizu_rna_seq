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
# BioC 3.0
### preparation of variables describing samples 

#####################################################################################################

RPath <- "/home/gosia/R/R_Shimizu_RNA_seq/Shimizu_analysis_2014-12-03/"
dataPath <- "/home/Shared/data/seq/Shimizu_RNA_seq/Data/"

analysisPath <- "Analysis_2014-12-03"
analysisPath <- paste0("/home/Shared/data/seq/Shimizu_RNA_seq/", analysisPath)
dir.create(analysisPath, showWarnings = FALSE)
setwd(analysisPath)



####################################################
### samps
####################################################


dir.create("Samples_out", recursive = TRUE)

### parse sample information
samps <- read.table(paste0(dataPath, "Samples/sample_list.csv"), sep=";", stringsAsFactors=FALSE, header=TRUE)
head(samps)
names(samps)

### CORRECT THE DATE!!! (old dates are wrong)
samps$year.month.day.old <- samps$year.month.day
samps$year.month.day <- paste(substr(samps$sample_name, nchar(samps$sample_name)-1, nchar(samps$sample_name)), substr(samps$sample_name, nchar(samps$sample_name)-3, nchar(samps$sample_name)-2), substr(samps$sample_name, nchar(samps$sample_name)-5, nchar(samps$sample_name)-4), sep=".")

head(samps[,c("sample_name","year.month.day", "year.month.day.old")])


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
ad <- read.table(paste0(dataPath, "Unique_days.csv"), sep=";")
all.days <- data.frame(days.ch = ad[,], days.nr = as.numeric(strptime(ad[,], "%d.%m.%y")))
head(all.days)
library(stringr)
month.days <- all.days[str_sub(all.days[,1], 1, 2)=="01",]

ad.short <- read.table(paste0(dataPath, "Unique_days_short.csv"), sep=";")
all.days.short <- data.frame(days.ch = ad.short[,], days.nr = as.numeric(strptime(ad.short[,], "%d.%m.%y")))
head(all.days.short)
library(stringr)
month.days.short <- all.days.short[str_sub(all.days.short[,1], 1, 2)=="01",]


####################################################
### Rainfall & water deficit
####################################################


wd <- read.table(paste0(dataPath, "MeteorologicalData/Lambir_meteorological_data_200710-201001mod.txt"), stringsAsFactors=FALSE, header=TRUE, sep = "\t")
head(wd)


wd$time_ch <- strptime(wd$Date, "%d/%m/%y")
wd$time_nr <- as.numeric(wd$time_ch)


pdf("Plots_Samples/Water_deficit_raw.pdf", width = 10, height = 5)
plot(wd$time_nr, wd$evaporation...rain..mm.30.days. , main="Water Deficit", xlab="Time", ylab="30-day moving total of water deficit (mm)", type="l", col="brown", lwd = 4)

plot(wd$time_nr, wd$evaporation...rain..mm.30.days. , main="Water Deficit", xlab="Time", ylab="30-day moving total of water deficit (mm)", type="l", col="brown", lwd = 4,  xaxt = "n" )
axis(side=1, at=month.days[,2], labels=month.days[,1])
abline(v=new.samps$time_nr, col="grey")

plot(wd$time_nr, wd$evaporation...rain..mm.30.days. , main="Water Deficit", xlab="Time", ylab="30-day moving total of water deficit (mm)", type="l", col="brown", lwd = 4 ,xlim=c(min(all.days.short[,2]), max(all.days.short[,2])),  xaxt = "n")
axis(side=1, at=month.days[,2], labels=month.days[,1])
abline(v=new.samps$time_nr, col="grey")

dev.off()



pdf("Plots_Samples/Rain_raw.pdf", width = 10, height = 5)

plot(wd$time_nr, wd$Rain.mm.30.days. , main="Rainfall", xlab="Time", ylab="30-day moving total of rainfall (mm)", type="l", col="darkblue", lwd = 4 )

plot(wd$time_nr, wd$Rain.mm.30.days. , main="Rainfall", xlab="Time", ylab="30-day moving total of rainfall (mm)", type="l", col="darkblue", lwd = 4,  xaxt = "n" )
axis(side=1, at=month.days[,2], labels=month.days[,1])
abline(v=new.samps$time_nr, col="grey")


plot(wd$time_nr, wd$Rain.mm.30.days. , main="Rainfall", xlab="Time", ylab="30-day moving total of rainfall (mm)", type="l", col="darkblue", lwd = 4, xlim=c(min(all.days.short[,2]), max(all.days.short[,2])), xaxt = "n" )
axis(side=1, at=month.days[,2], labels=month.days[,1])
abline(v=new.samps$time_nr, col="grey")

dev.off()


####################################################
### Soild Moisture
####################################################

### get soil moisture data
sm <- read.table(paste0(dataPath, "MeteorologicalData/Soil_Moisture.csv"), sep=";", stringsAsFactors=FALSE, header=TRUE)
head(sm)
sm <- sm[1:271,]

### time format change
sm$Year <- ifelse(sm$DayUniq <= 366, 2008, 2009)
sm$year.month.day <- paste(sm$Day,sm$Month, sm$Year, sep=".")
sm$time_ch <- strptime(sm$year.month.day, "%d.%m.%Y")
sm$time_nr <- as.numeric(sm$time_ch)


### plots of raw data
dir.create(path="Plots_Samples", showWarnings=FALSE, recursive=TRUE)

pdf("Plots_Samples/Soil_moisture.pdf", width = 10, height = 5)
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


pdf("Plots_Samples/Soil_moisture_DE970.pdf", width = 10, height = 5)
plot(sm$time_nr, sm$DE970, type="l", col=trees.order$color[trees.order$tree_ID=="970"], ylim=c(0,1.1), xlim=c(min(all.days.short[,2]), max(all.days.short[,2])), main="Tree 970", ylab="Soil moisture", xlab="Time",  xaxt = "n", lwd = 4, las=1, cex.lab=1.5)
lines(sm$time_nr, sm$C970, type="l", col=trees.order$color[trees.order$tree_ID=="970"], lty=3, lwd = 4)
axis(side=1, at=month.days[,2], labels=month.days[,1])
abline(v=new.samps$time_nr[new.samps$tree_ID==970], col="grey")
dev.off()


pdf("Plots_Samples/Soil_moisture_DE8266.pdf", width = 10, height = 5)
plot(sm$time_nr, sm$DE8266, type="l", col=trees.order$color[trees.order$tree_ID=="8266"], ylim=c(0,1.1), xlim=c(min(all.days.short[,2]), max(all.days.short[,2])), main="Tree 8266", ylab="Soil moisture", xlab="Time",  xaxt = "n", lwd=4, las=1, cex.lab=1.5)
lines(sm$time_nr, sm$C8266, col=trees.order$color[trees.order$tree_ID=="8266"], lty=3, lwd=4)
axis(side=1, at=month.days[,2], labels=month.days[,1])
abline(v=new.samps$time_nr[new.samps$tree_ID==8266], col="grey")
dev.off()



pdf("Plots_Samples/Soil_moisture_DE8212.pdf", width = 10, height = 5)
plot(sm$time_nr, sm$DE8212, type="l", col=trees.order$color[trees.order$tree_ID=="8212"], ylim=c(0,1.1), xlim=c(min(all.days.short[,2]), max(all.days.short[,2])), main="Tree 8212", ylab="Soil moisture", xlab="Time",  xaxt = "n", lwd=4, las=1, cex.lab=1.5)
lines(sm$time_nr, sm$C8212a, col=trees.order$color[trees.order$tree_ID=="8212"], lty=3, lwd=4)
lines(sm$time_nr, sm$C8212b, col=trees.order$color[trees.order$tree_ID=="8212"], lty=3, lwd=4)
axis(side=1, at=month.days[,2], labels=month.days[,1])
abline(v=new.samps$time_nr[new.samps$tree_ID==8212], col="grey")
dev.off()



# Interpolate and roll 

source(paste0(RPath, "Interpolate_and_roll.R"))

roll.days=c(14, 21, 28)

new.samps.DE970 <- interpolate.and.roll(intrp.x=sm$time_nr, intrp.y=sm$DE970, intrp.points=all.days[,2], table1=new.samps[new.samps$tree_ID == 970,], roll.days=roll.days, col.values="Soil.Moisture", plot.name="Soil_MoistureDE970", color=trees.order$color[trees.order$tree_ID=="970"], month.days=month.days, all.days.short=all.days.short)

new.samps.DE8212 <- interpolate.and.roll(intrp.x=sm$time_nr, intrp.y=sm$DE8212, intrp.points=all.days[,2], table1=new.samps[new.samps$tree_ID == 8212,], roll.days=roll.days, col.values="Soil.Moisture", plot.name="Soil_MoistureDE8212", color=trees.order$color[trees.order$tree_ID=="8212"], month.days=month.days, all.days.short=all.days.short)

new.samps.DE8266 <- interpolate.and.roll(intrp.x=sm$time_nr, intrp.y=sm$DE8266, intrp.points=all.days[,2], table1=new.samps[new.samps$tree_ID == 8266,], roll.days=roll.days, col.values="Soil.Moisture", plot.name="Soil_MoistureDE8266",  color=trees.order$color[trees.order$tree_ID=="8266"], month.days=month.days, all.days.short=all.days.short)


new.samps.control <- new.samps[new.samps$drough.control == "control",]
Soil.moisture.control <- data.frame(time_nr = new.samps.control$time_nr, matrix(NA, nrow = nrow(new.samps.control), ncol = (length(roll.days) + 1) ))
colnames(Soil.moisture.control) <- c("time_nr", paste0("Soil.Moisture", c("", roll.days)))

new.samps.control <- merge(new.samps.control, Soil.moisture.control, by="time_nr")

new.samps.sm <- rbind(new.samps.DE970, new.samps.DE8212, new.samps.DE8266, new.samps.control)

new.samps <- new.samps.sm

names(new.samps)


####################################################
### Water Potential
####################################################

### prepare Mean.Water.Potential var 

wp <- read.table(paste0(dataPath, "MeteorologicalData/Water_potential_from_Inoue.csv"), sep=";", header=T)
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
  
  pdf(paste("Plots_Samples/Water_potential_", exp.tree.ids[i], ".pdf", sep=""), width = 10, height = 5)
  plot(wp$time_nr[!is.na(wp[,tree.sampl])], wp[!is.na(wp[,tree.sampl]), tree.sampl], type="l", col=trees.order$color[trees.order$tree_ID==tree.ids[i]], ylim=c(0.2,0.8), xlim=c(min(all.days.short[,2]), max(all.days.short[,2])), main=paste0("Tree ",exp.tree.ids[i]), ylab="Mean Water Potential", xlab="Time", xaxt = "n", las=1, cex.lab=1.5, lwd=4)
  axis(side=1, at=month.days[,2], labels=month.days[,1])
  abline(v=new.samps$time_nr[new.samps$tree_ID==tree.ids[i]], col="grey")
  dev.off()
  
  
}


pdf(paste("Plots_Samples/Water_potential.pdf", sep=""), width = 10, height = 5)
i=1
tree.sampl <- paste(exp.tree.ids[i], "sMean", sep="")
plot(wp$time_nr[!is.na(wp[,tree.sampl])], wp[!is.na(wp[,tree.sampl]), tree.sampl], type="l", col=trees.order$color[trees.order$tree_ID==tree.ids[i]], ylim=c(0.2,0.8), xlim=c(min(wp$time_nr), max(wp$time_nr)), main="Water Potential", ylab="Mean Water Potential", xlab="Time", xaxt = "n", las=1, cex.lab=1.5, lwd=3)
axis(side=1, at=month.days[,2], labels=month.days[,1])
for(i in 2:length(tree.ids)){
  # i=5
  tree.sampl <- paste(exp.tree.ids[i], "sMean", sep="")
  lines(wp$time_nr[!is.na(wp[,tree.sampl])], wp[!is.na(wp[,tree.sampl]), tree.sampl], col=trees.order$color[trees.order$tree_ID==tree.ids[i]], lwd=3)
}
legend("bottomleft", trees.order$legend[1:6], col=trees.order$color[1:6], lty=rep(1, 6) , cex=0.8, lwd=4)
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
  new.samps.wp <- rbind(new.samps.wp, new.samps.wp.tmp <- interpolate.and.roll(intrp.x=wp$time_nr[!is.na(wp[,tree.sampl])], intrp.y=wp[!is.na(wp[,tree.sampl]),tree.sampl], intrp.points=all.days[,2], table1=new.samps[new.samps$tree_ID == tree.ids[i],] , roll.days=c(14, 21, 28), col.values="Water.Potential", plot.name=paste0("Water_Potential",tree.sampl), month.days=month.days, all.days.short=all.days.short, color=trees.order$color[trees.order$tree_ID==tree.ids[i]]))
  
}

names(new.samps.wp)
head(new.samps.wp)


new.samps <- new.samps.wp


####################################################
### Temperature
####################################################

tempr <- read.table(paste0(dataPath, "MeteorologicalData/Lambir_temperature_data.csv"), head=T, sep=";")
head(tempr)

library(stringr)

tempr$time_ch <- strptime(substring(tempr$Time, 1, 8), "%d.%m.%y")
tempr$time_nr <- as.numeric(tempr$time_ch)

tempr$time_hour <- substring(tempr$Time, 10, 11)
tempr$time_hour_nr <- as.numeric(tempr$time_hour)

pdf("Plots_Samples/Temperature_raw.pdf", width = 10, height = 5)
plot(tempr$time_nr, tempr$Temp, pch=".", main="Temperature", xlab="Time", ylab="Temperature")
#abline(v=new.samps$time_nr, col="grey")
dev.off()

pdf("Plots_Samples/Temperature_over_day.pdf", width = 10, height = 5)
boxplot(Temp ~ time_hour, data=tempr, xlab="Hour", ylab="Temperature")
dev.off()


### calculate average temp per day 
days <- unique(tempr$time_ch)
tempr.day <- data.frame(time_ch=days, time_nr=as.numeric(days) , TempAvg=0, TempMax=0, TempMin=0)

for(i in 1:nrow(tempr.day)){
  # i=1
  tempr.day[i , "TempAvg"] <- mean(tempr$Temp[tempr$time_ch==days[i] & tempr$time_hour_nr %in% 10:18])
  tempr.day[i , "TempMin"] <- min(tempr$Temp[tempr$time_ch==days[i] & tempr$time_hour_nr %in% 10:18])
}


### lowess 
tempr.l <- list()

tempr.l[["TempAvg"]] <- lowess(tempr.day[,"time_nr"], tempr.day[,"TempAvg"], f=0.1)

pdf("Plots_Samples/Temperature_Avg_lowess.pdf", width = 10, height = 5)
plot(tempr.day[,"time_nr"], tempr.day[,"TempAvg"], col="grey", pch=20, main="Control trees", ylab="Soil moisture", xlab="Time",  xaxt = "n")
axis(side=1, at=month.days[,2], labels=month.days[,1])
lines(tempr.l[["TempAvg"]], lwd=2, col="darkred")
dev.off()


### some plots of avg, max, min temperature
pdf("Plots_Samples/Temperature_Avg.pdf", width = 10, height = 5)
plot(tempr.day$time_nr, tempr.day$TempAvg, pch=20, main="Average Temperature", xlab="Time", ylab="Temperature",xlim=c(min(all.days.short[,2]), max(all.days.short[,2])),  xaxt = "n")
axis(side=1, at=month.days[,2], labels=month.days[,1])
lines(tempr.l[["TempAvg"]], lwd=2, col="darkred")
abline(v=new.samps$time_nr, col="grey")
dev.off()


tempr.l[["TempMin"]] <- lowess(tempr.day[,"time_nr"], tempr.day[,"TempMin"], f=0.1)

pdf("Plots_Samples/Temperature_Min_lowess.pdf", width = 10, height = 5)
plot(tempr.day[,"time_nr"], tempr.day[,"TempMin"], col="grey", pch=20, main="Control trees", ylab="Soil moisture", xlab="Time",  xaxt = "n")
axis(side=1, at=month.days[,2], labels=month.days[,1])
lines(tempr.l[["TempMin"]], lwd=2, col="darkred")
dev.off()


### some plots of avg, max, min temperature
pdf("Plots_Samples/Temperature_Min.pdf", width = 10, height = 5)
plot(tempr.day$time_nr, tempr.day$TempMin, pch=20, main="Min Temperature", xlab="Time", ylab="Temperature",xlim=c(min(all.days.short[,2]), max(all.days.short[,2])),  xaxt = "n")
axis(side=1, at=month.days[,2], labels=month.days[,1])
lines(tempr.l[["TempMin"]], lwd=2, col="darkred")
abline(v=new.samps$time_nr, col="grey")
dev.off()



tempr.ll <- data.frame(time_nr = tempr.l$TempAvg$x, TempAvg = tempr.l$TempAvg$y, TempMin = tempr.l$TempMin$y)
head(tempr.ll)

### rolled means of temperature

library(zoo)

roll.days = c(14, 21, 28)
variables = c("TempAvg", "TempMin")
plots.path = "Plots_Samples_Interpolate"

new.samps <- merge(new.samps, tempr.ll, by="time_nr", all.x = TRUE)

for(v in variables){
  # v="TempAvg"

  pdf( paste0(plots.path, "/InterpRoll_", v,".pdf", sep="" ) , width = 10, height = 5)
  
  for(r in roll.days){
    # r=14
    name <- paste(v, r, sep="")
    tempr.day[, name] <- NA
    tempr.day[r:length(tempr.day[, name]), name] <- rollmean(tempr.day[,v], r)

    plot(tempr.day$time_nr, tempr.day[,v], col=1, pch=20, ylab=v, xlab="Time", xlim=c(min(all.days.short[,2]), max(all.days.short[,2])), xaxt = "n", cex.lab=1.5, las=1, main=paste0("Roll over ", r, " days"))
    axis(side=1, at=month.days[,2], labels=month.days[,1])
    abline(v=new.samps$time_nr, col="grey")

    lines(tempr.l[[v]], lwd=1, col="darkred", lty=3)

      lines(tempr.day$time_nr, tempr.day[,name], col="darkred", lty=1, lwd=4)

  }
  
  
  dev.off()
  
  new.samps <- merge(new.samps, tempr.day[,c("time_nr", paste(v, roll.days, sep=""))], by="time_nr")
}

head(new.samps)


new.samps <- unique(new.samps)

dir.create(path="Samples_out", showWarnings=F, recursive=T)

# save(new.samps, file="Samples_out/new_samps_interpolation.RData")

write.table(new.samps, "Samples_out/new_samps_interpolation.xls", sep="\t",  row.names = F)





#####################################################################################################

### Load data 

#####################################################################################################


RPath <- "/home/gosia/R/R_Shimizu_RNA_seq/Shimizu_analysis_2014-12-03/"
dataPath <- "/home/Shared/data/seq/Shimizu_RNA_seq/Data/"

analysisPath <- "Analysis_2014-12-03"
analysisPath <- paste0("/home/Shared/data/seq/Shimizu_RNA_seq/", analysisPath)
dir.create(analysisPath)
setwd(analysisPath)


### samples & factors
new.samps <- read.table("Samples_out/new_samps_interpolation.xls", header=TRUE, sep="\t", stringsAsFactors=FALSE)

rownames(new.samps) <- new.samps$sample_name
head(new.samps)


### list of all unique days in 2008 and 2009
ad <- read.table(paste0(dataPath, "/Unique_days.csv"), sep=";")
all.days <- data.frame(days.ch = ad[,], days.nr = as.numeric(strptime(ad[,], "%d.%m.%y")))
head(all.days)
library(stringr)
month.days <- all.days[str_sub(all.days[,1], 1, 2)=="01",]

ad.short <- read.table(paste0(dataPath, "Unique_days_short.csv"), sep=";")
all.days.short <- data.frame(days.ch = ad.short[,], days.nr = as.numeric(strptime(ad.short[,], "%d.%m.%y")))
head(all.days.short)
library(stringr)
month.days.short <- all.days.short[str_sub(all.days.short[,1], 1, 2)=="01",]


# colors reprezening WP level 

colors.wp <- new.samps[, c("sample_name", "Water.Potential")]
colors.wp$Water.Potential <- 100*round(colors.wp$Water.Potential, digits=2)
colors.palette <- as.character( colorRampPalette(c("red","blue"))(diff(range(na.omit(colors.wp$Water.Potential)))+1))
colors.wp$colors.wp <- colors.palette[colors.wp$Water.Potential - min(na.omit(colors.wp$Water.Potential))+1]
colors.wp$colors.wp[is.na(colors.wp$colors.wp)] <- "#838B8B" 
rownames(colors.wp) <- colors.wp$sample_name
colors.wp <- colors.wp[,"colors.wp", drop=F]

new.samps <- merge(new.samps, colors.wp, by=0, all.x=T)[,-1]
rownames(new.samps) <- new.samps$sample_name



# colors reprezening WP level 
colors.palette <- as.character( colorRampPalette(c("orange","darkgreen"))(nrow(all.days.short)))
names(colors.palette) <- all.days.short$days.nr

colors.time <- new.samps[, c("sample_name", "time_nr")]

colors.time$colors.time <- colors.palette[as.character(colors.time$time_nr)]

colors.time <- colors.time[,"colors.time", drop=F]

new.samps <- merge(new.samps, colors.time, by=0, all.x=T)[,-1]
rownames(new.samps) <- new.samps$sample_name




### counts
x <- read.table(paste0(dataPath, "Data_2013-07-01/raw_data-clean.csv"), sep=",", header=T, row.names=1)
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







###### AT gene matches 

AT.id <- read.table(paste0(dataPath, "Data_2013-07-01/best_hit_blast_result_Sl_predicted_exons.txt"), sep=",", stringsAsFactors=FALSE)
head(AT.id)



###### gene descriptions

library(org.At.tair.db)

syms <- mget(AT.id[,2], org.At.tairSYMBOL)
syms <- sapply(syms, function(u) if(is.na(u[1])) NA else paste(u,collapse=","))
AT.id$at_symbol <- syms

tfd <- read.table(paste0(dataPath, "Data_2013-07-01/TAIR10_functional_descriptions"),sep="\t", header=TRUE, quote="", comment.char="", stringsAsFactors=FALSE)

t10id <- gsub("\\.[1-9]","",tfd$Model_name)
m <- match(AT.id[,2],t10id)

AT.id$description <- as.character(tfd$Short_description[m])


genes.description <- merge(data.frame(rownames(x)), AT.id, by=1, all.x=T)

genes.description <- genes.description[,c(1,2,6,7)]

colnames(genes.description) <- c( "ID", "AT_ID", "AT_symbol", "Description")

write.table(genes.description, paste0(dataPath, "Data_2013-07-01/genes_description.xls"), sep="\t", quote=F, row.names=F)


genes.description <- read.table(paste0(dataPath, "Data_2013-07-01/genes_description.xls"), sep="\t", header=T)
head(genes.description)



### Flowering genes - Table S4 from MolEcol paper
Athaliana.flowering.genes <- read.table(paste0(dataPath,"GeneControlSets/Flowering_Genes/gene_list_flowering_related_MolEcol_Athaliana_S4.csv"), sep=",", header=T, stringsAsFactors=FALSE, skip=1)
Athaliana.flowering.genes[,1] <- gsub(" ", "", Athaliana.flowering.genes[,1])
Athaliana.flowering.genes <- data.frame(AT_ID = Athaliana.flowering.genes[,1], Flowering = "TRUE")


### Drought responce genes 

UP.genes <- read.table(paste0(dataPath,"GeneControlSets/Drought_Genes/mDr_Day10_drought_up_regulated_genes_Harb_etal.csv"), sep=";", header=T, stringsAsFactors=FALSE)
DOWN.genes <- read.table(paste0(dataPath,"GeneControlSets/Drought_Genes/mDr_Day10_drought_down_regulated_genes_Harb_etal.csv"), sep=";", header=T, stringsAsFactors=FALSE)
Drought.genes <- rbind(UP.genes[,1:2], DOWN.genes[,1:2])
names(Drought.genes) <- c("AT_ID", "Drought_regulation")



genes.full.description <- merge(genes.description, Drought.genes, by="AT_ID" ,all.x=TRUE)
genes.full.description <- merge(genes.full.description, Athaliana.flowering.genes, by = "AT_ID", all.x = TRUE)

rownames(genes.full.description) <- genes.full.description$ID

write.table(genes.full.description, paste0(dataPath, "Data_2013-07-01/genes_description_full.xls"), sep="\t", quote=F, row.names=F)


save(new.samps, x, flowered.samps, control.samps, trees.order, all.days, month.days, all.days.short, month.days.short, genes.full.description ,file = "Shimizu_workspace.Rdata")












