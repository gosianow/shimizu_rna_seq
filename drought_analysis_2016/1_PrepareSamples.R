##############################################################################
## <<1_PrepareSamples.R>>

# BioC 3.3
# Created 22 June 2016

# Preparation of variables describing samples 

##############################################################################
Sys.time()
##############################################################################

library(stringr)
library(ggplot2)
library(reshape2)
library(limma)

##############################################################################
# Test arguments
##############################################################################

rwd='/home/Shared/data/seq/shimizu_rna_seq'
rcode='/home/gosia/R/shimizu_rna_seq/drought_analysis_2016'
data_dir='Data_2016'
analysis_dir='Analysis_2016-06-22'
out_dir='Plots_Samples'
Interpolate_and_roll_function='Interpolate_and_roll.R'

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

source(paste0(rcode, "/", Interpolate_and_roll_function))


####################################################
### Parse sample information
####################################################


### Read in sample_list
samps <- read.table(paste0(data_dir, "/Samples/orig/sample_list.txt"), sep="\t", stringsAsFactors=FALSE, header=TRUE)
head(samps)
names(samps)


### CORRECT THE DATE!!! (old dates are wrong, there is November 2009 and should be January 2009)
samps$year.month.day.old <- samps$year.month.day
samps$year.month.day <- paste(substr(samps$sample_name, nchar(samps$sample_name)-1, nchar(samps$sample_name)), substr(samps$sample_name, nchar(samps$sample_name)-3, nchar(samps$sample_name)-2), substr(samps$sample_name, nchar(samps$sample_name)-5, nchar(samps$sample_name)-4), sep=".")

samps[,c("sample_name","year.month.day", "year.month.day.old")]


### Columns 8-16 come from Lambir_meteorological_data_...xls; no gene variables 
new.samps <- samps[,c("sample_num", "sample_name","sample_ID","tree_ID","drough.control","developmental_stage","year.month.day")]

new.samps$flowered <- new.samps$tree_ID=="8266"

### define params for ploting: colors, legend
rownames(new.samps) <- new.samps$sample_name
new.samps$tree_ID <- as.factor(new.samps$tree_ID)
new.samps$tree_col <- new.samps$tree_ID
levels(new.samps$tree_col) <- c("orange", "cyan3", "green3", "blue", "red", "magenta3")
new.samps$tree_col <- as.character(new.samps$tree_col)
new.samps$tree_legend <- paste(new.samps$tree_ID, new.samps$drough.control, sep="-")
new.samps$short.name <- paste(new.samps$tree_ID, new.samps$year.month.day, sep="_")


### time format change 
new.samps$time_ch <- as.Date(new.samps$year.month.day, "%d.%m.%y")
new.samps$time_nr <- as.numeric(new.samps$time_ch)

write.table(new.samps, paste0(data_dir, "/Samples/samps.xls"), quote=FALSE, sep="\t", row.names=FALSE)


### Create an object with info to plot legends

trees.order <- unique(new.samps[, c("tree_legend", "tree_col", "tree_ID", "drough.control")])

rownames(trees.order) <- NULL
trees.order$condition <- ifelse(trees.order$drough.control == "control", "C", "D") 
trees.order <- trees.order[order(trees.order$condition, trees.order$tree_ID), ]

write.table(trees.order, paste0(data_dir, "/Samples/trees.xls"), quote=FALSE, sep="\t", row.names=FALSE)



### List of all unique days in 2008 and 2009 and dates of monts needed for nice x-axis labels

ad <- read.table(paste0(data_dir, "/Unique_days.csv"), sep=";")

all.days <- data.frame(days.ch = as.Date(ad[,1], "%d.%m.%y"), days.nr = as.numeric(as.Date(ad[,1], "%d.%m.%y")))

month.days <- all.days[strsplit2(all.days[, "days.ch"], "-")[, 3] == "01",]


# Shorter version

ad.short <- read.table(paste0(data_dir, "/Unique_days_short.csv"), sep=";")

all.days.short <- data.frame(days.ch = as.Date(ad.short[,1], "%d.%m.%y"), days.nr = as.numeric(as.Date(ad.short[,1], "%d.%m.%y")))

head(all.days.short)

write.table(all.days.short, paste0(data_dir, "/Unique_days_short_nr.xls"), quote=FALSE, sep="\t", row.names=FALSE)


month.days.short <- all.days.short[strsplit2(all.days.short[, "days.ch"], "-")[, 3] == "01",]

write.table(month.days.short, paste0(data_dir, "/Unique_days_short_nr_month.xls"), quote=FALSE, sep="\t", row.names=FALSE)



####################################################
### Rainfall & water deficit
####################################################


wd <- read.table(paste0(data_dir, "/MeteorologicalData/Lambir_meteorological_data_200710-201001mod.txt"), stringsAsFactors=FALSE, header=TRUE, sep = "\t")
head(wd)


wd$time_ch <- as.Date(wd$Date, "%d/%m/%y")
wd$time_nr <- as.numeric(wd$time_ch)

### Plot water deficit 

ggp <- ggplot(wd, aes(x = time_ch, y = evaporation...rain..mm.30.days.)) +
  geom_line(color = "brown", size = 2) +
  xlab("Time") +
  ylab("30-day moving total of water deficit (mm)") +
  ggtitle("Water Deficit") + 
  theme_bw() +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months")


pdf(paste0(out_dir, "/Water_deficit_raw.pdf"), width = 10, height = 5)
print(ggp)
dev.off()


ggp2 <- ggp +
  geom_vline(data = new.samps, aes(xintercept = time_nr), linetype = 3)

pdf(paste0(out_dir, "/Water_deficit_raw_experiment_dates.pdf"), width = 10, height = 5)
print(ggp2)
dev.off()


ggp3 <- ggp2 +
  scale_x_date(date_labels = "%b %d", date_breaks = "1 month", limits = c(min(all.days.short[,1]), max(all.days.short[,1])))


pdf(paste0(out_dir, "/Water_deficit_raw_experiment_dates_zoom.pdf"), width = 10, height = 5)
print(ggp3)
dev.off()



### Plot rainfall

ggp <- ggplot(wd, aes(x = time_ch, y = Rain.mm.30.days.)) +
  geom_line(color = "darkblue", size = 2) +
  xlab("Time") +
  ylab("30-day moving total of rainfall (mm)") +
  ggtitle("Rainfall") + 
  theme_bw() +
  scale_x_date(date_labels = "%b %Y",date_breaks = "3 months")


pdf(paste0(out_dir, "/Rainfall_raw.pdf"), width = 10, height = 5)
print(ggp)
dev.off()


ggp2 <- ggp +
  geom_vline(data = new.samps, aes(xintercept = time_nr), linetype = 3)

pdf(paste0(out_dir, "/Rainfall_raw_experiment_dates.pdf"), width = 10, height = 5)
print(ggp2)
dev.off()


ggp3 <- ggp2 +
  scale_x_date(date_labels = "%b %d", date_breaks = "1 month", limits = c(min(all.days.short[,1]), max(all.days.short[,1])))


pdf(paste0(out_dir, "/Rainfall_raw_experiment_dates_zoom.pdf"), width = 10, height = 5)
print(ggp3)
dev.off()



####################################################
### Soild Moisture
####################################################

### get soil moisture data
sm <- read.table(paste0(data_dir, "/MeteorologicalData/Soil_Moisture.csv"), sep=";", stringsAsFactors=FALSE, header=TRUE)
head(sm)


### time format change
sm$Year <- ifelse(sm$DayUniq <= 366, 2008, 2009)

sm$year.month.day <- paste(sm$Day,sm$Month, sm$Year, sep=".")
sm$time_ch <- as.Date(sm$year.month.day, "%d.%m.%Y")
sm$time_nr <- as.numeric(sm$time_ch)


### reshape sm data

smm <- melt(sm, id.vars = c("Month", "Day", "DayUniq", "Year", "year.month.day", "time_ch", "time_nr"), value.name = "soil_moisture", variable.name = "condition_tree", factorsAsStrings = FALSE)

### fix "time_ch"
smm$time_ch <- as.Date(smm$year.month.day, "%d.%m.%Y")

smm$tree <- factor(gsub("[[:alpha:]]", "", smm$condition_tree), levels = trees.order$tree_ID[trees.order$drough.control == "drought"])
smm$condition <- factor(substr(smm$condition_tree, start = 1, stop = 1), levels = c("D", "C"))


ggp <- ggplot(smm, aes(x = time_ch, y = soil_moisture, group = condition_tree, color = tree, linetype = condition)) +
  geom_line(size = 1) +
  xlab("Time") +
  ylab("Soil moisture") +
  ggtitle("Soil Moisture") + 
  theme_bw() +
  scale_x_date(date_labels = "%b %Y", date_breaks = "2 months") +
  scale_color_manual(values = trees.order$tree_col[trees.order$drough.control == "drought"])


pdf(paste0(out_dir, "/Soil_moisture.pdf"), width = 12, height = 5)
print(ggp)
dev.off()



new.samps_drought <- new.samps[new.samps$drough.control == "drought", ]
new.samps_drought$tree <- factor(new.samps_drought$tree_ID, levels = trees.order$tree_ID[trees.order$drough.control == "drought"])



ggp2 <- ggp +
  geom_vline(data = new.samps_drought, aes(xintercept = time_nr), linetype = 3) +
  facet_wrap(~ tree, ncol = 1)
  
  
pdf(paste0(out_dir, "/Soil_moisture_per_tree.pdf"), width = 12, height = 15)
print(ggp2)
dev.off()


### Interpolate and roll - add soil moisture info to new.samps

roll.days=c(14, 21, 28)

new.samps.DE970 <- interpolate.and.roll(intrp.x=sm$time_nr, intrp.y=sm$DE970, intrp.points=all.days[,2], table1=new.samps[new.samps$tree_ID == 970,], roll.days=roll.days, col.values="Soil.Moisture", plots.path = paste0(out_dir, "/Plots_Samples_Interpolate"), plot.name="Soil_MoistureDE970", color=trees.order$tree_col[trees.order$tree_ID=="970"], month.days=month.days, all.days.short=all.days.short)

new.samps.DE8212 <- interpolate.and.roll(intrp.x=sm$time_nr, intrp.y=sm$DE8212, intrp.points=all.days[,2], table1=new.samps[new.samps$tree_ID == 8212,], roll.days=roll.days, col.values="Soil.Moisture", plots.path = paste0(out_dir, "/Plots_Samples_Interpolate"),plot.name="Soil_MoistureDE8212", color=trees.order$tree_col[trees.order$tree_ID=="8212"], month.days=month.days, all.days.short=all.days.short)

new.samps.DE8266 <- interpolate.and.roll(intrp.x=sm$time_nr, intrp.y=sm$DE8266, intrp.points=all.days[,2], table1=new.samps[new.samps$tree_ID == 8266,], roll.days=roll.days, col.values="Soil.Moisture", plots.path = paste0(out_dir, "/Plots_Samples_Interpolate"), plot.name="Soil_MoistureDE8266",  color=trees.order$tree_col[trees.order$tree_ID=="8266"], month.days=month.days, all.days.short=all.days.short)


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

wp <- read.table(paste0(data_dir, "/MeteorologicalData/Water_potential_from_Inoue.csv"), sep=";", header=TRUE, as.is = TRUE)
head(wp)
names(wp)


### make wp values in (0,1) except SD
wp[,-c(1, grep("SD", colnames(wp)))] <- wp[,-c(1, grep("SD", colnames(wp)))] + 1

### change time format
wp$time_ch <- as.Date(wp$Date, "%y.%m.%d")
wp$time_nr <- as.numeric(wp$time_ch)


### reshape wp data

wpm <- wp[, grep("Date|sMean|time", colnames(wp))]

wpm <- melt(wpm, id.vars = c("Date", "time_ch", "time_nr"), value.name = "water_potential", variable.name = "condition_tree", factorsAsStrings = FALSE)

### fix "time_ch"
wpm$time_ch <- as.Date(wpm$Date, "%y.%m.%d")

wpm$tree <- factor(gsub("[[:alpha:]]", "", wpm$condition_tree), levels = trees.order$tree_ID)
wpm$condition <- factor(substr(wpm$condition_tree, start = 1, stop = 1), levels = c("D", "C"))
wpm <- wpm[complete.cases(wpm$water_potential), ]


ggp <- ggplot(wpm, aes(x = time_ch, y = water_potential, group = tree, color = tree)) +
  geom_line(size = 1) +
  xlab("Time") +
  ylab("Mean Water Potential") +
  ggtitle("Water Potential") + 
  theme_bw() +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  scale_color_manual(values = trees.order$tree_col)


pdf(paste0(out_dir, "/Water_potential.pdf"), width = 12, height = 5)
print(ggp)
dev.off()


new.samps$tree <- factor(new.samps$tree_ID, levels = trees.order$tree_ID)


ggp2 <- ggp +
  geom_vline(data = new.samps, aes(xintercept = time_nr), linetype = 3) +
  facet_wrap(~ tree, ncol = 1)


pdf(paste0(out_dir, "/Water_potential_per_tree.pdf"), width = 12, height = 30)
print(ggp2)
dev.off()





### add avg values to the May samples

wp <- wp[order(wp$time_nr), ]
new.samps <- new.samps[order(new.samps$time_nr), ]

tail(wp)

wp.temp <- as.data.frame(matrix(NA, 2, ncol(wp)))
colnames(wp.temp) <- colnames(wp)

wp.temp$time_nr <- c(1241992800, 1242079200)
wp.temp$time_ch <- as.Date(c("2009-05-11", "2009-05-12"), "%Y-%m-%d")

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

tree.ids <- as.character(trees.order$tree_ID)


### intepolation
for(i in 1:length(tree.ids)){
  # i=1
  tree.sampl <- grep(paste0(tree.ids[i], "sMean"), colnames(wp))
  
  new.samps.wp.tmp <- interpolate.and.roll(intrp.x=wp$time_nr[!is.na(wp[,tree.sampl])], intrp.y=wp[!is.na(wp[,tree.sampl]),tree.sampl], intrp.points=all.days[,2], table1=new.samps[new.samps$tree_ID == tree.ids[i],] , roll.days=c(14, 21, 28), col.values="Water.Potential", plots.path = paste0(out_dir, "/Plots_Samples_Interpolate"), plot.name=paste0("Water_Potential", tree.ids[i]), month.days=month.days, all.days.short=all.days.short, color=trees.order$tree_col[trees.order$tree_ID==tree.ids[i]])
  
  new.samps.wp <- rbind(new.samps.wp, new.samps.wp.tmp)
  
}

names(new.samps.wp)
head(new.samps.wp)


new.samps <- new.samps.wp


####################################################
### Temperature
####################################################

tempr <- read.table(paste0(data_dir, "/MeteorologicalData/Lambir_temperature_data.csv"), head=T, sep=";")
head(tempr)


tempr$time_ch <- as.Date(substring(tempr$Time, 1, 8), "%d.%m.%y")
tempr$time_nr <- as.numeric(tempr$time_ch)

tempr$time_hour <- substring(tempr$Time, 10, 11)
tempr$time_hour_nr <- as.numeric(tempr$time_hour)

pdf(paste0(out_dir, "/Temperature_raw.pdf"), width = 10, height = 5)
plot(tempr$time_nr, tempr$Temp, pch=".", main="Temperature", xlab="Time", ylab="Temperature")
#abline(v=new.samps$time_nr, col="grey")
dev.off()

pdf(paste0(out_dir, "/Temperature_over_day.pdf"), width = 10, height = 5)
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

pdf(paste0(out_dir, "/Temperature_Avg_lowess.pdf"), width = 10, height = 5)
plot(tempr.day[,"time_nr"], tempr.day[,"TempAvg"], col="grey", pch=20, main="Control trees", ylab="Soil moisture", xlab="Time",  xaxt = "n")
axis(side=1, at=month.days[,2], labels=month.days[,1])
lines(tempr.l[["TempAvg"]], lwd=2, col="darkred")
dev.off()


### some plots of avg, max, min temperature
pdf(paste0(out_dir, "/Temperature_Avg.pdf"), width = 10, height = 5)
plot(tempr.day$time_nr, tempr.day$TempAvg, pch=20, main="Average Temperature", xlab="Time", ylab="Temperature",xlim=c(min(all.days.short[,2]), max(all.days.short[,2])),  xaxt = "n")
axis(side=1, at=month.days[,2], labels=month.days[,1])
lines(tempr.l[["TempAvg"]], lwd=2, col="darkred")
abline(v=new.samps$time_nr, col="grey")
dev.off()


tempr.l[["TempMin"]] <- lowess(tempr.day[,"time_nr"], tempr.day[,"TempMin"], f=0.1)

pdf(paste0(out_dir, "/Temperature_Min_lowess.pdf"), width = 10, height = 5)
plot(tempr.day[,"time_nr"], tempr.day[,"TempMin"], col="grey", pch=20, main="Control trees", ylab="Soil moisture", xlab="Time",  xaxt = "n")
axis(side=1, at=month.days[,2], labels=month.days[,1])
lines(tempr.l[["TempMin"]], lwd=2, col="darkred")
dev.off()


### some plots of avg, max, min temperature
pdf(paste0(out_dir, "/Temperature_Min.pdf"), width = 10, height = 5)
plot(tempr.day$time_nr, tempr.day$TempMin, pch=20, main="Min Temperature", xlab="Time", ylab="Temperature",xlim=c(min(all.days.short[,2]), max(all.days.short[,2])),  xaxt = "n")
axis(side=1, at=month.days[,2], labels=month.days[,1])
lines(tempr.l[["TempMin"]], lwd=2, col="darkred")
abline(v=new.samps$time_nr, col="grey")
dev.off()



tempr.ll <- data.frame(time_nr = tempr.l$TempAvg$x, TempAvg = tempr.l$TempAvg$y, TempMin = tempr.l$TempMin$y)
head(tempr.ll)

### rolled means of temperature

roll.days = c(14, 21, 28)
variables = c("TempAvg", "TempMin")
plots.path = paste0(out_dir, "/Plots_Samples_Interpolate")

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


####################################################
### Save new.samps
####################################################


write.table(new.samps, paste0(data_dir, "/Samples/new_samps_interpolation.xls"), quote=FALSE, sep="\t", row.names=FALSE)





























