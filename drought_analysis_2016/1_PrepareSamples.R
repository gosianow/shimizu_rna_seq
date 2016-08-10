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

### Clean up column names 
colnames(samps) <- tolower(gsub("\\.", "_", colnames(samps)))


### CORRECT THE DATE!!! (old dates are wrong, there is November 2009 and should be January 2009)
samps$year_month_day_old <- samps$year_month_day
samps$day_month_year <- paste(substr(samps$sample_name, nchar(samps$sample_name)-1, nchar(samps$sample_name)), substr(samps$sample_name, nchar(samps$sample_name)-3, nchar(samps$sample_name)-2), substr(samps$sample_name, nchar(samps$sample_name)-5, nchar(samps$sample_name)-4), sep=".")



samps[,c("sample_name","day_month_year", "year_month_day_old")]


### Columns 8-16 come from Lambir_meteorological_data_...xls; no gene variables 
samps_new <- samps[,c("sample_num", "sample_name","sample_id","tree_id","drough_control","developmental_stage","day_month_year")]

samps_new$flowered <- samps_new$tree_id=="8266"

### define params for ploting: colors, legend
rownames(samps_new) <- samps_new$sample_name
samps_new$tree_id <- as.factor(samps_new$tree_id)
samps_new$tree_col <- samps_new$tree_id
levels(samps_new$tree_col) <- c("orange", "cyan3", "green3", "blue", "red", "magenta3")
samps_new$tree_col <- as.character(samps_new$tree_col)
samps_new$tree_legend <- paste(samps_new$tree_id, samps_new$drough_control, sep="-")
samps_new$sample_name_short <- paste(samps_new$tree_id, samps_new$day_month_year, sep="_")


### time format change 
samps_new$time_ch <- as.Date(samps_new$day_month_year, "%d.%m.%y")
samps_new$time_nr <- as.numeric(samps_new$time_ch)

write.table(samps_new, paste0(data_dir, "/Samples/samps.xls"), quote=FALSE, sep="\t", row.names=FALSE)



####################################################
### Create an object with info to plot legends
####################################################


trees_order <- unique(samps_new[, c("tree_legend", "tree_col", "tree_id", "drough_control")])

rownames(trees_order) <- NULL
trees_order$condition <- ifelse(trees_order$drough_control == "control", "C", "D") 
trees_order <- trees_order[order(trees_order$condition, trees_order$tree_id), ]

write.table(trees_order, paste0(data_dir, "/Samples/trees.xls"), quote=FALSE, sep="\t", row.names=FALSE)




####################################################
### Plot the time that sampels were taken for each tree
####################################################

ggdf <- samps_new[samps_new$developmental_stage == "leaf_bud", c("time_ch", "tree_id")]

ggdf <- data.frame(table(ggdf$time_ch, ggdf$tree_id))

ggdf <- ggdf[ggdf$Freq > 0, ]

ggdf$time_ch <- as.Date(ggdf$Var1, "%Y-%m-%d")
ggdf$time_nr <- as.numeric(ggdf$time_ch)


ggdf$tree_id <- factor(ggdf$Var2, levels = trees_order$tree_id)

ggdf$Freq <- factor(ggdf$Freq)


ggp <- ggplot(ggdf, aes(x = time_ch, y = tree_id, color = tree_id, shape = Freq)) +
  geom_vline(aes(xintercept = time_nr), linetype = 3) +
  geom_point(size = 3) +
  xlab("Time") +
  ylab("Tree") +
  ggtitle("Time when the samples were taken") + 
  theme_bw() +
  scale_x_date(date_labels = "%d %b %Y", date_breaks = "1 month") +
  scale_color_manual(values = trees_order$tree_col) +
  scale_shape_manual(values = c(17, 18)) 
  



pdf(paste0(out_dir, "/Timing.pdf"), width = 10, height = 5)
print(ggp)
dev.off()




####################################################
### List of all unique days in 2008 and 2009 and dates of monts needed for nice x-axis labels
####################################################

ad <- read.table(paste0(data_dir, "/Unique_days.csv"), sep=";")

all.days <- data.frame(days_ch = as.Date(ad[,1], "%d.%m.%y"), days_nr = as.numeric(as.Date(ad[,1], "%d.%m.%y")))

month.days <- all.days[strsplit2(all.days[, "days_ch"], "-")[, 3] == "01",]


# Shorter version

ad.short <- read.table(paste0(data_dir, "/Unique_days_short.csv"), sep=";")

all.days.short <- data.frame(days_ch = as.Date(ad.short[,1], "%d.%m.%y"), days_nr = as.numeric(as.Date(ad.short[,1], "%d.%m.%y")))

head(all.days.short)

write.table(all.days.short, paste0(data_dir, "/Unique_days_short_nr.xls"), quote=FALSE, sep="\t", row.names=FALSE)


month.days.short <- all.days.short[strsplit2(all.days.short[, "days_ch"], "-")[, 3] == "01",]

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
  scale_x_date(date_labels = "%d %b %Y", date_breaks = "3 months")


pdf(paste0(out_dir, "/Water_deficit_raw.pdf"), width = 10, height = 5)
print(ggp)
dev.off()


ggp2 <- ggp +
  geom_vline(data = samps_new, aes(xintercept = time_nr), linetype = 3)

pdf(paste0(out_dir, "/Water_deficit_raw_experiment_dates.pdf"), width = 10, height = 5)
print(ggp2)
dev.off()


ggp3 <- ggp2 +
  scale_x_date(date_labels = "%d %b %Y", date_breaks = "1 month", limits = c(min(all.days.short[,1]), max(all.days.short[,1])))


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
  scale_x_date(date_labels = "%d %b %Y",date_breaks = "3 months")


pdf(paste0(out_dir, "/Rainfall_raw.pdf"), width = 10, height = 5)
print(ggp)
dev.off()


ggp2 <- ggp +
  geom_vline(data = samps_new, aes(xintercept = time_nr), linetype = 3)

pdf(paste0(out_dir, "/Rainfall_raw_experiment_dates.pdf"), width = 10, height = 5)
print(ggp2)
dev.off()


ggp3 <- ggp2 +
  scale_x_date(date_labels = "%d %b %Y", date_breaks = "1 month", limits = c(min(all.days.short[,1]), max(all.days.short[,1])))


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

sm$day_month_year <- paste(sm$Day,sm$Month, sm$Year, sep=".")
sm$time_ch <- as.Date(sm$day_month_year, "%d.%m.%Y")
sm$time_nr <- as.numeric(sm$time_ch)


### reshape sm data

smm <- melt(sm, id.vars = c("Month", "Day", "DayUniq", "Year", "day_month_year", "time_ch", "time_nr"), value.name = "soil_moisture", variable.name = "condition_tree", factorsAsStrings = FALSE)

### fix "time_ch"
smm$time_ch <- as.Date(smm$day_month_year, "%d.%m.%Y")

smm$tree <- factor(gsub("[[:alpha:]]", "", smm$condition_tree), levels = trees_order$tree_id[trees_order$drough_control == "drought"])
smm$condition <- factor(substr(smm$condition_tree, start = 1, stop = 1), levels = c("D", "C"))


ggp <- ggplot(smm, aes(x = time_ch, y = soil_moisture, group = condition_tree, color = tree, linetype = condition)) +
  geom_line(size = 1) +
  xlab("Time") +
  ylab("Soil moisture") +
  ggtitle("Soil Moisture") + 
  theme_bw() +
  scale_x_date(date_labels = "%d %b %Y", date_breaks = "2 months") +
  scale_color_manual(values = trees_order$tree_col[trees_order$drough_control == "drought"])


pdf(paste0(out_dir, "/Soil_moisture.pdf"), width = 12, height = 5)
print(ggp)
dev.off()



samps_new_drought <- samps_new[samps_new$drough_control == "drought", ]
samps_new_drought$tree <- factor(samps_new_drought$tree_id, levels = trees_order$tree_id[trees_order$drough_control == "drought"])



ggp2 <- ggp +
  geom_vline(data = samps_new_drought, aes(xintercept = time_nr), linetype = 3) +
  facet_wrap(~ tree, ncol = 1)
  
  
pdf(paste0(out_dir, "/Soil_moisture_per_tree.pdf"), width = 12, height = 15)
print(ggp2)
dev.off()


### Interpolate and roll - add soil moisture info to samps_new

roll.days=c(14, 21, 28)

samps_new.DE970 <- interpolate.and.roll(intrp.x=sm$time_nr, intrp.y=sm$DE970, intrp.points=all.days[,2], table1=samps_new[samps_new$tree_id == 970,], roll.days=roll.days, col.values="Soil.Moisture", plots.path = paste0(out_dir, "/Plots_Samples_Interpolate"), plot.name="Soil_MoistureDE970", color=trees_order$tree_col[trees_order$tree_id=="970"], month.days=month.days, all.days.short=all.days.short)

samps_new.DE8212 <- interpolate.and.roll(intrp.x=sm$time_nr, intrp.y=sm$DE8212, intrp.points=all.days[,2], table1=samps_new[samps_new$tree_id == 8212,], roll.days=roll.days, col.values="Soil.Moisture", plots.path = paste0(out_dir, "/Plots_Samples_Interpolate"),plot.name="Soil_MoistureDE8212", color=trees_order$tree_col[trees_order$tree_id=="8212"], month.days=month.days, all.days.short=all.days.short)

samps_new.DE8266 <- interpolate.and.roll(intrp.x=sm$time_nr, intrp.y=sm$DE8266, intrp.points=all.days[,2], table1=samps_new[samps_new$tree_id == 8266,], roll.days=roll.days, col.values="Soil.Moisture", plots.path = paste0(out_dir, "/Plots_Samples_Interpolate"), plot.name="Soil_MoistureDE8266",  color=trees_order$tree_col[trees_order$tree_id=="8266"], month.days=month.days, all.days.short=all.days.short)


samps_new.control <- samps_new[samps_new$drough_control == "control",]
Soil.moisture.control <- data.frame(time_nr = samps_new.control$time_nr, matrix(NA, nrow = nrow(samps_new.control), ncol = (length(roll.days) + 1) ))
colnames(Soil.moisture.control) <- c("time_nr", paste0("Soil.Moisture", c("", roll.days)))

samps_new.control <- merge(samps_new.control, Soil.moisture.control, by="time_nr")

samps_new.sm <- rbind(samps_new.DE970, samps_new.DE8212, samps_new.DE8266, samps_new.control)

samps_new <- samps_new.sm

names(samps_new)


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

wpm$tree <- factor(gsub("[[:alpha:]]", "", wpm$condition_tree), levels = trees_order$tree_id)
wpm$condition <- factor(substr(wpm$condition_tree, start = 1, stop = 1), levels = c("D", "C"))
wpm <- wpm[complete.cases(wpm$water_potential), ]


ggp <- ggplot(wpm, aes(x = time_ch, y = water_potential, group = tree, color = tree)) +
  geom_line(size = 1) +
  xlab("Time") +
  ylab("Mean Water Potential") +
  ggtitle("Water Potential") + 
  theme_bw() +
  scale_x_date(date_labels = "%d %b %Y", date_breaks = "1 month") +
  scale_color_manual(values = trees_order$tree_col)


pdf(paste0(out_dir, "/Water_potential.pdf"), width = 12, height = 5)
print(ggp)
dev.off()


samps_new$tree <- factor(samps_new$tree_id, levels = trees_order$tree_id)


ggp2 <- ggp +
  geom_vline(data = samps_new, aes(xintercept = time_nr), linetype = 3) +
  facet_wrap(~ tree, ncol = 1)


pdf(paste0(out_dir, "/Water_potential_per_tree.pdf"), width = 12, height = 30)
print(ggp2)
dev.off()





### add avg values to the May samples

wp <- wp[order(wp$time_nr), ]
samps_new <- samps_new[order(samps_new$time_nr), ]

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
samps_new.wp <- NULL

tree.ids <- as.character(trees_order$tree_id)


### intepolation
for(i in 1:length(tree.ids)){
  # i=1
  tree.sampl <- grep(paste0(tree.ids[i], "sMean"), colnames(wp))
  
  samps_new.wp.tmp <- interpolate.and.roll(intrp.x=wp$time_nr[!is.na(wp[,tree.sampl])], intrp.y=wp[!is.na(wp[,tree.sampl]),tree.sampl], intrp.points=all.days[,2], table1=samps_new[samps_new$tree_id == tree.ids[i],] , roll.days=c(14, 21, 28), col.values="Water.Potential", plots.path = paste0(out_dir, "/Plots_Samples_Interpolate"), plot.name=paste0("Water_Potential", tree.ids[i]), month.days=month.days, all.days.short=all.days.short, color=trees_order$tree_col[trees_order$tree_id==tree.ids[i]])
  
  samps_new.wp <- rbind(samps_new.wp, samps_new.wp.tmp)
  
}

names(samps_new.wp)
head(samps_new.wp)


samps_new <- samps_new.wp


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
#abline(v=samps_new$time_nr, col="grey")
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
abline(v=samps_new$time_nr, col="grey")
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
abline(v=samps_new$time_nr, col="grey")
dev.off()



tempr.ll <- data.frame(time_nr = tempr.l$TempAvg$x, TempAvg = tempr.l$TempAvg$y, TempMin = tempr.l$TempMin$y)
head(tempr.ll)

### rolled means of temperature

roll.days = c(14, 21, 28)
variables = c("TempAvg", "TempMin")
plots.path = paste0(out_dir, "/Plots_Samples_Interpolate")

samps_new <- merge(samps_new, tempr.ll, by="time_nr", all.x = TRUE)

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
    abline(v=samps_new$time_nr, col="grey")

    lines(tempr.l[[v]], lwd=1, col="darkred", lty=3)

      lines(tempr.day$time_nr, tempr.day[,name], col="darkred", lty=1, lwd=4)

  }
  
  
  dev.off()
  
  samps_new <- merge(samps_new, tempr.day[,c("time_nr", paste(v, roll.days, sep=""))], by="time_nr")
}

head(samps_new)



samps_new <- unique(samps_new)


####################################################
### Save samps_new
####################################################


### Clean up column names 
colnames(samps_new) <- tolower(gsub("\\.", "_", colnames(samps_new)))


write.table(samps_new, paste0(data_dir, "/Samples/new_samps_interpolation.xls"), quote=FALSE, sep="\t", row.names=FALSE)


####################################################
### Create a time grouping variable
####################################################


samps_new <- read.table(paste0(data_dir, "/Samples/new_samps_interpolation.xls"), header = TRUE, sep = "\t", as.is = TRUE)


samps_new$time_ch <- as.Date(samps_new$time_ch, "%Y-%m-%d")

samps_new$time_group <- NA

samps_new$time_group[samps_new$time_ch <= as.Date("2008-11-01", "%Y-%m-%d")] <- "time1"

samps_new$time_group[samps_new$time_ch > as.Date("2008-11-01", "%Y-%m-%d") & samps_new$time_ch <= as.Date("2008-12-03", "%Y-%m-%d")] <- "time2"

samps_new$time_group[samps_new$time_ch > as.Date("2008-12-03", "%Y-%m-%d") & samps_new$time_ch <= as.Date("2008-12-09", "%Y-%m-%d")] <- "time3"

samps_new$time_group[samps_new$time_ch > as.Date("2008-12-09", "%Y-%m-%d") & samps_new$time_ch <= as.Date("2008-12-22", "%Y-%m-%d")] <- "time4"

samps_new$time_group[samps_new$time_ch > as.Date("2008-12-22", "%Y-%m-%d") & samps_new$time_ch <= as.Date("2009-01-05", "%Y-%m-%d")] <- "time5"

samps_new$time_group[samps_new$time_ch > as.Date("2009-01-05", "%Y-%m-%d") & samps_new$time_ch <= as.Date("2009-01-19", "%Y-%m-%d")] <- "time6"


samps_new$time_group[samps_new$time_ch > as.Date("2009-01-19", "%Y-%m-%d") & samps_new$time_ch <= as.Date("2009-02-12", "%Y-%m-%d")] <- "time7"

samps_new$time_group[samps_new$time_ch > as.Date("2009-02-12", "%Y-%m-%d") & samps_new$time_ch <= as.Date("2009-02-24", "%Y-%m-%d")] <- "time8"

samps_new$time_group[samps_new$time_ch > as.Date("2009-02-24", "%Y-%m-%d") & samps_new$time_ch <= as.Date("2009-03-08", "%Y-%m-%d")] <- "time9"


samps_new$time_group[samps_new$time_ch > as.Date("2009-03-08", "%Y-%m-%d") & samps_new$time_ch <= as.Date("2009-04-16", "%Y-%m-%d")] <- "time10"


samps_new$time_group[samps_new$time_ch > as.Date("2009-04-16", "%Y-%m-%d") & samps_new$time_ch <= as.Date("2009-05-12", "%Y-%m-%d")] <- "time11"


samps_new$time_group_last_day <- NA

samps_new$time_group_last_day[samps_new$time_ch <= as.Date("2008-11-01", "%Y-%m-%d")] <- "2008-11-01"

samps_new$time_group_last_day[samps_new$time_ch > as.Date("2008-11-01", "%Y-%m-%d") & samps_new$time_ch <= as.Date("2008-12-03", "%Y-%m-%d")] <- "2008-12-03"

samps_new$time_group_last_day[samps_new$time_ch > as.Date("2008-12-03", "%Y-%m-%d") & samps_new$time_ch <= as.Date("2008-12-09", "%Y-%m-%d")] <- "2008-12-09"

samps_new$time_group_last_day[samps_new$time_ch > as.Date("2008-12-09", "%Y-%m-%d") & samps_new$time_ch <= as.Date("2008-12-22", "%Y-%m-%d")] <- "2008-12-22"

samps_new$time_group_last_day[samps_new$time_ch > as.Date("2008-12-22", "%Y-%m-%d") & samps_new$time_ch <= as.Date("2009-01-05", "%Y-%m-%d")] <- "2009-01-05"

samps_new$time_group_last_day[samps_new$time_ch > as.Date("2009-01-05", "%Y-%m-%d") & samps_new$time_ch <= as.Date("2009-01-19", "%Y-%m-%d")] <- "2009-01-19"


samps_new$time_group_last_day[samps_new$time_ch > as.Date("2009-01-19", "%Y-%m-%d") & samps_new$time_ch <= as.Date("2009-02-12", "%Y-%m-%d")] <- "2009-02-12"

samps_new$time_group_last_day[samps_new$time_ch > as.Date("2009-02-12", "%Y-%m-%d") & samps_new$time_ch <= as.Date("2009-02-24", "%Y-%m-%d")] <- "2009-02-24"

samps_new$time_group_last_day[samps_new$time_ch > as.Date("2009-02-24", "%Y-%m-%d") & samps_new$time_ch <= as.Date("2009-03-08", "%Y-%m-%d")] <- "2009-03-08"


samps_new$time_group_last_day[samps_new$time_ch > as.Date("2009-03-08", "%Y-%m-%d") & samps_new$time_ch <= as.Date("2009-04-16", "%Y-%m-%d")] <- "2009-04-16"


samps_new$time_group_last_day[samps_new$time_ch > as.Date("2009-04-16", "%Y-%m-%d") & samps_new$time_ch <= as.Date("2009-05-12", "%Y-%m-%d")] <- "2009-05-12"


write.table(samps_new, paste0(data_dir, "/Samples/new_samps_interpolation.xls"), quote=FALSE, sep="\t", row.names=FALSE)

####################################################
### Plot the time that sampels were taken for each tree
####################################################

ggdf <- samps_new[samps_new$developmental_stage == "leaf_bud", c("time_ch", "tree_id")]

ggdf <- data.frame(table(ggdf$time_ch, ggdf$tree_id))

ggdf <- ggdf[ggdf$Freq > 0, ]

ggdf$time_ch <- as.Date(ggdf$Var1, "%Y-%m-%d")
ggdf$time_nr <- as.numeric(ggdf$time_ch)

ggdf$tree_id <- factor(ggdf$Var2, levels = trees_order$tree_id)

ggdf$Freq <- factor(ggdf$Freq)



samps_new$time_group_last_day <- as.Date(samps_new$time_group_last_day, "%Y-%m-%d")
samps_new$time_group_last_day_nr <- as.numeric(samps_new$time_group_last_day + 1)



ggp <- ggplot(ggdf, aes(x = time_ch, y = tree_id, color = tree_id, shape = Freq)) +
  geom_vline(aes(xintercept = time_nr), linetype = 3) +
  geom_vline(data = samps_new, aes(xintercept = time_group_last_day_nr), linetype = 1, size = 0.1, color = 2) +
  geom_point(size = 3) +
  xlab("Time") +
  ylab("Tree") +
  ggtitle("Time when the samples were taken") + 
  theme_bw() +
  scale_x_date(date_labels = "%d %b %Y", date_breaks = "1 month") +
  scale_color_manual(values = trees_order$tree_col) +
  scale_shape_manual(values = c(17, 18)) 




pdf(paste0(out_dir, "/Timing_groups.pdf"), width = 10, height = 5)
print(ggp)
dev.off()





























