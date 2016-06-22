
#####################################################################################################

### preparation of variables describing samples - Soild Moisture & Water Potential & Temperature
### Model fitting : run.egdeR(); different combinations of models; comparison with UP & DOWN regulated genes

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

#new.samps <- samps[, c(1:16)]

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
  
  pdf( paste(plots.path, "/Interp2_", str_replace_all(paste(col.values, plot.name, sep=""),"[[:punct:]]", "_"),".pdf", sep=""), width = 10, height = 5)
  plot(intrp.x, intrp.y, col=1, pch=20, ylab=col.values, xlab="Time")
  abline(v=table1[,"time_nr"], col="grey")
  lines(interp.intrp.points, col=6)
  points(table1[,"time_nr"], table1[,col.values], col=4)
  dev.off()
  
  pdf( paste(plots.path, "/InterpRoll_", str_replace_all(paste(col.values, plot.name, sep=""),"[[:punct:]]", "_"),".pdf", sep=""), width = 10, height = 5)
  plot(intrp.x, intrp.y, col=1, pch=20, ylab=col.values, xlab="Time")
  abline(v=table1[,"time_nr"], col="grey")
  lines(table2[,"time_nr"], table2[,col.values], col=6)
  lines(table2[,"time_nr"], table2[,col.values14], col=6,  lty=2)
  lines(table2[,"time_nr"], table2[,col.values28], col=6, lty=3)
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

# p = wp[,"DE8212sMean"]
# plot(p, log(p/(1-p)))


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
  
  pdf(paste("Plots_Raw/Water_potential_", exp.tree.ids[i], ".pdf", sep=""), width = 10, height = 5)
  plot(wp$time_nr[!is.na(wp[,tree.sampl])], wp[!is.na(wp[,tree.sampl]), tree.sampl], type="l", col=3, ylim=c(0,1), xlim=c(1.220e+09, 1257289200), main=exp.tree.ids[i], ylab="Mean Water Potential", xlab="Time")
  abline(v=new.samps$time_nr[new.samps$tree_ID==tree.ids[i]], col="grey")
  dev.off()
  
  pdf(paste("Plots_Raw/Water_potential_", exp.tree.ids[i], "_interpolation.pdf", sep=""), width = 10, height = 5)
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

#####################################################################################################

### GLM // gene expression vs Soil Moisture & Water Potential

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


AT.id <- read.table("Data/genes_descr_control/best_hit_blast_result_Sl_predicted_exons.txt", sep=",", stringsAsFactors=FALSE)
head(AT.id)

UP.genes <- read.table("Data/genes_descr_control/mDr_Day10_drought_up_regulated_genes_Harb_etal.csv", sep=";", header=T, stringsAsFactors=FALSE)
DOWN.genes <- read.table("Data/genes_descr_control/mDr_Day10_drought_down_regulated_genes_Harb_etal.csv", sep=";", header=T, stringsAsFactors=FALSE)

# Drought.genes <- read.table("Data/genes_descr_control/mDr_Day10_drought_regulated_genes_Harb_etal.csv", sep=";", header=T, stringsAsFactors=FALSE)


MolEcol.genes.c1 <- read.table("Data/genes_descr_control/gene_list_in_clusterI_Mol_Ecol.csv", sep=";", header=T, stringsAsFactors=FALSE)
MolEcol.genes.c2 <- read.table("Data/genes_descr_control/gene_list_in_clusterII_Mol_Ecol.csv", sep=";", header=T, stringsAsFactors=FALSE)

# MolEcol.DE.genes <- read.table("Data/genes_descr_control/gene_list_in_all_clusters_DEgenes_Mol_Ecol.csv", sep=";", header=T, stringsAsFactors=FALSE)


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


run.edgeR <- function(x, new.samps, model.formula, varialbs, elim.samps, n.top=200, plots.path="PlotsTest", plot.name="", AT.id, UP.genes, DOWN.genes, varialbs.plot, fit.clr = 3, LRTcoef=2)
  { 
  
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
    
    plot(Top.fitting.list[[i]][,1], Top.fitting.list[[i]][,2], main=paste(i, "/", AT, "\n", "Coefffs:",paste(coeffs, collapse=", "), "\n FDR:", top.tags$table[i, "FDR"], "\n UP:", if.UP, "DOWN:", if.DOWN), xlab = varialbs.plot, ylab="Log(Counts in CPM)", cex.main=0.7, pch=ifelse(new.samps$drough.control=="drought", 5, 1))
    points(Top.fitting.list[[i]][,1], Top.fitting.list[[i]][,3], col=fit.clr)
    
  }
  
  dev.off()
  
  #save(Top.fitting.list, file=paste(plots.path, "/Top_fitting_", model.char, "_cpm", plot.name,".RData", sep=""))
  
  UP.coeffs <- UP.coeffs[complete.cases(UP.coeffs), ]  
  DOWN.coeffs <- DOWN.coeffs[complete.cases(DOWN.coeffs), ]
  
  invisible(list(Top.genes=top.tags$table, Top.fitting.list=Top.fitting.list, UP.coeffs=UP.coeffs, DOWN.coeffs=DOWN.coeffs))
  
}


### RUN run.egdeR()


model.var <- c("Soil.Moisture28", "Water.Potential28", "Temperature28")


for (i in model.var){
# i="Soil.Moisture28"
run.edgeR(x, new.samps, model.formula = as.formula(paste("~", i)), varialbs=c(i), elim.samps=flowered.samps, n.top=200, plots.path="PlotsRUNedgeR_22Aug", plot.name="", AT.id, UP.genes, DOWN.genes,  varialbs.plot = i, fit.clr = 4, LRTcoef=2)

run.edgeR(x, new.samps, model.formula = as.formula(paste("~", i)), varialbs=c(i), elim.samps=c(control.samps, flowered.samps), n.top=200, plots.path="PlotsRUNedgeR_22Aug", plot.name="NoContrl", AT.id, UP.genes, DOWN.genes,  varialbs.plot = i, fit.clr = 4, LRTcoef=2)


run.edgeR(x, new.samps, model.formula = as.formula(paste("~", i , "+ drough.control")) , varialbs=c(i), elim.samps=c(flowered.samps), n.top=200, plots.path="PlotsRUNedgeR_22Aug", plot.name="", AT.id, UP.genes, DOWN.genes,  varialbs.plot = i, fit.clr = 4, LRTcoef=2)


run.edgeR(x, new.samps, model.formula = as.formula(paste("~", i , "+ drough.control")) , varialbs=c(i), elim.samps=c(flowered.samps), n.top=200, plots.path="PlotsRUNedgeR_22Aug", plot.name="LRTcoef_2_3", AT.id, UP.genes, DOWN.genes,  varialbs.plot = i, fit.clr =4, LRTcoef=c(2,3))


run.edgeR(x, new.samps, model.formula =as.formula(paste("~", i , ": drough.control")) , varialbs=c(i), elim.samps=c(flowered.samps), n.top=200, plots.path="PlotsRUNedgeR_22Aug", plot.name="LRTcoef_3", AT.id, UP.genes, DOWN.genes, varialbs.plot = i, fit.clr =4, LRTcoef=c(3))

run.edgeR(x, new.samps, model.formula = as.formula(paste("~", i , ": drough.control")) , varialbs=c(i), elim.samps=c(flowered.samps), n.top=200, plots.path="PlotsRUNedgeR_22Aug", plot.name="LRTcoef_2_3", AT.id, UP.genes, DOWN.genes, varialbs.plot = i, fit.clr = 4, LRTcoef=c(2,3))


run.edgeR(x, new.samps, model.formula = as.formula(paste("~", i , "* drough.control"))  , varialbs=c(i), elim.samps=c(flowered.samps), n.top=200, plots.path="PlotsRUNedgeR_22Aug", plot.name="LRTcoef_3_4", AT.id, UP.genes, DOWN.genes,  varialbs.plot = i, fit.clr =4, LRTcoef=c(3,4))

run.edgeR(x, new.samps, model.formula =as.formula(paste("~", i , "* drough.control")) , varialbs=c(i), elim.samps=c(flowered.samps), n.top=200, plots.path="PlotsRUNedgeR_22Aug", plot.name="LRTcoef_2_3_4", AT.id, UP.genes, DOWN.genes,  varialbs.plot = i, fit.clr =4, LRTcoef=c(2,3,4))


}

### full models

run.edgeR(x, new.samps, model.formula = ~ Soil.Moisture + Water.Potential + Temperature, varialbs=c("Soil.Moisture", "Water.Potential", "Temperature"), elim.samps=flowered.samps, n.top=200, plots.path="PlotsRUNedgeR_22Aug/Full_Model", plot.name="", AT.id, UP.genes, DOWN.genes,  varialbs.plot = "Water.Potential", fit.clr = 6, LRTcoef=2)

run.edgeR(x, new.samps, model.formula = ~ Soil.Moisture + Water.Potential + Temperature, varialbs=c("Soil.Moisture", "Water.Potential", "Temperature"), elim.samps=c(control.samps, flowered.samps), n.top=200, plots.path="PlotsRUNedgeR_22Aug/Full_Model", plot.name="NoContrl", AT.id, UP.genes, DOWN.genes,  varialbs.plot = "Water.Potential", fit.clr = 6, LRTcoef=2)


run.edgeR(x, new.samps, model.formula = ~ Soil.Moisture + Water.Potential + Temperature + drough.control, varialbs=c("Soil.Moisture", "Water.Potential", "Temperature"), elim.samps=c(flowered.samps), n.top=200, plots.path="PlotsRUNedgeR_22Aug/Full_Model", plot.name="", AT.id, UP.genes, DOWN.genes,  varialbs.plot = "Water.Potential", fit.clr = 6, LRTcoef=2)





