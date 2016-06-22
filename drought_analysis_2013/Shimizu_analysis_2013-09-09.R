####################################################

# preparation of variables describing samples:
# I
# -Soild Moisture: logit + loess
# -Water Potential: extrapolation for May samps; logit + loess 
# -Temperature: exact values;  max, min, avg
# II
# -Soild Moisture: C: loess + interp; DE: interp
# -Water Potential: avg for samples from May, interpolation
# -Temperature: exact values;  max, min, avg


# Correlation between factors
# plots of expression for flowering genes

## MODEL SELECTION: backwise (p-values) and min AIC criteria

# clustering
# GO - topGO

####################################################


#####################################################################################################

### preparation of variables describing samples 

#####################################################################################################

setwd("/Users/gosia/Analysis/Shimizu_RNA_seq/")

####################################################
### samps
####################################################

### parse sample information
samps <- read.table("data/sample_list.csv", sep=";", stringsAsFactors=FALSE, header=TRUE)
head(samps)
names(samps)

# no genes variables &  no samps from November 2009
new.samps <- samps[!samps$sample_name %in% c("N9_970_20090114","M1_1099_20090114","I5_8212_20090119","F8_8266_20090119"), c(1:16)]

# no genes variables
#new.samps <- samps[, c(1:16)]
write.table(new.samps,"Samples_out/samps.csv", quote=FALSE, sep=";", row.names=FALSE)

### time format change 
new.samps$time_ch <- strptime(new.samps$year.month.day, "%d.%m.%y")
new.samps$time_nr <- as.numeric(new.samps$time_ch)

new.samps$flowered <- new.samps$tree_ID=="8266"

names(new.samps)

### list of all unique days in 2008 and 2009
ad <- read.table("Data/Unique_days_short.csv", sep=";")
all.days <- as.numeric(strptime(ad[,], "%d.%m.%y"))


# columns 8-16 come from Lambir_meteorological_data_...xls


# FUN that matches two tables 
match.func <- function(DF1, col.match1, DF2, col.match2, col.values2){
  
  DF1[, col.values2] <- NA
  
  for(i in unique(DF1[,col.match1])){
    # i=unique(DF1[,col.match1])[1]
    DF1[DF1[,col.match1]==i, col.values2] <- ifelse(length(DF2[DF2[,col.match2]==i, col.values2])==0, NA, DF2[DF2[,col.match2]==i, col.values2][1])
    
  }
  invisible(DF1)
}


### FUN that do logit transformation, fits loess and calculates rolled means over roll.days days
logit.loess.and.roll <- function(intrp.x, intrp.y, intrp.points, table1, roll.days=c(14, 28), col.values="Soil.Moisture", plots.path="Plots_Loess", plot.name="", span=0.2, logit=TRUE, ylim=c(0,1))
  {
  # intrp.x, intrp.points: values of "time_nr" type
  # table1: has to have "time_nr" column
  # RETURNS: table1 with interpolated intrp.y values at table1[,"time_nr"] time
  
  # FUN that matches two tables 
  match.func <- function(DF1, col.match1, DF2, col.match2, col.values2){
    
    DF1[, col.values2] <- NA
    
    for(i in unique(DF1[,col.match1])){
      # i=unique(DF1[,col.match1])[1]
      DF1[DF1[,col.match1]==i, col.values2] <- ifelse(length(DF2[DF2[,col.match2]==i, col.values2])==0, NA, DF2[DF2[,col.match2]==i, col.values2][1])
      
    }
    invisible(DF1)
  }
  
  intrp.points <- sort(intrp.points)
  
  if(logit){
    intrp.y <- ifelse(intrp.y==0, min(intrp.y[intrp.y!=0]), intrp.y)
    intrp.y <- ifelse(intrp.y==1, max(intrp.y[intrp.y!=1]), intrp.y)
    intrp.y <- log(intrp.y/(1-intrp.y))
  }
  
  data <- data.frame(intrp.x, intrp.y)
  data <- data[!is.na(data[,2]),]
  data <- data[order(data[,1]), ]
  colnames(data) <- c("time_nr", col.values)
  
  ### loess
  table2.loess <- loess(as.formula(paste(col.values, "~", "time_nr")), data, span = span, control = loess.control(surface = "direct"))
  table2.pred <- predict(table2.loess , intrp.points, se = TRUE)
  
  if(logit){
    table2.pred$fit <- exp(table2.pred$fit)/(exp(table2.pred$fit)+1)
    intrp.y <- exp(intrp.y)/(exp(intrp.y)+1)
  }
  table2 <- data.frame(intrp.points, table2.pred$fit)
  colnames(table2) <- c("time_nr", col.values)
  
  ### rollmean over 14 and 28 days
  library(zoo)
  
  table1 <- match.func(DF1=table1, col.match1="time_nr", DF2=table2, col.match2="time_nr", col.values2=col.values)
  
  for(r in roll.days){
  
  col.values.r <- paste(col.values, r, sep="")
  table2[,col.values.r] <- NA
  table2[r:nrow(table2), col.values.r] <- rollmean(table2[,col.values], r)
  
  table1 <- match.func(DF1=table1, col.match1="time_nr", DF2=table2, col.match2="time_nr", col.values2=col.values.r)
  
  }
  
  dir.create(plots.path, recursive=T, showWarnings=FALSE)
  
  library(stringr)
  
  pdf( paste(plots.path, "/Loess_", str_replace_all(paste(col.values, plot.name, sep=""),"[[:punct:]]", "_"),".pdf", sep=""), width = 10, height = 5)
  plot(intrp.x, intrp.y, pch=20, col=1, main=plot.name ,ylab=col.values, xlab="Time", ylim=ylim, xlim=c(1217541600,1243720800))
  abline(v=table1[,"time_nr"], col="grey")
  lines(intrp.points, table2.pred$fit, col=3)
  points(table1[,"time_nr"], table1[,col.values], col=4)
  dev.off()
  
  invisible(table1)
  
}


# intrp.x=sm$time_nr
# intrp.y=sm$DE8266
# intrp.points=all.days
# table1=new.samps[new.samps$tree_ID == 8266,]
# col.values="Soil.Moisture"
# plots.path="Plots_Interpolate"
# plot.name="DE8266"
# roll.days=c(14, 28)
# ylim=c(0,1)

### FUN that interpolates and calculates rolled means 
interpolate.and.roll <- function(intrp.x, intrp.y, intrp.points, table1, roll.days=c(14, 28), col.values="Soil.Moisture", plots.path="Plots_Interpolate", plot.name="", ylim=c(0,1)){
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
  
  table1 <- match.func(DF1=table1, col.match1="time_nr", DF2=table2, col.match2="time_nr", col.values2=col.values)
  
  ### rollmean over 14 and 28 days
  for(r in roll.days){
  
  col.values.r <- paste(col.values, r , sep="")
  table2[,col.values.r] <- NA
  table2[r:nrow(table2), col.values.r] <- rollmean(table2[,col.values], r)
  
  table1 <- match.func(DF1=table1, col.match1="time_nr", DF2=table2, col.match2="time_nr", col.values2=col.values.r)

  }
  
  dir.create(plots.path, recursive=T, showWarnings=FALSE)
  
  library(stringr)
  
#   pdf( paste(plots.path, "/Interp2_", str_replace_all(paste(col.values, plot.name, sep=""),"[[:punct:]]", "_"),".pdf", sep=""), width = 10, height = 5)
#   plot(intrp.x, intrp.y, col=1, pch=20, ylab=col.values, xlab="Time")
#   abline(v=table1[,"time_nr"], col="grey")
#   lines(interp.intrp.points, col=6)
#   points(table1[,"time_nr"], table1[,col.values], col=4)
#   dev.off()
  
  pdf( paste(plots.path, "/InterpRoll_", str_replace_all(paste(col.values, plot.name, sep=""),"[[:punct:]]", "_"),".pdf", sep=""), width = 10, height = 5)
  plot(intrp.x, intrp.y, col=1, pch=20, ylab=col.values, xlab="Time", ylim=ylim, xlim=c(1217541600,1243720800))
  abline(v=table1[,"time_nr"], col="grey")
  lines(table2[,"time_nr"], table2[,col.values], col=6)
  for(r in roll.days){
    col.values.r <- paste(col.values, r , sep="")
    lines(table2[,"time_nr"], table2[,col.values.r], col=6,  lty=(which(roll.days==r)+1))   
  }
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
# dir.create(path="Plots_Raw", showWarnings=FALSE, recursive=TRUE)
# 
# pdf("Plots_Raw/Soil_moisture.pdf", width = 10, height = 5)
# plot(sm$time_nr, sm$DE970, type="l", col=1, ylim=c(-0.2,1.1), main="Soil Moisture", ylab="Soil moisture", xlab="Time")
# lines(sm$time_nr, sm$C970, col=1, lty=2)
# lines(sm$time_nr, sm$DE8266, col=3)
# lines(sm$time_nr, sm$C8266, col=3, lty=2)
# lines(sm$time_nr, sm$DE8212, col=4)
# lines(sm$time_nr, sm$C8212a, col=4, lty=2)
# lines(sm$time_nr, sm$C8212b, col=4, lty=3)
# legend("bottomleft", c("DE970", "C970","DE8266","C8266" ,"DE8212","C8212a" ,"C8212b"), col=c(1, 1, 3, 3, 4, 4, 4), lty=c(1,2,1,2,1,2,3), cex=0.5)
# dev.off()
# 
# 
# pdf("Plots_Raw/Soil_moisture_Control.pdf", width = 10, height = 5)
# plot(sm$time_nr, sm$C970, type="l", col=2, ylim=c(0,1.1), xlim=c(1.220e+09, 1257289200), main="Control trees", ylab="Soil moisture", xlab="Time")
# lines(sm$time_nr, sm$C8266, col=3)
# lines(sm$time_nr, sm$C8212a, col=4)
# lines(sm$time_nr, sm$C8212b, col=5)
# abline(v=new.samps$time_nr[new.samps$tree_ID %in% c(990, 1099, 1377)], col="grey")
# dev.off()
# 
# pdf("Plots_Raw/Soil_moisture_DE970.pdf", width = 10, height = 5)
# plot(sm$time_nr, sm$DE970, type="l", col=2, ylim=c(0,1.1), xlim=c(1.220e+09, 1257289200), main="DE970", ylab="Soil moisture", xlab="Time")
# abline(v=new.samps$time_nr[new.samps$tree_ID==970], col="grey")
# dev.off()
# 
# pdf("Plots_Raw/Soil_moisture_DE8266.pdf", width = 10, height = 5)
# plot(sm$time_nr, sm$DE8266, type="l", col=2, ylim=c(0,1.1), xlim=c(1.220e+09, 1257289200), main="DE8266", ylab="Soil moisture", xlab="Time")
# abline(v=new.samps$time_nr[new.samps$tree_ID==8266], col="grey")
# dev.off()
# 
# pdf("Plots_Raw/Soil_moisture_DE8212.pdf", width = 10, height = 5)
# plot(sm$time_nr, sm$DE8212, type="l", col=2, ylim=c(0,1.1), xlim=c(1.220e+09, 1257289200), main="DE8212", ylab="Soil moisture", xlab="Time")
# abline(v=new.samps$time_nr[new.samps$tree_ID==8212], col="grey")
# dev.off()


# merge control samples
control <- rbind(as.matrix(sm[!is.na(sm$C970), c("time_nr", "C970")]),
                 as.matrix(sm[!is.na(sm$C8266), c("time_nr", "C8266")]),
                 as.matrix(sm[!is.na(sm$C8212a), c("time_nr", "C8212a")]),
                 as.matrix(sm[!is.na(sm$C8212b), c("time_nr", "C8212b")]))
control <- as.data.frame(control)
colnames(control) <- c("time_nr", "Soil.Moisture")
control <- control[order(control$time_nr), ]


### logit.loess.and.roll

# new.samps.control <- logit.loess.and.roll(intrp.x=control[,1], intrp.y=control[,2], intrp.points=all.days, table1=new.samps[new.samps$drough.control=="control",], roll.days=c(14, 28),col.values="Soil.Moisture", plots.path="Plots_Loess", plot.name="Control", span=3/nrow(control)*6, logit=TRUE, ylim=c(0,1))
# 
# new.samps.DE970 <- logit.loess.and.roll(intrp.x=sm$time_nr, intrp.y=sm$DE970, intrp.points=all.days, table1=new.samps[new.samps$tree_ID == 970,], roll.days=c(14, 28),col.values="Soil.Moisture", plots.path="Plots_Loess", plot.name="DE970", span=3/nrow(sm)*8, logit=TRUE, ylim=c(0,1))
# 
# new.samps.DE8212 <- logit.loess.and.roll(intrp.x=sm$time_nr, intrp.y=sm$DE8212, intrp.points=all.days, table1=new.samps[new.samps$tree_ID == 8212,], roll.days=c(14, 28),col.values="Soil.Moisture", plots.path="Plots_Loess", plot.name="DE8212", span=3/nrow(sm)*8, logit=TRUE, ylim=c(0,1))
# 
# new.samps.DE8266 <- logit.loess.and.roll(intrp.x=sm$time_nr, intrp.y=sm$DE8266 ,intrp.points=all.days, table1=new.samps[new.samps$tree_ID == 8266,], roll.days=c(14, 28),col.values="Soil.Moisture", plots.path="Plots_Loess", plot.name="DE8266", span=3/nrow(sm)*8, logit=TRUE, ylim=c(0,1))
# 
# new.samps <- rbind(new.samps.control, new.samps.DE970, new.samps.DE8212, new.samps.DE8266)


### lowess for control
control.l <- lowess(control[,1], control[,2], f=0.01)

pdf("Plots_Interpolate/Soil_moisture_Control_lowess.pdf", width = 10, height = 5)
plot(control, col=2, ylim=c(0,1.1), pch=20, main="Control trees", ylab="Soil moisture", xlab="Time")
lines(control.l)
dev.off()



new.samps.control <- interpolate.and.roll(intrp.x=control.l$x, intrp.y=control.l$y, intrp.points=all.days, table1=new.samps[new.samps$drough.control=="control",], roll.days=c(14, 28), col.values="Soil.Moisture", plot.name="Control")

new.samps.DE970 <- interpolate.and.roll(intrp.x=sm$time_nr, intrp.y=sm$DE970, intrp.points=all.days, table1=new.samps[new.samps$tree_ID == 970,], roll.days=c(14, 28), col.values="Soil.Moisture", plot.name="DE970")

new.samps.DE8212 <- interpolate.and.roll(intrp.x=sm$time_nr, intrp.y=sm$DE8212, intrp.points=all.days, table1=new.samps[new.samps$tree_ID == 8212,], roll.days=c(14, 28), col.values="Soil.Moisture", plot.name="DE8212")

new.samps.DE8266 <- interpolate.and.roll(intrp.x=sm$time_nr, intrp.y=sm$DE8266, intrp.points=all.days, table1=new.samps[new.samps$tree_ID == 8266,], roll.days=c(14, 28), col.values="Soil.Moisture", plot.name="DE8266")

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


#tree.ids <- unique(as.character(new.samps[,"tree_ID"]))
tree.ids <- c("970" , "1099", "8266", "8212", "1377", "990" )
exp.tree.ids <-  c("DE970" , "C1099", "DE8266", "DE8212", "C1377", "C990" )

### PLOTS of raw WP & sampling points 
for(i in 1:length(tree.ids)){
  # i=5 
  tree.sampl <- paste(exp.tree.ids[i], "sMean", sep="")
  
  pdf(paste("Plots_Raw/Water_potential_", exp.tree.ids[i], ".pdf", sep=""), width = 10, height = 5)
  plot(wp$time_nr[!is.na(wp[,tree.sampl])], wp[!is.na(wp[,tree.sampl]), tree.sampl], type="l", col=3, ylim=c(0,1), xlim=c(1.220e+09, 1257289200), main=exp.tree.ids[i], ylab="Mean Water Potential", xlab="Time")
  abline(v=new.samps$time_nr[new.samps$tree_ID==tree.ids[i]], col="grey")
  dev.off()
  
  
}


pdf(paste("Plots_Raw/Water_potential.pdf", sep=""), width = 10, height = 5)
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

#tree.ids <- unique(as.character(new.samps[,"tree_ID"]))
tree.ids <- c("970" , "1099", "8266", "8212", "1377", "990" )
exp.tree.ids <-  c("DE970" , "C1099", "DE8266", "DE8212", "C1377", "C990" )

### loess
tree.ids <- c("970" , "1099", "8266", "8212", "1377", "990" )
span <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.6)


new.samps.wp <- NULL

for(i in 1:length(tree.ids)){
  # i=6
  tree.sampl <- paste(exp.tree.ids[i], "sMean", sep="")
  
  new.samps.wp <- rbind(new.samps.wp, new.samps.wp.tmp <- logit.loess.and.roll(intrp.x=wp$time_nr[!is.na(wp[,tree.sampl])], intrp.y=wp[!is.na(wp[,tree.sampl]),tree.sampl], intrp.points=all.days, table1=new.samps[new.samps$tree_ID == tree.ids[i],], roll.days=c(14, 28), col.values="Water.Potential", plots.path="Plots_Loess", plot.name=paste(tree.sampl, "sp07", sep=""), span=0.7, logit=TRUE, ylim=c(0,1)))
  
}

### cubic splines smoothing
# for(i in 1:length(tree.ids)){
#   # i=3
#   tree.sampl <- paste(exp.tree.ids[i], "sMean", sep="")
#   
#   intrp.x=wp$time_nr[!is.na(wp[,tree.sampl])]
#   intrp.y=wp[!is.na(wp[,tree.sampl]),tree.sampl]
#   intrp.points=all.days
#   
#   ss <- smooth.spline(intrp.x, intrp.y, spar=0.2)
#   
#   pred <- predict(ss, intrp.points)
#   
#   plot(intrp.x, intrp.y, xlim=c(min(all.days), max(all.days)), ylim=c(0,1))
#   lines(pred)
#   
# }

new.samps <- new.samps.wp



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

tree.ids <- c("970" , "1099", "8266", "8212", "1377", "990" )
exp.tree.ids <-  c("DE970" , "C1099", "DE8266", "DE8212", "C1377", "C990" )

new.samps.wp <- NULL

### intepolation
for(i in 1:length(tree.ids)){
  # i=1
  tree.sampl <- paste(exp.tree.ids[i], "sMean", sep="")
  new.samps.wp <- rbind(new.samps.wp, new.samps.wp.tmp <- interpolate.and.roll(intrp.x=wp$time_nr[!is.na(wp[,tree.sampl])], intrp.y=wp[!is.na(wp[,tree.sampl]),tree.sampl], intrp.points=all.days, table1=new.samps[new.samps$tree_ID == tree.ids[i],] , roll.days=c(14, 28), col.values="Water.Potential", plot.name=tree.sampl))
  
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
abline(v=new.samps$time_nr, col="grey")
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

### some plots of avg, max, min temperature

#plot(tempr.day$TempAvg, tempr.day$TempMax)

pdf("Plots_Raw/Temperature_Avg.pdf", width = 10, height = 5)
plot(tempr.day$time_nr, tempr.day$TempAvg, pch=20, main="Average Temperature", xlab="Time", ylab="Temperature")
abline(v=new.samps$time_nr, col="grey")
dev.off()

pdf("Plots_Raw/Temperature_Max.pdf", width = 10, height = 5)
plot(tempr.day$time_nr, tempr.day$TempMax, pch=20, main="Max Temperature between 10 and 18", xlab="Time", ylab="Temperature")
abline(v=new.samps$time_nr, col="grey")
dev.off()

pdf("Plots_Raw/Temperature_Min.pdf", width = 10, height = 5)
plot(tempr.day$time_nr, tempr.day$TempMin, pch=20, main="Min Temperature between 10 and 18", xlab="Time", ylab="Temperature")
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

write.table(new.samps, "Samples_out/new_samps_interpolation.csv", sep=";",  row.names = F)

####################################################
### Correlation between factors
####################################################

setwd("/Users/gosia/Analysis/Shimizu_RNA_seq/")
load(file="Samples_out/new_samps_interpolation.RData")

names(new.samps)

library(graphics)

new.samps$tree_ID_nr <- as.factor(new.samps$tree_ID)
levels(new.samps$tree_ID_nr) <- 1:length(levels(new.samps$tree_ID_nr))

pdf("Plots_Interpolate/Correlation.pdf", width = 10, height = 10)
pairs(new.samps[,c("Soil.Moisture", "Water.Potential", "TempAvg")], col=new.samps$tree_ID_nr, pch=19)
dev.off()

pdf("Plots_Interpolate/Correlation28.pdf", width = 10, height = 10)
pairs(new.samps[,c("Soil.Moisture28", "Water.Potential28", "TempAvg28")], col=new.samps$tree_ID_nr, pch=19)
dev.off()

pdf("Plots_Interpolate/Correlation14.pdf", width = 10, height = 10)
pairs(new.samps[,c("Soil.Moisture14", "Water.Potential14", "TempAvg14")], col=new.samps$tree_ID_nr, pch=19)
dev.off()



### Partial correlation coefficient
### a,b=variables to be correlated, c=variable to be partialled out of both
# pcor = function(a,b,c)
# {
#   (cor(a,b)-cor(a,c)*cor(b,c))/sqrt((1-cor(a,c)^2)*(1-cor(b,c)^2))
# }


#####################################################################################################

### Load data 

#####################################################################################################

# setwd("/home/gosia/Shimizu_RNA_seq/")

### Robust packages
# install.packages(pkgs="R/R_packages_Robust_versions/edgeR_4.0.33.tar.gz", lib="R/Library/")
# install.packages(pkgs="R/R_packages_Robust_versions/limma_4.17.12.tar.gz", lib="R/Library/")

# library(package=limma, lib.loc="R/Library/")
# library(package=edgeR, lib.loc="R/Library/")

setwd("/Users/gosia/Analysis/Shimizu_RNA_seq/")

### samples & factors
#load(file="Samples_out/new_samps_interpolation.RData")
new.samps <- read.table("Samples_out/new_samps_interpolation_November.csv", header=TRUE, sep=";", stringsAsFactors=FALSE)

rownames(new.samps) <- new.samps$sample_name
new.samps$tree_ID <- as.factor(new.samps$tree_ID)
### counts
x <- read.table("Data/raw_data-clean.csv", sep=",", header=T, row.names=1)
x <- x[, new.samps$sample_name]


flowered.samps <- c("E7_8266_20090416") # flowered sample
control.samps <- rownames(new.samps[new.samps$drough.control=="control",])
November.samps <- c("N9_970_20090114","M1_1099_20090114","I5_8212_20090119","F8_8266_20090119")

###########################
### files with control genes
###########################
AT.id <- read.table("Data/genes_descr_control/best_hit_blast_result_Sl_predicted_exons.txt", sep=",", stringsAsFactors=FALSE)

# UP.genes <- read.table("Data/genes_descr_control/mDr_Day10_drought_up_regulated_genes_Harb_etal.csv", sep=";", header=T, stringsAsFactors=FALSE)
# DOWN.genes <- read.table("Data/genes_descr_control/mDr_Day10_drought_down_regulated_genes_Harb_etal.csv", sep=";", header=T, stringsAsFactors=FALSE)
# Drought.genes <- read.table("Data/genes_descr_control/mDr_Day10_drought_regulated_genes_Harb_etal.csv", sep=";", header=T, stringsAsFactors=FALSE)

### DE unigenes for S.beccariana overlap with BLAST A.thaliana, all clusters
# MolEcol.DE.genes <- read.table("Data/genes_descr_control/gene_list_in_all_clusters_DEgenes_Mol_Ecol.csv", sep=";", header=T, stringsAsFactors=FALSE)


### Flowering genes - Table S4
Athaliana.flowering.genes <- read.table("Data/genes_descr_control/gene_list_flowering_related_MolEcol_Athaliana_S4.csv", sep=";", header=T, stringsAsFactors=FALSE, skip=1)

Athaliana.flowering.genes[,1] <- gsub(" ", "", Athaliana.flowering.genes[,1])

# "AT1G65480" - FT gene, "AT2G22540" - SVP gene
Athaliana.flowering.genes.FT.SVP <- Athaliana.flowering.genes[Athaliana.flowering.genes[,1] %in% c("AT1G65480", "AT2G22540"), ]



#####################################################################################################

### plots of expression for flowering genes

#####################################################################################################

#AT.genes <- Athaliana.flowering.genes.FT.SVP
AT.genes <- Athaliana.flowering.genes

#elim.samps=c(flowered.samps)
elim.samps=NULL

new.samps <- new.samps[order(new.samps$time_nr), ]

x <- x[, new.samps$sample_name] 
x <- x[,!names(x) %in% elim.samps]
new.samps <- new.samps[!new.samps$sample_name %in% elim.samps, ]

new.samps$tree_ID <- as.factor(new.samps$tree_ID)

library(edgeR)
d.org <- DGEList(x, group=new.samps$tree_ID)
d.org <- calcNormFactors(d.org)


d.cpm <- cpm(d.org, normalized.lib.sizes=TRUE)
d.cpm.l <- log(d.cpm + min(d.cpm[d.cpm != 0]))

####################################
# files per AT gene
####################################

for(g in 1:nrow(AT.genes)){
  # g=1
  genes <- AT.id[AT.id[,2] == AT.genes[g, 1], 1]
  
  pdf(paste0("Plots_of_flowering_genes/" , AT.genes[g, 1] ,".pdf"), h=5, w=10)
  #par(mfrow=c(6,1))

  for(j in 1:length(genes)){
    # j=1
    
    plot(0, type="n", main=paste0(AT.genes[g, 1] ," - ", AT.genes[g, 2] , "\n",  genes[j]) ,xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(d.cpm.l[genes[j], ])), max(na.omit(d.cpm.l[genes[j], ]))), xlab="Time", ylab="Gene Expression in cpm")
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[246])
    for(t in levels(new.samps$tree_ID)){
      lines(new.samps$time_nr[new.samps$tree_ID==t], d.cpm.l[genes[j], new.samps$tree_ID==t] , col=which(levels(new.samps$tree_ID)==t), type="b", pch=ifelse(new.samps$drough.control[new.samps$tree_ID==t]=="drought", 18, 16)) 
    }
    legend("topleft", legend = levels(new.samps$tree_ID), col=1:length(levels(new.samps$tree_ID)), cex=1, text.col= 1:length(levels(new.samps$tree_ID)))
     
  }
  
  
  plot(0, type="n", xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(new.samps$Soil.Moisture)), max(na.omit(new.samps$Soil.Moisture))), xlab="Time", ylab="Soil.Moisture")
  for(t in levels(new.samps$tree_ID)){
    lines(new.samps$time_nr[new.samps$tree_ID==t], new.samps$Soil.Moisture[new.samps$tree_ID==t], col=which(levels(new.samps$tree_ID)==t), type="b", pch=ifelse(new.samps$drough.control[new.samps$tree_ID==t]=="drought", 18, 16)) 
  }
  legend("topleft", legend = levels(new.samps$tree_ID), col=1:length(levels(new.samps$tree_ID)), cex=0.5, text.col= 1:length(levels(new.samps$tree_ID)))
  
  plot(0, type="n", xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(new.samps$Water.Potential)), max(na.omit(new.samps$Water.Potential))), xlab="Time", ylab="Water.Potential")
  for(t in levels(new.samps$tree_ID)){
    lines(new.samps$time_nr[new.samps$tree_ID==t], new.samps$Water.Potential[new.samps$tree_ID==t], col=which(levels(new.samps$tree_ID)==t), type="b", pch=ifelse(new.samps$drough.control[new.samps$tree_ID==t]=="drought", 18, 16)) 
  }
  legend("topleft", legend = levels(new.samps$tree_ID), col=1:length(levels(new.samps$tree_ID)), cex=0.5, text.col= 1:length(levels(new.samps$tree_ID)))
  
  plot(0, type="n", xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(new.samps$TempAvg28)), max(na.omit(new.samps$TempAvg28))), xlab="Time", ylab="TempAvg28")
  for(t in levels(new.samps$tree_ID)){
    lines(new.samps$time_nr[new.samps$tree_ID==t], new.samps$TempAvg28[new.samps$tree_ID==t], col=which(levels(new.samps$tree_ID)==t), type="b", pch=ifelse(new.samps$drough.control[new.samps$tree_ID==t]=="drought", 18, 16)) 
  }
  legend("topleft", legend = levels(new.samps$tree_ID), col=1:length(levels(new.samps$tree_ID)), cex=0.5, text.col= 1:length(levels(new.samps$tree_ID)))
  
  plot(0, type="n", xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(new.samps$TempAvg)), max(na.omit(new.samps$TempAvg))), xlab="Time", ylab="TempAvg")
  for(t in levels(new.samps$tree_ID)){
    lines(new.samps$time_nr[new.samps$tree_ID==t], new.samps$TempAvg[new.samps$tree_ID==t], col=which(levels(new.samps$tree_ID)==t), type="b", pch=ifelse(new.samps$drough.control[new.samps$tree_ID==t]=="drought", 18, 16)) 
  }
  legend("topleft", legend = levels(new.samps$tree_ID), col=1:length(levels(new.samps$tree_ID)), cex=0.5, text.col= 1:length(levels(new.samps$tree_ID)))
  
  #     plot(0, type="n", xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(new.samps$TempMax28)), max(na.omit(new.samps$TempMax28))), xlab="Time", ylab="TempMax28")
  #     for(t in levels(new.samps$tree_ID)){
  #       lines(new.samps$time_nr[new.samps$tree_ID==t], new.samps$TempMax28[new.samps$tree_ID==t], col=which(levels(new.samps$tree_ID)==t), type="b", pch=ifelse(new.samps$drough.control[new.samps$tree_ID==t]=="drought", 18, 16)) 
  #     }
  #     legend("topleft", legend = levels(new.samps$tree_ID), col=1:length(levels(new.samps$tree_ID)), cex=0.5, text.col= 1:length(levels(new.samps$tree_ID)))
  
  dev.off()
  
  
  
}


####################################
### all in one file
####################################

pdf(paste0("Plots_of_flowering_genes/" , "Flowering_genes_from_S4_table_Nov_fl" ,".pdf"), h=5, w=10)

info.table <- NULL

for(g in 1:nrow(AT.genes)){
  # g=8
  cat(paste(g, ", "))
  genes <- AT.id[AT.id[,2] == AT.genes[g, 1], 1]
  
  if(length(genes!=0)){
  for(j in 1:length(genes)){
    # j=1
    
    info.table <- rbind(info.table, data.frame(AT.ID=AT.genes[g, 1], AT.description=AT.genes[g, 2] , ID=genes[j]))
    
    plot(0, type="n", main=paste0(AT.genes[g, 1] ," - ", AT.genes[g, 2] , "\n",  genes[j]) ,xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(d.cpm[genes[j], ])), max(na.omit(d.cpm[genes[j], ]))), xlab="Time", ylab="Gene Expression in cpm")
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[246])
    for(t in levels(new.samps$tree_ID)){
      lines(new.samps$time_nr[new.samps$tree_ID==t], d.cpm[genes[j], new.samps$tree_ID==t] , col=which(levels(new.samps$tree_ID)==t), type="b", pch=ifelse(new.samps$drough.control[new.samps$tree_ID==t]=="drought", 18, 16)) 
    }
    legend("topleft", legend = levels(new.samps$tree_ID), col=1:length(levels(new.samps$tree_ID)), cex=1, text.col= 1:length(levels(new.samps$tree_ID)))
    
  }
  
  }
}

plot.vars <- names(new.samps)[20:34]

for(v in plot.vars){
  
  plot(0, type="n", xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(new.samps[,v])), max(na.omit(new.samps[,v]))), xlab="Time", ylab=v)
  for(t in levels(new.samps$tree_ID)){
    lines(new.samps$time_nr[new.samps$tree_ID==t], new.samps[new.samps$tree_ID==t, v], col=which(levels(new.samps$tree_ID)==t), type="b", pch=ifelse(new.samps$drough.control[new.samps$tree_ID==t]=="drought", 18, 16)) 
  }
  legend("topleft", legend = levels(new.samps$tree_ID), col=1:length(levels(new.samps$tree_ID)), cex=0.5, text.col= 1:length(levels(new.samps$tree_ID)))
  
}

write.table(info.table,"Plots_of_flowering_genes/info_table.xls", quote=FALSE, sep="\t", row.names=FALSE)

dev.off()


### log cpm

pdf(paste0("Plots_of_flowering_genes/" , "Flowering_genes_from_S4_table_Nov_fl_log" ,".pdf"), h=5, w=10)

info.table <- NULL

for(g in 1:nrow(AT.genes)){
  # g=8
  cat(paste(g, ", "))
  genes <- AT.id[AT.id[,2] == AT.genes[g, 1], 1]
  
  if(length(genes!=0)){
    for(j in 1:length(genes)){
      # j=1
      
      info.table <- rbind(info.table, data.frame(AT.ID=AT.genes[g, 1], AT.description=AT.genes[g, 2] , ID=genes[j]))
      
      plot(0, type="n", main=paste0(AT.genes[g, 1] ," - ", AT.genes[g, 2] , "\n",  genes[j]) ,xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(d.cpm.l[genes[j], ])), max(na.omit(d.cpm.l[genes[j], ]))), xlab="Time", ylab="Gene Expression in log(cpm)")
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[246])
      for(t in levels(new.samps$tree_ID)){
        lines(new.samps$time_nr[new.samps$tree_ID==t], d.cpm.l[genes[j], new.samps$tree_ID==t] , col=which(levels(new.samps$tree_ID)==t), type="b", pch=ifelse(new.samps$drough.control[new.samps$tree_ID==t]=="drought", 18, 16)) 
      }
      legend("topleft", legend = levels(new.samps$tree_ID), col=1:length(levels(new.samps$tree_ID)), cex=1, text.col= 1:length(levels(new.samps$tree_ID)))
      
    }
    
  }
}

plot.vars <- names(new.samps)[20:34]

for(v in plot.vars){
  
  plot(0, type="n", xlim=c(min(new.samps$time_nr), max(new.samps$time_nr)), ylim=c(min(na.omit(new.samps[,v])), max(na.omit(new.samps[,v]))), xlab="Time", ylab=v)
  for(t in levels(new.samps$tree_ID)){
    lines(new.samps$time_nr[new.samps$tree_ID==t], new.samps[new.samps$tree_ID==t, v], col=which(levels(new.samps$tree_ID)==t), type="b", pch=ifelse(new.samps$drough.control[new.samps$tree_ID==t]=="drought", 18, 16)) 
  }
  legend("topleft", legend = levels(new.samps$tree_ID), col=1:length(levels(new.samps$tree_ID)), cex=0.5, text.col= 1:length(levels(new.samps$tree_ID)))
  
}

write.table(info.table,"Plots_of_flowering_genes/info_table.xls", quote=FALSE, sep="\t", row.names=FALSE)

dev.off()



#####################################################################################################

### fitting models

#####################################################################################################


### FUN run.edgeR.models.fit()

run.edgeR.models.fit <- function(x, new.samps, varialbs, elim.samps, out.path="Models", out.name="Models_fitting"){ 
  # fits all possible models build from variabls
  
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
  
  models.results <- list()
  
  for(c in length(varialbs):1){
    # c=2
    models.comb <- combn(varialbs, c)
    
    for(m in 1:ncol(models.comb)){
      # m=1
      
      model.formula <- as.formula(paste0("~", paste(models.comb[,m], collapse="+")))
      cat("Fitting: ", paste(models.comb[,m], collapse="+"), "\n")
      design <- model.matrix(model.formula, data=new.samps)
      
      # estimate dispersion  
      d <- estimateGLMCommonDisp(d,design)
      d <- estimateGLMTrendedDisp(d,design)
      d <- estimateGLMTagwiseDisp(d,design)
      
      # glmFit
      fit <- glmFit(d,design)
      
      # make contrast to test any difference from control
      #mc <- makeContrasts(contrasts=c("Intercept", "Soil.Moisture"),levels=c("Intercept", "Soil.Moisture"))
      
      mc <- diag(rep(1, length(models.comb[,m])+1))    
      rownames(mc) <- colnames(mc) <- c("Intercept", models.comb[,m])
      
      LRT <- list()
      for(v in 1:length(models.comb[,m])){
        # v=1
        cat(" Testing: ", colnames(mc)[1+v], "\n")
        lrt <- glmLRT(fit, contrast=mc[,1+v]) 
        LRT[[colnames(mc)[1+v]]] <- lrt[[c("table")]]  
      }
      
      models.results[[paste(models.comb[,m], collapse="+")]] <- list(COEFFS = fit$coefficients, DEV = fit$deviance, DES = fit$design, LRT = LRT)
      
      #str(models.results)
      
    }   
  }
  
  dir.create(out.path, recursive=T, showWarnings = FALSE)
  save(models.results, file=paste(out.path, "/", out.name,".RData", sep=""))
  invisible(models.results)
  
}


Mf <- run.edgeR.models.fit(x, new.samps, varialbs=c("Soil.Moisture", "Water.Potential", "TempAvg"), elim.samps=c(control.samps, flowered.samps), out.path="Models", out.name="Models_fitting_NoControl")


### FUN run.edgeR.model.select

run.edgeR.model.select <- function(Mf, tr=0.01, AT.id, out.path="Models", out.name="Models_selecting"){
  # Mf - output from run.edgeR.models.fit
  # perform backwise selection (nested models)
  
  full.model <- names(Mf)[1] 
  genes <- row.names(Mf[[1]]$LRT[[1]])
  gene.models <- list()
  
  for(g in 1:length(genes)){
    # g <- which(genes=="GID062074_3539215")   
    models.tree <- vector()
    dev <- vector()
    Chi.S <- vector()
    Chi.pv <- vector()
    
    m <- 1
    models.tree[1] <- full.model
    model.var <- unlist(strsplit(models.tree[m], split="+", fixed=TRUE))
    lrt.pvs <- matrix(NA, length(model.var), length(model.var))
    colnames(lrt.pvs) <- model.var
    
    dev[m] <- Mf[[models.tree[m]]]$DEV[g]  
    lrt.pvs.temp <- sapply(model.var, function(v) Mf[[models.tree[1]]]$LRT[[v]][g,"PValue"])
    lrt.pvs[m, names(lrt.pvs.temp)] <- lrt.pvs.temp
    
    while(max(na.omit(lrt.pvs[m,])) >= tr && m < ncol(lrt.pvs)){
      
      var.elim <- model.var[ which.max(na.omit(lrt.pvs[m,])) ]
      
      model.var <- setdiff(model.var, var.elim)
      models.tree[m+1] <- paste(model.var, collapse="+", sep="")
      
      dev[m+1] <- Mf[[models.tree[m+1]]]$DEV[g]
      
      Chi.S[m] <- dev[m+1]-dev[m]
      Chi.pv[m] <- pchisq(q=Chi.S[m], df=1, lower.tail = FALSE)
      
      lrt.pvs.temp <- sapply(model.var, function(v) Mf[[models.tree[m+1]]]$LRT[[v]][g,"PValue"])
      lrt.pvs[m+1, names(lrt.pvs.temp)] <- lrt.pvs.temp
      
      m <- m+1
    }
    
    if(m < ncol(lrt.pvs)){
      gene.models[[genes[g]]] <- list(lrt.pvs=lrt.pvs, Chi.pv=Chi.pv, final.model=models.tree[m])
    }else{
      if(na.omit(lrt.pvs[m,]) <= tr){
        gene.models[[genes[g]]] <- list(lrt.pvs=lrt.pvs, Chi.pv=Chi.pv, final.model=models.tree[m])
      }else{
        gene.models[[genes[g]]] <- list(lrt.pvs=lrt.pvs, Chi.pv=Chi.pv, final.model="Intercept")
      }
      
    }
    
  } # end for
  
  ### genes clustered by model 
  genes <- names(gene.models)
  models <- c(names(models.results), "Intercept")
  
  clusters <- vector("list", length(models))
  names(clusters) <- models
  
  clusters.AT <- vector("list", length(models))
  names(clusters.AT) <- models
  
  for(g in genes){
    # g="GID000003_777"
    clusters[[gene.models[[g]]$final.model]] <- c(clusters[[gene.models[[g]]$final.model]], g)
    
    if(g %in% AT.id[,1])
      clusters.AT[[gene.models[[g]]$final.model]] <- c(clusters.AT[[gene.models[[g]]$final.model]], AT.id[g == AT.id[,1],2])
    
  }
  
  save(gene.models, clusters, clusters.AT, file=paste(out.path, "/", out.name,".RData", sep=""))
  
  invisible(list(gene.models=gene.models, clusters=clusters, clusters.AT=clusters.AT))
  
}


Ms <- run.edgeR.model.select(Mf, tr=0.01, AT.id, out.path="Models", out.name="Models_selecting_NoControl")


### FUN run.edgeR.models.fit.AIC()

run.edgeR.models.fit.AIC <- function(x, new.samps, var.roll, elim.samps, out.path="ModelsAIC", out.name="Models_fit_AIC", nr.cores=2){ 
  # fits all possible models build from variabls
  library(parallel)
  x <- x[, new.samps$sample_name] 
  x <- x[,!names(x) %in% elim.samps]
  new.samps <- new.samps[!new.samps$sample_name %in% elim.samps, ]
  
  library(edgeR)
  d.org <- DGEList(x, group=new.samps$tree_ID)
  d.org <- calcNormFactors(d.org)
  
  ### make sure a gene is expressed (CPM > 1) in more than 2 samples
  d.cpm <- cpm(d.org, normalized.lib.sizes=TRUE)
  d.org <- d.org[ rowSums(d.cpm>1) > 2, ]
  
  genes <- rownames(d.org$counts)
  
  models.results <- NULL
  
  for(c in length(var.roll):1){
    # c=2
    models.comb <- combn(names(var.roll), c)
    
    for(m in 1:ncol(models.comb)){
      # m=1      
      models.all <- as.matrix(expand.grid(var.roll[models.comb[,m]]))

      models.results <- c( models.results, mclapply(1:nrow(models.all), function(ma){  
        # ma=1
        model <- paste(models.all[ma,], collapse="+")
        model.formula <- as.formula( paste0("~", model))
        cat("Fitting: ", model, "\n")
        design <- model.matrix(model.formula, data=new.samps)   
        d <- d.org
        d <- d[, rownames(design)]
        
        # estimate dispersion  
        d <- estimateGLMCommonDisp(d,design)
        d <- estimateGLMTrendedDisp(d,design)
        d <- estimateGLMTagwiseDisp(d,design)
        
        # glmFit
        fit <- glmFit(d,design)
        fit$model <- model
        
        # AIC
        fit$AIC <- numeric(length = length(genes))
        fit$BIC <- numeric(length = length(genes))
        
        for( g in 1:length(genes)){
          # g=2        
          fit$AIC[g] <- 2*c - 2* sum(log(dnbinom(x = fit$counts[genes[g],], size = 1/fit$dispersion[g], mu = fit$fitted.values[genes[g],])))  
          fit$BIC[g] <- log(length(fit$counts[genes[g],]))*c - 2* sum(log(dnbinom(x = fit$counts[genes[g],], size = 1/fit$dispersion[g], mu = fit$fitted.values[genes[g],]))) 
        }  
        
        fit$counts <- NULL
        return(fit)
        
#         return(list(model=fit$model ,coefficients=fit$coefficients, fitted.values=fit$fitted.values, deviance=fit$deviance, design=fit$design, offset=fit$offset, dispersion=fit$dispersion, AIC=fit$AIC))
             
      }, mc.cores=nr.cores ) )
            
    }
  }
  
  dir.create(out.path, recursive=T, showWarnings = FALSE)
  save(models.results, d.org, file=paste(out.path, "/", out.name,".RData", sep=""))
  invisible(list(models.results=models.results, d.org=d.org))
  
}


MAIC <- run.edgeR.models.fit.AIC(x, new.samps, var.roll=list(SM=c("Soil.Moisture", "Soil.Moisture14", "Soil.Moisture28"), WP=c("Water.Potential", "Water.Potential14", "Water.Potential28"), Temp=c("TempAvg", "TempAvg14", "TempAvg28", "TempMax", "TempMax14", "TempMax28", "TempMin", "TempMin14", "TempMin28")), elim.samps=c(control.samps, flowered.samps), out.path="ModelsAIC", out.name="Models_fit_AIC_BIC", nr.cores=15)



##### FUN run.edgeR.models.select.AIC

setwd("/Users/gosia/Analysis/Shimizu_RNA_seq/")
# setwd("/home/gosia/Shimizu_RNA_seq/")

load("ModelsAIC/Models_fit_AIC_BIC.RData")


run.edgeR.models.select.AIC <- function(models.results, d.org, out.path="ModelsAIC", out.name="Best_Models_fit_AIC", nr.cores=10){
  
  AIC <- NULL
  BIC <- NULL
  models <- NULL
  
  for(i in 1:length(models.results)){
    models  <- c(models, models.results[[i]]$model)
    AIC <- cbind(AIC, models.results[[i]]$AIC)
    BIC <- cbind(BIC, models.results[[i]]$BIC)
  }
  
  genes <- rownames(d.org$counts)
  rownames(AIC) <- genes
  rownames(BIC) <- genes
  
  model.minAIC <- apply(AIC, 1, which.min)
  model.minBIC <- apply(BIC, 1, which.min)
  
  genes.modelsAIC <- data.frame(genes = genes, models = models[model.minAIC])
  #table(genes.modelsAIC$models)
  
  names(models.results) <- unlist(lapply(models.results, function(f) f$model))
  
  genes.modelsAIC$testV1 <- NA
  genes.modelsAIC$testV2 <- NA
  genes.modelsAIC$testV3 <- NA
  
  library(parallel)
  library(edgeR)

  genes.modelsAIC.mc <- mclapply(1:length(genes.modelsAIC[,"genes"]), function(g){
    cat("Testing gene nr", g, "out of ", length(genes.modelsAIC[,"genes"]), "\n")
    
    fit <- models.results[[as.character(genes.modelsAIC[g,"models"])]]
    fit$counts <- d.org$counts[, colnames(fit$fitted.values)]
    
    mc <- diag(rep(1, ncol(fit$design)))    
    rownames(mc) <- colnames(mc) <- colnames(fit$design)
    
    LRT <- list()
    for(v in 1:(ncol(fit$design)-1)){
      # v=1
      # cat(" Testing: ", colnames(mc)[1+v], "\n")
      lrt <- glmLRT(fit, contrast=mc[,1+v]) 
      LRT[[colnames(mc)[1+v]]] <- lrt$table[genes.modelsAIC[g,"genes"], "PValue"] 
    }
    
    #LRT
    lrt.pvs <- unlist(LRT)
    testV <- c("testV1", "testV2", "testV3")
    
    for(i in 1:length(lrt.pvs)){
      genes.modelsAIC[g, testV[i]] <- lrt.pvs[i]    
    }    
    return(genes.modelsAIC[g,])
    
  }, mc.cores=nr.cores)

  save(genes.modelsAIC.mc, file=paste0(out.path, "/", out.name, ".Rdata"))
  
  genes.modelsAIC.mc <- do.call("rbind", genes.modelsAIC.mc) 
  
  write.table(genes.modelsAIC.mc, paste0(out.path, "/", out.name, ".xls"), sep="\t", row.names = FALSE, quote = FALSE)
  
  invisible(genes.modelsAIC.mc)
  
}



run.edgeR.models.select.AIC(models.results, d.org, out.path="ModelsAIC", out.name="Best_Models_fit_AIC", nr.cores=5)



#####################################################################################################

### clustering

#####################################################################################################

############################################
### PAM only log(cpm) + euclidean dist
# does not distinguish the change of expression
# beecause no normalization per gene
# only parallel plots
############################################

elim.samps=c(flowered.samps, "K4_990_20081205", "K5_990_20090511")
elim.samps=NULL

x <- x[, new.samps$sample_name] 
x <- x[,!names(x) %in% elim.samps]
new.samps <- new.samps[!new.samps$sample_name %in% elim.samps, ]

library(edgeR)
d.org <- DGEList(x, group=new.samps$tree_ID)
d.org <- calcNormFactors(d.org)

### make sure a gene is expressed (CPM > 1) in more than 2 samples
d.cpm <- cpm(d.org, normalized.lib.sizes=TRUE)
d.org <- d.org[ rowSums(d.cpm>1) > 2, ]

d.org$counts <- d.org$counts + 1
genes <- rownames(d.org$counts)

new.samps$tree_ID <- as.factor(new.samps$tree_ID)

design <- model.matrix(~-1+tree_ID, data=new.samps)

# estimate dispersion  
d.org <- estimateGLMCommonDisp(d.org,design)
d.org <- estimateGLMTrendedDisp(d.org,design)
d.org <- estimateGLMTagwiseDisp(d.org,design)

# glmFit
fit <- glmFit(d.org,design)

# make contrast to test any difference from control
mc <- makeContrasts(contrasts=c("(tree_ID970+tree_ID8212)/2-(tree_ID1099+tree_ID1377)/2", 
                                "tree_ID8266-(tree_ID1099+tree_ID1377)/2"),levels=colnames(design))
colnames(mc) <- c("sheet_vs_control","flower_vs_control")

# glmLRT conducts likelihood ratio tests for one or more coefficients in the linear model
lrt <- glmLRT(fit, contrast = mc[,"sheet_vs_control"])
head(lrt$table)
tt <- topTags(lrt, n=nrow(d.org))
ttt <- tt$table

tr=0.1

de.genes <- rownames(tt$table[tt$table[,"FDR"] <= tr,])

### work on cpm and de genes
d.org.cpm <- cpm(d.org[de.genes, ], normalized.lib.sizes=TRUE)

#plot(log(d.org.cpm[,1:2]))

library(cluster)
library(clusterSim)
library(lattice)

pam.x <- log(d.org.cpm[, new.samps[order(new.samps$time_nr), "sample_name"]])
pam.df <- as.data.frame(pam.x)
pam.dist <- dist(x=pam.x , method = "euclidean")

pam.obj <- list()


for(i in 2:15){
  # i=3
  pam.obj[[i]] <- pam(x=pam.dist, k=i) 
  
  pam.cl <- pam.obj[[i]]$clustering
  pam.sill <- silhouette(x=pam.cl, dist=pam.dist)
  
  pam.obj[[i]]$index.S <- index.S(d=pam.dist, cl=pam.cl)
  pam.obj[[i]]$index.G1 <- index.G1(x=pam.x, d=pam.dist, cl=pam.cl, centrotypes="medoids")
  
  pdf(paste0("Clustering/Sillhouette_plot_k" ,i, ".pdf"), w=15, h=10)
  plot(pam.sill, col=1:i)
  dev.off()
    
  # pdf(paste0("Clustering/Parcoord_plot_k" ,i, ".pdf"), w=15, h=10)
  # for(j in 1:i){
  #   parcoord(x=pam.x[pam.cl==j, ], col = 1, lty = 1, var.label = FALSE)
  # }
  # dev.off()
  
  pam.df$clustering <- pam.cl
  
#   png(paste0("Clustering/Parallel1_plot_k" ,i, ".png"), w=1000, h=1000)
#   parallelplot(~pam.df[, -ncol(pam.df)], pam.df, groups=clustering)
#   dev.off()
  
  png(paste0("Clustering/Parallel1_plot_k" ,i, ".png"), w=1000, h=1000)
  plot(parallel(~pam.df[, -ncol(pam.df)], pam.df, groups=clustering))
  dev.off()
  
#   png(paste0("Clustering/Parallel2_plot_k" ,i, ".png"), w=1000, h=1000)
#   parallelplot(~pam.df[, -ncol(pam.df)] | clustering, pam.df)
#   dev.off()
  
  png(paste0("Clustering/Parallel2_plot_k" ,i, ".png"), w=1000, h=1000)
  plot(parallel(~pam.df[, -ncol(pam.df)] | clustering, pam.df))
  dev.off()
  
  
}

save(pam.obj, file="Clustering/PAM_log_cpm/pam_obj_de_genes.Rdata")


############################################
### apply normalization per gene
############################################

#elim.samps=c(flowered.samps, "K4_990_20081205", "K5_990_20090511")
elim.samps=NULL

x <- x[, new.samps$sample_name] 
x <- x[,!names(x) %in% elim.samps]
new.samps <- new.samps[!new.samps$sample_name %in% elim.samps, ]
new.samps$tree_ID <- as.factor(new.samps$tree_ID)

library(limma)
library(edgeR)
d.org <- DGEList(x, group=new.samps$tree_ID)
d.org <- calcNormFactors(d.org)

### make sure a gene is expressed (CPM > 1) in more than 2 samples
d.cpm <- cpm(d.org, normalized.lib.sizes=TRUE)
d.org <- d.org[ rowSums(d.cpm>1) > 2, ]
genes <- rownames(d.org$counts)

# estimate dispersion  
design <- model.matrix(~-1+tree_ID, data=new.samps)
d.org <- estimateGLMCommonDisp(d.org,design)
d.org <- estimateGLMTrendedDisp(d.org,design)
d.org <- estimateGLMTagwiseDisp(d.org,design)

# glmFit
fit <- glmFit(d.org,design)

# make contrast to test any difference from control
mc <- makeContrasts(contrasts=c("(tree_ID970+tree_ID8212)/2-(tree_ID1099+tree_ID1377)/2", 
                                "tree_ID8266-(tree_ID1099+tree_ID1377)/2"),levels=colnames(design))
colnames(mc) <- c("sheet_vs_control","flower_vs_control")

# glmLRT conducts likelihood ratio tests for one or more coefficients in the linear model
lrt <- glmLRT(fit, contrast = mc[,"sheet_vs_control"])
tt <- topTags(lrt, n=nrow(d.org))
tr=0.1
de.genes <- rownames(tt$table[tt$table[,"FDR"] <= tr,])


### work on DE genes
d.sel <- d.org[de.genes, ]

### work on cpm and de genes // cmp for plots in log scale (no zeros)
d.cpm <- cpm(d.sel, normalized.lib.sizes=TRUE)
d.cpm.l <- log(d.cpm + min(d.cpm[d.cpm != 0]))

# plot(d.sel$counts[1,])
# points(d.sel$counts[5,], col=2)
# 
# plot(d.cpm[1,])
# points(d.cpm[5,], col=2)
# 
# plot(d.cpm.l[1, ])
# points(d.cpm.l[5,], col=2)

### normalize to mean=0 and var=1 per gene

normalize.counts <- function(counts, norm.method=c("norm", "01")[1]){
  
  t(apply(counts, 1, function(g){
    
    if(norm.method=="norm"){
      m <- mean(g)
      sd <- sd(g) 
      return((g-m)/sd)
    } else if(norm.method=="01"){
      return((g-min(g))/(max(g)-min(g)))
    }
  }))
}


d.cpm.l.n <- normalize.counts(counts=d.cpm.l, norm.method="01")
d.cpm.n <- normalize.counts(counts=d.cpm, norm.method="01")


# plot(d.cpm.l.n[1,])
# points(d.cpm.l.n[5,], col=2)



### assings weights based on mean-variance trend
### voom ?

design <- model.matrix(~-1+tree_ID, data=new.samps)

d.voom <- voom(counts=d.sel, design = design, plot = TRUE, span=0.7)

d.voom$weights[,1]

plot(d.voom$E[1,])
points(d.voom$E[5,], col=2)

### variance stabilization
### DESeq ?

library(DESeq)
cds <- newCountDataSet(countData=d.org$counts, conditions=new.samps$tree_ID)
cds <- estimateSizeFactors( cds )
cds <- estimateDispersions( cds, method="blind" )
vsd <- getVarianceStabilizedData( cds )
cols <- conditions(cds) == "8266"
plot( rank( rowMeans( vsd[,cols] ) ), genefilter::rowVars( vsd[,cols] ) )



############################################
# PAM
######################

library(cluster)
library(clusterSim)
library(lattice)

library(gplots) 
library(RColorBrewer)
source("~//R_UZH//heatmap2.r")


PAM.clustering <- function(pam.x, new.samps , out.path="Clustering/PAM/", out.name=""){
  
  dir.create(out.path, showWarnings=FALSE, recursive=TRUE)
  
  samps.order <- row.names(new.samps[order(new.samps$drough.control, new.samps$tree_ID, -new.samps$time_nr), ])
  pam.x <- pam.x[, samps.order]

  pam.df <- as.data.frame(pam.x)
  pam.dist <- dist(x=pam.x , method = "euclidean")
#   ylimL <- min(na.omit(pam.x))
#   ylimU <- max(na.omit(pam.x))
  pam.obj <- list()

  for(i in 2:25){
    # i=3
    pam.obj[[i]] <- pam(x=pam.dist, k=i) 
    
    pam.cl <- pam.obj[[i]]$clustering
    pam.sill <- silhouette(x=pam.cl, dist=pam.dist)
    
    pam.obj[[i]]$index.S <- index.S(d=pam.dist, cl=pam.cl)
    pam.obj[[i]]$index.G1 <- index.G1(x=pam.x, d=pam.dist, cl=pam.cl, centrotypes="medoids")
    
    pam.df$clustering <- pam.cl
    
#     pdf(paste0(out.path, "/Sillhouette_plot_k" ,i, ".pdf"), w=15, h=10)
#     plot(pam.sill, col=colorRampPalette(brewer.pal(12,"Set3"))(i))
#     dev.off()
    
#     png(paste0(out.path,"/Parallel1_plot_k" ,i, ".png"), w=600, h=1000)
#     plot(parallel(~pam.df[, -ncol(pam.df)], pam.df, groups=clustering))
#     dev.off()
#     
#     png(paste0(out.path,"/Parallel2_plot_k" ,i, ".png"), w=1000, h=1000)
#     plot(parallel(~pam.df[, -ncol(pam.df)] | clustering, pam.df))
#     dev.off()
    
#     pdf(paste0(out.path, "/Parcoord_plot_k" ,i, ".pdf"), w=10, h=5)
#     for(j in 1:i){
#       parcoord(x=pam.x[pam.cl==j, ], col = colorRampPalette(brewer.pal(12,"Set3"))(i)[j], lty = 1, var.label = FALSE)
#     }
#     dev.off()
    
#     pdf(paste0(out.path, "/Parcoord_plot_k" ,i, ".pdf"), w=10, h=5)
#     for(j in 1:i){
#      plot(1:ncol(pam.x), rep(1, ncol(pam.x)), type="n", ylim=c(ylimL,ylimU))
#       apply(pam.x[pam.cl==j, ], 1, function(g){
#         lines(1:ncol(pam.x), g, col=colorRampPalette(brewer.pal(12,"Set3"))(i)[j])
#       }) 
#     }
#     dev.off()
    
    pam.df.sort <- pam.df[order(pam.df$clustering), ]
    
    png(paste0(out.path,"/Heatmap_plot_k" ,i, ".png"), w=1000, h=1000)
    heatmap.2axis(x=as.matrix(pam.df.sort[,-ncol(pam.df.sort)]), Rowv=NA, Colv=NA, dendrogram="none", labRow=NA, labCol=NULL, scale="none", trace="none", col=colorRampPalette(brewer.pal(9,"Blues"))(100), key=T, keysize =c(0.5, 0.5), density.info="none", margins = c(5,5), xlab = "Samples", ylab = "Clusters of DE Genes", RowSideColors=colorRampPalette(brewer.pal(12,"Set3"))(i)[pam.df.sort[, ncol(pam.df.sort)]], ColSideColors=colorRampPalette(brewer.pal(8,"Pastel2"))(length(levels(new.samps$tree_ID)))[new.samps[samps.order, "tree_ID"]]) 
    dev.off()
    
    pam.medoids <- pam.obj[[i]]$medoids
    
#     png(paste0(out.path,"/Heatmap_plot_k" ,i, "b.png"), w=1000, h=1000)
#     heatmap.2axis(x=as.matrix(pam.df.sort[pam.medoids,-ncol(pam.df.sort)]), Rowv=NA, Colv=NA, dendrogram="none", labRow=NA, labCol=NULL, scale="none", trace="none", col=colorRampPalette(brewer.pal(9,"Blues"))(100), key=T, keysize =c(0.5, 0.5), density.info="none", margins = c(5,5), xlab = "Samples", ylab = "Medoids", RowSideColors=colorRampPalette(brewer.pal(12,"Set3"))(i)[pam.df.sort[pam.medoids, ncol(pam.df.sort)]], ColSideColors=colorRampPalette(brewer.pal(8,"Pastel2"))(length(levels(new.samps$tree_ID)))[new.samps[samps.order, "tree_ID"]]) 
#     dev.off()
    
  }
  
  save(pam.obj, file=paste0(out.path,"/pam_obj_norm.Rdata"))
  
}



PAM.clustering(pam.x=d.cpm.l.n, new.samps , out.path="Clustering/PAM_log_01/", out.name="")


PAM.clustering(pam.x=d.cpm.n, new.samps , out.path="Clustering/PAM_01/", out.name="")


############################################
# hclust
######################

library(fastcluster)
library(gplots) 
library(RColorBrewer)
source("~//R_UZH//heatmap2.r") # function heatmap.2axis


HC.clustering <- function(x.obj, cl.method="ward", out.path="Clustering/HC/",  out.name=""){
  
  dir.create(out.path, showWarnings=F, recursive=T)
  
  samps.order <- row.names(new.samps[order(new.samps$drough.control, new.samps$tree_ID, new.samps$time_nr), ])
  
  x.obj <- x.obj[, samps.order]
  new.samps <- new.samps[samps.order, ]
  
  #dist.obj.g <- dist(x.obj, method = "euclidean")
#   sim.obj.g <- cor(t(x.obj), method="spearman")
#   sim.obj.s <- cor(x.obj, method="spearman")
#   
#   dist.obj.g <- as.dist(1-sim.obj.g)
#   dist.obj.s <- as.dist(1-sim.obj.s)
  
  dist.obj.g <- dist(x.obj, method = "euclidean")
  dist.obj.s <- dist(t(x.obj), method = "euclidean")
  
  #hc.obj <- hclust.vector(X=x.obj,  method="single",  metric="euclidean")
  #hc.obj.ss <- stats::hclust(d=dist.obj.s, method="complete")
  hc.obj.g <- fastcluster::hclust(d=dist.obj.g, method=cl.method)
  hc.obj.s <- fastcluster::hclust(d=dist.obj.s, method=cl.method)
  
  
  # png(paste0(out.path,"/DendoS", out.name ,".png"), w=1000, h=700)
  # plot(hc.obj.s, hang=-1, col=1:nrow(new.samps))
  # dev.off()
  
  #library(ggplot2)
  #library(ggdendro)
  
  # cols <- c("black","blue","orange","darkgreen","red","salmon")
  # cols[new.samps[hc.obj.s$order, "tree_ID"]]
  # colorRampPalette(brewer.pal(8,"Pastel2"))(length(levels(new.samps$tree_ID)))[new.samps[hc.obj.s$order, "tree_ID"]]
  
  # png(paste0(out.path,"/DendoS2", out.name ,".png"), w=700, h=1000)
  # ggdendrogram(hc.obj.ss, rotate=T, size=4, theme_dendro=F, color=colorRampPalette(brewer.pal(8,"Pastel2"))(length(levels(new.samps$tree_ID)))[new.samps[hc.obj.s$order, "tree_ID"]])
  # dev.off()
  
  dend.obj.g <- as.dendrogram(hc.obj.g) 
  dend.obj.s <- as.dendrogram(hc.obj.s) 
  
  png(paste0(out.path,"/hmG", out.name ,".png"), w=1000, h=1000)
  heatmap.2axis(x.obj, Rowv=dend.obj.g, Colv=NA, dendrogram="row", labRow=NA, labCol=NULL, scale="none", trace="none", col=colorRampPalette(brewer.pal(9,"Blues"))(100), key=T, keysize =c(0.5, 1), density.info="none", margins = c(5,5), xlab = "Samples", ylab = "DE Genes", ColSideColors=colorRampPalette(brewer.pal(8,"Pastel2"))(length(levels(new.samps$tree_ID)))[new.samps[samps.order, "tree_ID"]]) 
  dev.off()
  
  
  png(paste0(out.path,"/hmGS", out.name ,".png"), w=1000, h=1000)
  heatmap.2axis(x.obj, Rowv=dend.obj.g, Colv=dend.obj.s, dendrogram="both", labRow=NA, labCol=NULL, scale="none", trace="none", col=colorRampPalette(brewer.pal(9,"Blues"))(100), key=T, keysize =c(1, 1), density.info="none", margins = c(4,4), xlab = "Samples", ylab = "DE Genes", ColSideColors=colorRampPalette(brewer.pal(8,"Pastel2"))(length(levels(new.samps$tree_ID)))[new.samps[samps.order, "tree_ID"]]) 
  dev.off()
  
}


HC.clustering(x.obj=d.cpm.l.n, cl.method="ward", out.path="Clustering/HC_log_01/",  out.name="_ward")
HC.clustering(x.obj=d.cpm.l.n, cl.method="complete", out.path="Clustering/HC_log_01/",  out.name="_comp")
HC.clustering(x.obj=d.cpm.l.n, cl.method="average", out.path="Clustering/HC_log_01/",  out.name="_avg")
HC.clustering(x.obj=d.cpm.l.n, cl.method="centroid", out.path="Clustering/HC_log_01/",  out.name="_cent")
HC.clustering(x.obj=d.cpm.l.n, cl.method="mcquitty", out.path="Clustering/HC_log_01/",  out.name="_mcq")

HC.clustering(x.obj=d.cpm.n, cl.method="ward", out.path="Clustering/HC_01/",  out.name="_ward")
HC.clustering(x.obj=d.cpm.n, cl.method="complete", out.path="Clustering/HC_01/",  out.name="_comp")
HC.clustering(x.obj=d.cpm.n, cl.method="average", out.path="Clustering/HC_01/",  out.name="_avg")
HC.clustering(x.obj=d.cpm.n, cl.method="centroid", out.path="Clustering/HC_01/",  out.name="_cent")
HC.clustering(x.obj=d.cpm.n, cl.method="mcquitty", out.path="Clustering/HC_01/",  out.name="_mcq")



#####################################################################################################

### GO analysis with topGO

#####################################################################################################

load("Analysis_till_Sep09/Models/Models_fitting_NoControl.RData")
load("Analysis_till_Sep09/Models/Models_selecting_NoControl.RData")
load("Analysis_till_Sep09/Models/Models_clusters_NoControl.RData")

load("Analysis_till_Sep09/Models/Models_fitting_28_NoControl.RData")
load("Analysis_till_Sep09/Models/Models_selecting_28_NoControl.RData")
load("Analysis_till_Sep09/Models/Models_clusters_28_NoControl.RData")

genes <- names(gene.models)
models <- c(names(models.results), "Intercept")


gene.models.table <- data.frame(genes=genes, model=unlist(lapply(gene.models, function(gm) gm$final.model)))

gene.descr <- read.table("Data/genes_descr_control/genes_description.xls", sep="\t", stringsAsFactors=FALSE, head=TRUE)

head(gene.descr)

gene.models.table.descr <- merge(gene.models.table, gene.descr, by=1, all.x=TRUE)

head(gene.models.table.descr)

write.table(gene.models.table.descr, file="Analysis_till_Sep09/Models/Table_genes_with_models.csv", sep=";", row.names=FALSE, quote=FALSE)


write.table(gene.models.table.descr, file="Analysis_till_Sep09/Models/Table_genes_with_models28.csv", sep=";", row.names=FALSE, quote=FALSE)


MF2 <- read.table("ModelsAIC/Best_Models_fit_AIC.xls", sep="\t",header=TRUE)


gene.models.table.descr.MF2 <- merge(gene.models.table.descr, MF2, by=1, all.x=TRUE)


write.table(gene.models.table.descr.MF2, file="Analysis_till_Sep09/Models/Table_genes_with_models28_plus_AICmodels.csv", sep=";", row.names=FALSE, quote=FALSE)

###############################################################

length(clusters[["Intercept"]])

length(clusters[["Water.Potential"]])
length(clusters[["Soil.Moisture"]])
length(clusters[["Temperature"]])
length(unique(clusters.AT[["Water.Potential"]]))
length(unique(clusters.AT[["Soil.Moisture"]]))
length(unique(clusters.AT[["Temperature"]]))


length(clusters[["Water.Potential28"]])
length(clusters[["Soil.Moisture28"]])
length(clusters[["Temperature28"]])
length(unique(clusters.AT[["Water.Potential28"]]))
length(unique(clusters.AT[["Soil.Moisture28"]]))
length(unique(clusters.AT[["Temperature28"]]))


head(MolEcol.DE.genes)

pdf("Plot_Hist/Hist_clust_WP.pdf", width = 10, height = 10)
hist(MolEcol.DE.genes[ MolEcol.DE.genes[,1] %in% clusters.AT[["Water.Potential"]], 2], breaks=7, main="Water.Potential")
dev.off()

pdf("Plot_Hist/Hist_clust_SM.pdf", width = 10, height = 10)
hist(MolEcol.DE.genes[ MolEcol.DE.genes[,1] %in% clusters.AT[["Soil.Moisture"]], 2], breaks=7, main="Soil.Moisture")
dev.off()

pdf("Plot_Hist/Hist_clust_T.pdf", width = 10, height = 10)
hist(MolEcol.DE.genes[ MolEcol.DE.genes[,1] %in% clusters.AT[["Temperature"]], 2], breaks=7, main="Temperature")
dev.off()


pdf("Plot_Hist/Hist_clust_WP28.pdf", width = 10, height = 10)
hist(MolEcol.DE.genes[ MolEcol.DE.genes[,1] %in% clusters.AT[["Water.Potential28"]], 2], breaks=7, main="Water.Potential28")
dev.off()

pdf("Plot_Hist/Hist_clust_SM28.pdf", width = 10, height = 10)
hist(MolEcol.DE.genes[ MolEcol.DE.genes[,1] %in% clusters.AT[["Soil.Moisture28"]], 2], breaks=7, main="Soil.Moisture28")
dev.off()

pdf("Plot_Hist/Hist_clust_T28.pdf", width = 10, height = 10)
hist(MolEcol.DE.genes[ MolEcol.DE.genes[,1] %in% clusters.AT[["Temperature28"]], 2], breaks=7, main="Temperature28")
dev.off()

############################
### GO
############################

library(topGO)
library(org.At.tair.db)

fun.gene.sel <- function(gene.vector) {
  return(gene.vector <- ifelse(gene.vector==0, FALSE, TRUE))
}

assayed.genes <- unique(unlist(clusters.AT, use.names = FALSE))

for(m in models[c(5,6,7)]){
  
  de.genes <- unique(clusters.AT[[m]])
  
  gene.vector=as.numeric(assayed.genes %in% de.genes)
  names(gene.vector) <- assayed.genes
  names(assayed.genes) <- gene.vector
  
  for(go in c("BP","MF","CC")){
    
    sampleGOdata <- new("topGOdata", description = "Simple session", ontology = go, allGenes = gene.vector, geneSel = fun.gene.sel , nodeSize = 10, annot = annFUN.org, mapping = "org.At.tair.db")
    
    resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
    resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
    resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")
    
    allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,classicKS = resultKS, elimKS = resultKS.elim, orderBy = "elimKS", ranksOf = "elimKS", topNodes = 10)
    
    write.table(allRes, paste("GO/GO_", m , go, ".csv", sep=""), sep=";")
    
    pdf(paste("Plots_GO/GO_", m , go, ".pdf", sep=""), width = 10, height = 10)
    showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')
    dev.off()
    
  }
  
}





