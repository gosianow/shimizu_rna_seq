### FUN that interpolates and calculates rolled means 

library(zoo)
library(stringr)

# intrp.x=sm$time_nr
# intrp.y=sm$DE8266
# intrp.points=all.days[,2]
# table1=new.samps[new.samps$tree_ID == 8266,]
# col.values="Soil.Moisture"
# plots.path="Plots_Interpolate"
# plot.name="DE8266"
# roll.days=c(7, 14, 28)
# ylim=c(0,1)
# color=6

interpolate.and.roll <- function(intrp.x, intrp.y, intrp.points, table1, roll.days=c(7, 14, 28), col.values="Soil.Moisture", plots.path="Plots_Interpolate", plot.name="", ylim=c(0,1), month.days=month.days, all.days.short=all.days.short, color=6){
  #table1 has to have "time_nr" column
   
  interp.intrp.points <- approx(intrp.x[!is.na(intrp.y)], intrp.y[!is.na(intrp.y)], xout=intrp.points)
  
  table2<-  data.frame(time_nr=interp.intrp.points$x[!is.na(interp.intrp.points$y)], interp=interp.intrp.points$y[!is.na(interp.intrp.points$y)])
  
  colnames(table2) <- c("time_nr", col.values)
  
  ### rollmean over 14 and 28 days
  for(r in roll.days){
    
    col.values.r <- paste(col.values, r , sep="")
    table2[,col.values.r] <- NA
    table2[r:nrow(table2), col.values.r] <- rollmean(table2[,col.values], r)
     
  }
  
  table1 <-  merge(table1, table2, by="time_nr", all.x=T)
  
  dir.create(plots.path, recursive=T, showWarnings=FALSE)
  
#   pdf( paste(plots.path, "/InterpRoll_", str_replace_all(paste(col.values, plot.name, sep=""),"[[:punct:]]", "_"),".pdf", sep=""), width = 10, height = 5)
#   plot(intrp.x, intrp.y, col=1, pch=20, ylab=col.values, xlab="Time", ylim=ylim, xlim=c(min(all.days.short[,2]), max(all.days.short[,2])), xaxt = "n")
#   axis(side=1, at=month.days[,2], labels=month.days[,1])
#   abline(v=table1[,"time_nr"], col="grey")
#   lines(table2[,"time_nr"], table2[,col.values], col=color)
#   for(r in roll.days){
#     col.values.r <- paste(col.values, r , sep="")
#     lines(table2[,"time_nr"], table2[,col.values.r], col=color,  lty=(which(roll.days==r)+1))   
#   }
#   dev.off()
  
  
  invisible(table1)
  
}

