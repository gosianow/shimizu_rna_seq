### FUN that interpolates and calculates rolled means 

library(zoo)


# intrp.x=sm$time_nr
# intrp.y=sm$DE8266
# intrp.points=all.days[,2]
# table1=new.samps[new.samps$tree_ID == 8266,]
# col.values="Soil.Moisture"
# plots.path="Plots_Interpolate"
# plot.name="DE8266"
# roll.days=c(14, 28)
# ylim=c(0,1)
# color=6

interpolate.and.roll <- function(intrp.x, intrp.y, intrp.points, table1, roll.days=c(7, 14, 28), col.values="Soil.Moisture", plots.path="Plots_Samples_Interpolate", plot.name="", ylim=c(0,1), month.days, all.days.short, color=6){
  
  #table1 has to have "time_nr" column
   
  intrp.points.tmp <- stats::approx(intrp.x[!is.na(intrp.y)], intrp.y[!is.na(intrp.y)], xout=intrp.points)
  
  table2 <-  data.frame(time_nr=intrp.points.tmp$x[!is.na(intrp.points.tmp$y)], interp=intrp.points.tmp$y[!is.na(intrp.points.tmp$y)])
  
  colnames(table2) <- c("time_nr", col.values)
  
  ### rollmean over 14 and 28 days
  for(r in roll.days){
    
    col.values.r <- paste(col.values, r , sep="")
    table2[,col.values.r] <- NA
    table2[r:nrow(table2), col.values.r] <- zoo::rollmean(table2[,col.values], r)
     
  }
  
  table1 <-  merge(table1, table2, by="time_nr", all.x=T)
  
  dir.create(plots.path, recursive=T, showWarnings=FALSE)
  
  pdf( paste0(plots.path, "/InterpRoll_", plot.name,".pdf", sep="" ) , width = 10, height = 5)
	  for(r in roll.days){
  plot(intrp.x, intrp.y, col=1, pch=20, ylab=col.values, xlab="Time", ylim=ylim, xlim=c(min(all.days.short[,2]), max(all.days.short[,2])), xaxt = "n", cex.lab=1.5, las=1, main=paste0("Roll over ", r, " days"))
  axis(side=1, at=month.days[,2], labels=month.days[,1])
  abline(v=table1[,"time_nr"], col="grey")
  lines(table2[,"time_nr"], table2[,col.values], col=color, lwd=1, lty=3)
	

    col.values.r <- paste(col.values, r , sep="")
    # lines(table2[,"time_nr"], table2[,col.values.r], col=color,  lty=(which(roll.days==r)+1), lwd=3)
		lines(table2[,"time_nr"], table2[,col.values.r], col=color, lty=1, lwd=4)
  }
	
    dev.off()
  
  invisible(table1)
  
}

