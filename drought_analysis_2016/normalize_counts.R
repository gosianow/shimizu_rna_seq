normalize.counts <- function(counts, norm.method="norm", c=0, d=1){
  
  t(apply(counts, 1, function(g){
    
    if(norm.method=="norm"){
      m <- mean(g, na.rm=T)
      sd <- sd(g, na.rm=T) 
      return((g-m)/sd)
    }else if(norm.method=="mean0"){
      m <- mean(g, na.rm=T)
      return(g-m)
    }else if(norm.method=="01"){
      return((g-min(g,na.rm=T))/(max(g,na.rm=T)-min(g,na.rm=T)))
    }
    else if(norm.method=="range"){
      return(c*(1-(g-min(g,na.rm=T))/(max(g,na.rm=T)-min(g,na.rm=T)))+d*((g-min(g,na.rm=T))/(max(g,na.rm=T)-min(g,na.rm=T))))
    }
  }))
}

