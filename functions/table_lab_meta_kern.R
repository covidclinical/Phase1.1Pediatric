se.sum <- function(n1,n2,sd1,sd2,mean1,mean2){
  set_2 <- which(is.na(sd2) & ! is.na(sd1))
  set_1 <- which(is.na(sd1) & ! is.na(sd2))
  if (length(set_1) > 0){
    n1[set_1] <- 0
    sd1[set_1] <- 0
    mean1[set_1] <- 0
  }
  if (length(set_2) > 0){
    n2[set_2] <- 0
    sd2[set_2] <- 0
    mean2[set_2] <- 0
  }
  n0 = n1 + n2
  mean0 = (mean1 * n1 + mean2 * n2) / n0
  v1=sd1^2 * (n1 - 1) + mean1^2 * n1
  v2=sd2^2 * (n2 - 1) + mean2^2 * n2
  se.sum=sqrt((v1 + v2 - n0 * mean0^2) / (n0 - 1) / n0)
  return(list(mean = mean0, se = se.sum))
}



table_lab_meta_kern=function(dat.clean.lab, nm.lab, day=0, y.scale="original", setting="diff",
                             method = 'random effect'){
  lab.dat.list=dat.clean.lab[[y.scale]]
  n.site2=dim(lab.dat.list[[1]][["all"]]$labmu.mat)[2]
  site.lab.all2=colnames(lab.dat.list[[1]][["all"]]$labmu.mat)
  tt0=lab.dat.list[[nm.lab]][["t0"]]
  id.keep=which(tt0==day)
  
  if (setting == 'sum'){
    nmax.sum=lab.dat.list[[nm.lab]][["ever"]]$nmax.all+lab.dat.list[[nm.lab]][["never"]]$nmax.all
    
    n.ever.mat <- lab.dat.list[[nm.lab]][["ever"]]$labn.mat[id.keep,]
    n.never.mat <- lab.dat.list[[nm.lab]][["never"]]$labn.mat[id.keep,]
    n.ever.mat[is.na(n.ever.mat)] <- 0
    n.never.mat[is.na(n.never.mat)] <- 0
    n.mat <- n.ever.mat + n.never.mat
    
    #rowSums(n.ever.mat, na.rm = T)
   
    nmax.all = sum(nmax.sum, na.rm = T)
    labn.mat = lab.dat.list[[nm.lab]][["ever"]]$labn.mat+lab.dat.list[[nm.lab]][["never"]]$labn.mat
    n.collect = sum(labn.mat, na.rm = T)
    
    labmu.mat = do.call(rbind,lapply(1:dim(lab.dat.list[[nm.lab]][["ever"]]$labn.mat)[1], 
                                     function(ll) se.sum(lab.dat.list[[nm.lab]][["ever"]]$labn.mat[ll,],
                                                         lab.dat.list[[nm.lab]][["never"]]$labn.mat[ll,],
                                                         lab.dat.list[[nm.lab]][["ever"]]$labsd.mat[ll,],
                                                         lab.dat.list[[nm.lab]][["never"]]$labsd.mat[ll,],
                                                         lab.dat.list[[nm.lab]][["ever"]]$labmu.mat[ll,],
                                                         lab.dat.list[[nm.lab]][["never"]]$labmu.mat[ll,])$mean))
    
    labse.mat = do.call(rbind,lapply(1:dim(labmu.mat)[1], 
                                     function(ll) se.sum(lab.dat.list[[nm.lab]][["ever"]]$labn.mat[ll,],
                                                         lab.dat.list[[nm.lab]][["never"]]$labn.mat[ll,],
                                                         lab.dat.list[[nm.lab]][["ever"]]$labsd.mat[ll,],
                                                         lab.dat.list[[nm.lab]][["never"]]$labsd.mat[ll,],
                                                         lab.dat.list[[nm.lab]][["ever"]]$labmu.mat[ll,],
                                                         lab.dat.list[[nm.lab]][["never"]]$labmu.mat[ll,])$se))
    
    mu.plot=labmu.mat; n.plot = labn.mat; se.plot=labse.mat
    mu.plot[mu.plot==Inf|mu.plot==-Inf] = NA; 
    n.plot[mu.plot==Inf|mu.plot==-Inf] = NA
    se.plot[mu.plot==Inf|mu.plot==-Inf] = NA
    
    mu.table=mu.plot[id.keep,]
    n.table=n.plot[id.keep,]
    se.table=se.plot[id.keep,]
    
  }
  
  if(setting=="diff"){
    nmax.sum=lab.dat.list[[nm.lab]][["ever"]]$nmax.all+lab.dat.list[[nm.lab]][["never"]]$nmax.all
    nmax.all=paste(lab.dat.list[[nm.lab]][["ever"]]$nmax.all,lab.dat.list[[nm.lab]][["never"]]$nmax.all,sep=":")
    labmu.mat = lab.dat.list[[nm.lab]][["ever"]]$labmu.mat-lab.dat.list[[nm.lab]][["never"]]$labmu.mat
    labn.ever.mat = lab.dat.list[[nm.lab]][["ever"]]$labn.mat
    labn.never.mat=lab.dat.list[[nm.lab]][["never"]]$labn.mat
    labse.mat = do.call(rbind,lapply(1:dim(labmu.mat)[1], function(ll) se.diff(lab.dat.list[[nm.lab]][["ever"]]$labn.mat[ll,],lab.dat.list[[nm.lab]][["never"]]$labn.mat[ll,],lab.dat.list[[nm.lab]][["ever"]]$labsd.mat[ll,],lab.dat.list[[nm.lab]][["never"]]$labsd.mat[ll,])))

    mu.table=labmu.mat[id.keep,];
    n.table.print = paste(labn.ever.mat[id.keep,], labn.never.mat[id.keep,],sep=":")
    n.table= labn.ever.mat[id.keep,]+labn.never.mat[id.keep,]
    }
  if(setting == "ever" | setting == "never"){
    labmu.mat = lab.dat.list[[nm.lab]][[setting]]$labmu.mat
    labsd.mat = lab.dat.list[[nm.lab]][[setting]]$labsd.mat
    labn.mat = lab.dat.list[[nm.lab]][[setting]]$labn.mat
    labse.mat = lab.dat.list[[nm.lab]][[setting]]$labsd.mat/sqrt(labn.mat)
    
    mu.table=labmu.mat[id.keep,];
    n.table.print = labn.mat[id.keep,]
    n.table= n.table.print
  }
  
  if (setting!="sum"){
    se.table=labse.mat[id.keep,]
  }
  
  if (method == 'random effect'){
    mu.table[mu.table==Inf|mu.table==-Inf|n.table<3|is.na(se.table)|se.table<=0] = NA; 
    n.table[is.na(mu.table)] = NA
    se.table[is.na(mu.table)] = NA
    se.table[se.table==0] = NA
    
    junk=rma(yi=mu.table, sei=se.table, method="DL")
    tab=data.frame(junk[c("b", "se", "ci.lb", "ci.ub", "pval")])
    names(tab)=c("est", "se", "lb", "ub", "pvalue")
  }
  
  if (method == 'pooled'){
    
    mu.table[mu.table==Inf|mu.table==-Inf|is.na(n.table)|n.table<1|is.na(se.table)|se.table<=0] = NA; 
    n.table[is.na(mu.table)] = NA
    se.table[is.na(mu.table)] = NA
    se.table[se.table==0] = NA
    
    n0 = sum(n.table, na.rm=T)
    mu0 = sum(mu.table * n.table, na.rm=T) / n0
    sdi <- se.table * sqrt(n.table)
    sd0 <- sqrt((sum(sdi^2 * (n0 - 1) + mu.table^2 * n.table, na.rm = T) - n0 * mu0^2) / (n0 - 1))
    se0 <- sd0 / sqrt(n0)
    tab=c(mu0, se0, mu0 - qnorm(0.975) * se0, mu0 + qnorm(0.975) * se0, 2 * pnorm(- abs(mu0) / se0), n0)
    tab <- as.list(tab)
    names(tab)=c("est", "se", "lb", "ub", "pvalue", 'n0')
  }
  
  if (method == 'fixed effect'){
    
    id.keep = which(is.na(mu.table)!=1 & n.table>=3 & se.table != Inf & is.na(se.table)!= 1)

    junk=rma(yi=mu.table[id.keep], sei=se.table[id.keep], method="FE")
    tab=data.frame(junk[c("b", "se", "ci.lb", "ci.ub", "pval")])
    names(tab)=c("est", "se", "lb", "ub", "pvalue")
  }
  
  tab
}

