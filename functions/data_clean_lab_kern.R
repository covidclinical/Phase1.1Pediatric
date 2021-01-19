data_clean_lab_kern=function(site.read, code.dict, y.scale="original", tab.obs){
site.list=site.read$SiteID
county.list=site.read$Country
combine.set=c("48065-7","48066-5")
combine.nm=paste(combine.set, collapse=":")

dat=NULL
for(site.nm in site.list){
    dat.site=read.csv(paste0(dir.dat, "Labs-", site.nm, ".csv"))
    colnames(dat.site)=tolower(colnames(dat.site))
    obs <- tab.obs$obfuscation[which(tab.obs$siteid == site.nm)]
    impute99 <- NA
    #if (length(obs) > 0){
    #  if (obs != 'none'){
    #    impute99 <- 0.5 * as.integer(as.character(obs))
    #  }
    #}
    dat.site[dat.site==-99]=impute99;dat.site[dat.site==-999]=NA;dat.site[dat.site==-Inf]=NA;dat.site[dat.site==Inf]=NA
    
    ## combine two d-dimers
    dat.add=dat.site[dat.site$loinc%in%combine.set,]
    if(dim(dat.add)[1]!=0){
    dat.add.new=dat_clean_lab_merge_ddimer_fun(dat.add)
    dat.site=rbind(dat.site, dat.add.new)
    }
    dat=rbind(dat,dat.site)
}

dat$siteid=toupper(dat$siteid)
nm.site="siteid"
nm.day="days_since_admission"
nm.loinc="loinc"
nm.common=c("num_patients", "mean_value", "stdev_value")
nm.log.common=c("num_patients", "mean_log_value", "stdev_log_value")

nm.all=paste0(nm.common,"_all")
nm.ever=paste0(nm.common,"_ever_severe")
nm.never=paste0(nm.common,"_never_severe")

nm.log.all=paste0(nm.log.common,"_all")
nm.log.ever=paste0(nm.log.common,"_ever_severe")
nm.log.never=paste0(nm.log.common,"_never_severe")

nm.col=unique(c(nm.site, nm.day,nm.loinc, nm.all, nm.ever, nm.log.all, nm.log.ever))
names(nm.common)=names(nm.all)=names(nm.ever)=names(nm.never)=
  names(nm.log.common)=names(nm.log.all)=names(nm.log.ever)=names(nm.log.never)=c("num", "mean", "stdev")

comb.lab=dat[,nm.col]
comb.lab$loinc=trimws(comb.lab$loinc, which = c("both", "left", "right"))
comb.lab=left_join(comb.lab, data.frame(code.dict), by="loinc")
comb.lab = comb.lab[!is.na(comb.lab$mean_value_all),]
nm.lab.all = unique(comb.lab$labname)
site.lab.all = unique(comb.lab$siteid); n.site = length(site.lab.all)

dat.list = as.list(nm.lab.all); names(dat.list) = nm.lab.all

for(nm.lab in nm.lab.all){
if(y.scale=="original"){nm.all.use=nm.all; nm.ever.use=nm.ever}
if(y.scale=="log"){nm.all.use=nm.log.all; nm.ever.use=nm.log.ever}

  tt0 = sort(unique(comb.lab[,nm.day][comb.lab$labname==nm.lab]))
  n.site = length(site.lab.all); 
  labsd.mat.all = labn.mat.all = labmu.mat.all = matrix(NA, nrow=length(tt0), ncol=n.site); 
  labsd.mat.ever = labn.mat.ever = labmu.mat.ever = matrix(NA, nrow=length(tt0), ncol=n.site); 
  labsd.mat.never = labn.mat.never = labmu.mat.never = matrix(NA, nrow=length(tt0), ncol=n.site); 
  
  colnames(labsd.mat.all) = colnames(labn.mat.all) = colnames(labmu.mat.all) = site.lab.all
  colnames(labsd.mat.ever) = colnames(labn.mat.ever) = colnames(labmu.mat.ever) = site.lab.all
  colnames(labsd.mat.never) = colnames(labn.mat.never) = colnames(labmu.mat.never) = site.lab.all
  
  for(site.ID in site.lab.all){
    tmpdat = comb.lab[comb.lab$siteid==site.ID&comb.lab$labname==nm.lab,]
    if(nrow(tmpdat)>0){
      n.lab = paste(range(tmpdat[,nm.all.use["num"]]),collapse="_")
      labmu.mat.all[match(tmpdat[,nm.day], tt0),site.ID] = tmpdat[,nm.all.use["mean"]]
      labsd.mat.all[match(tmpdat[,nm.day], tt0),site.ID] = tmpdat[,nm.all.use["stdev"]]
      labn.mat.all[match(tmpdat[,nm.day], tt0),site.ID] = tmpdat[,nm.all.use["num"]]
      
      labmu.mat.ever[match(tmpdat[,nm.day], tt0),site.ID] = tmpdat[,nm.ever.use["mean"]]
      labsd.mat.ever[match(tmpdat[,nm.day], tt0),site.ID] = tmpdat[,nm.ever.use["stdev"]]
      labn.mat.ever[match(tmpdat[,nm.day], tt0),site.ID] = tmpdat[,nm.ever.use["num"]]
    }
  }
  
  tmp.list=as.list(c("t0", "all", "ever", "never"))
  names(tmp.list)=c("t0", "all", "ever", "never")
  for(setting in c("all", "ever")){
    if(setting=="all"){labmu.mat=labmu.mat.all; labsd.mat=labsd.mat.all; labn.mat=labn.mat.all}
    if(setting=="ever"){labmu.mat=labmu.mat.ever; labsd.mat=labsd.mat.ever; labn.mat=labn.mat.ever}
    
    labn.mat[labn.mat <0] = NA
    labsd.mat[labsd.mat <0] = NA;
    if(y.scale=="log"){ 
      labmu.mat[labmu.mat %in%c(-99, -999, -Inf)] = NA
    }
    if(y.scale!="log"){ 
      labmu.mat[labmu.mat<0] = NA
    }  
    
    labsd.mat[labsd.mat=="NULL"] = NA; mode(labsd.mat) = "numeric"
    if(is.children==F){
      site.ICSM = Find.Strings("ICSM",colnames(labmu.mat))
      labmu.mat = cbind(labmu.mat,"ICSM" = apply((labmu.mat*labn.mat)[,site.ICSM],1,sum,na.rm=T)/apply(labn.mat[,site.ICSM],1,sum,na.rm=T))
      tmpsd = labsd.mat[,site.ICSM]; tmpmu = labmu.mat[,site.ICSM]; tmpn = labn.mat[,site.ICSM]; tmpsd[tmpn==1] = 0; 
      
      sd.comb = sapply(1:nrow(tmpsd),function(kk){sd.pool(tmpmu[kk,],tmpsd[kk,],tmpn[kk,])})
      labsd.mat = cbind(labsd.mat,"ICSM" = sd.comb)
      labn.mat = cbind(labn.mat, "ICSM" = apply(labn.mat[,site.ICSM],1,sum,na.rm=T))
      site.lab.all2 = c("ICSM", setdiff(site.lab.all, site.ICSM))
      country.lab.all2=left_join(data.frame("SiteID"=site.lab.all2), site.read, by="SiteID")[,"Country"]
      country.lab.all2[1]="ITALY"
    }else{
      site.lab.all2=site.lab.all
      country.lab.all2=left_join(data.frame("SiteID"=site.lab.all2), site.read, by="SiteID")[,"Country"]
    }
    site.lab.all2 = site.lab.all2[order(country.lab.all2)]; country.lab.all2 = country.lab.all2[order(country.lab.all2)]
    
    labmu.mat = labmu.mat[,site.lab.all2]; labn.mat = labn.mat[,site.lab.all2]; labsd.mat = labsd.mat[,site.lab.all2]; n.site2 = length(site.lab.all2)
    labmu.mat[labmu.mat==0] = NA
    labn.mat[is.na(labmu.mat)]=0
    labsd.mat[is.na(labmu.mat)]=NA
    nmax.all = apply(labn.mat,2,max, na.rm=T); nmax.all[nmax.all == -Inf] = NA; ##nmax.all = pmax(nmax.all,0)
    site.outlier = abs(apply(labmu.mat,2,median,na.rm=T)-median(labmu.mat,na.rm=T)) > 10*IQR(labmu.mat,na.rm=T); 
    site.outlier[is.na(site.outlier)] = T
    tmp.list[[setting]]$labmu.mat=labmu.mat
    tmp.list[[setting]]$labsd.mat=labsd.mat
    tmp.list[[setting]]$labn.mat=labn.mat
    tmp.list[[setting]]$nmax.all=nmax.all
    tmp.list[[setting]]$site.outlier=site.outlier
  }
  
  tmp.list[["never"]]$nmax.all=tmp.list[["all"]]$nmax.all-tmp.list[["ever"]]$nmax.all
  tmp.list[["never"]]$labn.mat=tmp.list[["all"]]$labn.mat-tmp.list[["ever"]]$labn.mat
  tmp.list[["never"]]$labmu.mat=do.call(rbind,lapply(1:dim(tmp.list[["all"]]$labmu.mat)[1], 
                                                     function(ll) mu.never(mu.all=tmp.list[["all"]]$labmu.mat[ll,], 
                                                                           n.all=tmp.list[["all"]]$labn.mat[ll,], 
                                                                           mu.ever=tmp.list[["ever"]]$labmu.mat[ll,],
                                                                           n.ever=tmp.list[["ever"]]$labn.mat[ll,])))
  
  tmp.list[["never"]]$labsd.mat=do.call(rbind,
                                        lapply(1:dim(tmp.list[["all"]]$labsd.mat)[1], 
                                               function(ll)
                                                 sd.subtract(mu0=tmp.list[["all"]]$labmu.mat[ll,], 
                                                             mu1=tmp.list[["ever"]]$labmu.mat[ll,], 
                                                             mu2=tmp.list[["never"]]$labmu.mat[ll,], 
                                                             n0= tmp.list[["all"]]$labn.mat[ll,], 
                                                             n1= tmp.list[["ever"]]$labn.mat[ll,], 
                                                             n2= tmp.list[["never"]]$labn.mat[ll,], 
                                                             sd0=tmp.list[["all"]]$labsd.mat[ll,], 
                                                             sd1=tmp.list[["ever"]]$labsd.mat[ll,])))
  tmp.list[["never"]]$labsd.mat[tmp.list[["never"]]$labsd.mat==Inf]=NA
  #check for NAs:
  #mu0=1.173
  #mu1=1.009
  #mu2=1.337000
  #n0=30
  #n1=15
  #n2=15
  #sd0=0.267
  #sd1=0.448
  #sd.subtract(mu0,mu1,mu2,n0,n1,n2,sd0,sd1)
  ## correct for MGB and BIDMC
  if(is.children==F){
    for(site.correct in c("MGB", "BIDMC")){
      tmp.correct=read.csv(paste0("../../Phase1.1/enchance/Labs_enhance_", site.correct,".csv"))
      if(nm.lab=="D-dimer"){
      tmp.correct1=tmp.correct[as.character(tmp.correct$loinc)==as.character(code.dict[which(code.dict$labname=="D-dimer (FEU)"),"loinc"]),]
      tmp.correct2=tmp.correct[as.character(tmp.correct$loinc)==as.character(code.dict[which(code.dict$labname=="D-dimer (DDU)"),"loinc"]),]
      dat.add=rbind(tmp.correct1, tmp.correct2)
      dat.add.new=dat_clean_lab_merge_ddimer_fun(dat.add)
      tmp.correct=dat.add.new
      }
      tmp.correct=tmp.correct[as.character(tmp.correct$loinc)==as.character(code.dict[which(code.dict$labname==nm.lab),"loinc"]),c("days_since_admission", gsub("_ever_", "_never_", nm.ever.use))]
      if(dim(tmp.correct)[1]!=0){
        tmp.correct[tmp.correct==-99|tmp.correct==-999|tmp.correct==-Inf|tmp.correct==Inf]=NA
        tmp.correct=left_join(data.frame(days_since_admission=tt0), tmp.correct, by="days_since_admission")
        tmp.list[["never"]]$nmax.all[site.correct]=max(tmp.correct[,gsub("_ever_", "_never_", nm.ever.use)[1]],na.rm=T)
        tmp.list[["never"]]$labn.mat[,site.correct]=tmp.correct[,gsub("_ever_", "_never_", nm.ever.use)[1]]
        tmp.list[["never"]]$labmu.mat[,site.correct]=tmp.correct[,gsub("_ever_", "_never_", nm.ever.use)[2]]
        tmp.list[["never"]]$labsd.mat[,site.correct]=tmp.correct[,gsub("_ever_", "_never_", nm.ever.use)[3]]
      }
    }
  }
  
  if(is.children==T){
    for(site.correct in c("MGBPED")){
      tmp.correct=read.csv(paste0("../../Phase1.1/enchance/Labs_enhanced_", site.correct,".csv"))
      if(nm.lab=="D-dimer"){
        tmp.correct1=tmp.correct[as.character(tmp.correct$loinc)==as.character(code.dict[which(code.dict$labname=="D-dimer (FEU)"),"loinc"]),]
        tmp.correct2=tmp.correct[as.character(tmp.correct$loinc)==as.character(code.dict[which(code.dict$labname=="D-dimer (DDU)"),"loinc"]),]
        dat.add=rbind(tmp.correct1, tmp.correct2)
        dat.add.new=dat_clean_lab_merge_ddimer_fun(dat.add)
        tmp.correct=dat.add.new
      }
      tmp.correct=tmp.correct[as.character(tmp.correct$loinc)==as.character(code.dict[which(code.dict$labname==nm.lab),"loinc"]),c("days_since_admission", gsub("_ever_", "_never_", nm.ever.use))]
      if(dim(tmp.correct)[1]!=0){
        tmp.correct[tmp.correct==-99|tmp.correct==-999|tmp.correct==-Inf|tmp.correct==Inf]=NA
        tmp.correct=left_join(data.frame(days_since_admission=tt0), tmp.correct, by="days_since_admission")
        tmp.list[["never"]]$nmax.all[site.correct]=max(tmp.correct[,gsub("_ever_", "_never_", nm.ever.use)[1]],na.rm=T)
        tmp.list[["never"]]$labn.mat[,site.correct]=tmp.correct[,gsub("_ever_", "_never_", nm.ever.use)[1]]
        tmp.list[["never"]]$labmu.mat[,site.correct]=tmp.correct[,gsub("_ever_", "_never_", nm.ever.use)[2]]
        tmp.list[["never"]]$labsd.mat[,site.correct]=tmp.correct[,gsub("_ever_", "_never_", nm.ever.use)[3]]
      }
    }
  }
  
  
  dat.list[[nm.lab]] = as.list(c("t0", "all", "ever", "never")); names(dat.list[[nm.lab]]) =c("t0", "all", "ever", "never")
  dat.list[[nm.lab]] = list("t0"=tt0, "all"=tmp.list[["all"]], "ever"=tmp.list[["ever"]], "never"=tmp.list[["never"]])
}
names(dat.list)=nm.lab.all
dat.list
}



