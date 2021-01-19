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


plot_lab_meta_kern_final=function(dat.clean.lab, nm.lab, lab.range, lab.qc, daymax, 
                                  nm.day="days_since_admission", y.scale, 
                                  method = 'random effect'){
#ylim.all=read.csv("xx.csv")

lab.dat.list=dat.clean.lab[[y.scale]]
n.site2=dim(lab.dat.list[[1]][["all"]]$labmu.mat)[2]
site.lab.all2=colnames(lab.dat.list[[1]][["all"]]$labmu.mat)
tt0=lab.dat.list[[nm.lab]][["t0"]]
id.keep=which(tt0>=0 & tt0<=daymax)
tt0=tt0[id.keep]

mycolors = rainbow(n.site2)
mypch = c(t(cbind(1:n.site2,1:n.site2+12)))[1:n.site2]

setting.all=c("sum")
mu.plot.list=sd.plot.list=se.plot.list=n.plot.list=tmpci.rma.list=as.list(setting.all)
names(mu.plot.list) = names(sd.plot.list)=names(n.plot.list)=names(tmpci.rma.list)=setting.all

for(setting in setting.all){
if (setting == 'sum'){
  lab.dat.list[[nm.lab]][["ever"]]$nmax.all[which(lab.dat.list[[nm.lab]][["ever"]]$nmax.all < 0)] <- 0
  lab.dat.list[[nm.lab]][["never"]]$nmax.all[which(lab.dat.list[[nm.lab]][["never"]]$nmax.all < 0)] <- 0
  nmax.sum=lab.dat.list[[nm.lab]][["ever"]]$nmax.all+lab.dat.list[[nm.lab]][["never"]]$nmax.all
  print(nmax.sum)
  
  n.ever.mat <- lab.dat.list[[nm.lab]][["ever"]]$labn.mat[id.keep,]
  n.never.mat <- lab.dat.list[[nm.lab]][["never"]]$labn.mat[id.keep,]
  n.ever.mat[is.na(n.ever.mat)] <- 0
  n.never.mat[is.na(n.never.mat)] <- 0
  n.mat <- n.ever.mat + n.never.mat
  
  #rowSums(n.ever.mat, na.rm = T)
  nmin <- unlist(lapply(1:dim(n.mat)[2], function(ll) min(n.mat[,ll], na.rm = T)))
  nmin.all <- sum(nmin, na.rm = T)
  
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
  
  mu.plot=mu.plot[id.keep,]
  n.plot=n.plot[id.keep,]
  se.plot=se.plot[id.keep,]
  
}
  
  mu.plot.list[[setting]]=mu.plot
  n.plot.list[[setting]]=n.plot
  se.plot.list[[setting]]=se.plot
  tmpci.rma.list[[setting]]=data.frame(do.call(rbind,(lapply(1:dim(mu.plot)[1], 
                                                             function(ll){
  yi=mu.plot[ll,]
  sei=se.plot[ll,]
  ni=n.plot[ll,]
  yi[sei<=0]=NA
  sdi=sei * sqrt(ni)
    
  if (method == 'random effect'){
    id.keep=which(is.na(yi)!=1 & ni>=3 & sei != Inf & is.na(sei)!= 1)
    return(tryCatch(
      as.numeric(rma(yi=yi[id.keep], sei=sei[id.keep], method="DL")[c("b", "se", "ci.lb", "ci.ub")]), 
      error=function(e) rep(NA, 4)))
  }
  if (method == 'fixed effect'){
    id.keep = which(is.na(yi)!=1 & ni>=3 & sei != Inf & is.na(sei)!= 1)
    return(tryCatch(
      as.numeric(rma(yi = yi[id.keep], sei=sei[id.keep], method="FE")[c("b", "se", "ci.lb", "ci.ub")]), 
      error=function(e) rep(NA, 4)))
  }
  
  if (method == 'pooled'){
    
    id.keep=which(is.na(yi)!=1 & ni>=1 & sei != Inf & is.na(sei)!= 1)
    n0 = sum(ni[id.keep],na.rm=T)
    mu0 = sum(yi[id.keep] * ni[id.keep], na.rm=T) / n0
    sdi.keep <- sei[id.keep] * sqrt(ni[id.keep])
    sd0 <- sqrt((sum(sdi.keep^2 * (ni[id.keep] - 1) + 
                       yi[id.keep]^2 * ni[id.keep], na.rm = T) - n0 * mu0^2) / (n0 - 1))
    se0 <- sd0 / sqrt(n0)
    if (n0 >= 2){
      return(c(mu0, se0, mu0 - 1.96 * se0, mu0 + 1.96 * se0))
    }else{
      return(rep(NA, 4))
    }
  }
    
}))))
colnames(tmpci.rma.list[[setting]]) = c("mean","se","ci_95L","ci_95U")
}









ylim=range(tmpci.rma.list[["sum"]][,c("mean", "ci_95L", "ci_95U")], na.rm=T)
nm.lab.print=nm.lab
unit.lab=lab.range[lab.range$Name==nm.lab,"Units"]
nm.lab.print[nm.lab.print=="C-reactive protein (CRP) (Normal Sensitivity)"]="C-reactive Protein (CRP)"
nm.lab.print[nm.lab.print=="albumin"]="Albumin"
nm.lab.print[nm.lab.print=="alanine aminotransferase (ALT)"]="Alanine Aminotransferase (ALT)"
nm.lab.print[nm.lab.print=="aspartate aminotransferase (AST)"]="Aspartate Aminotransferase (AST)"
nm.lab.print[nm.lab.print=="lactate dehydrogenase (LDH)"]="Lactate Dehydrogenase (LDH)"
nm.lab.print[nm.lab.print=="neutrophil count"]="Neutrophil Count"
nm.lab.print[nm.lab.print=="cardiac troponin (High Sensitivity)"]="Cardiac Troponin (High Sensitivity)"
nm.lab.print[nm.lab.print=="creatinine"]="Creatinine"
nm.lab.print[nm.lab.print=="lymphocyte count"]="Lymphocyte Count"
nm.lab.print[nm.lab.print=="procalcitonin"]="Procalcitonin"
nm.lab.print[nm.lab.print=="prothrombin time (PT)"]="Prothrombin Time (PT)"
nm.lab.print[nm.lab.print=="total bilirubin"]="Total Bilirubin"
#nm.lab.print[nm.lab.print=="white blood cell count (Leukocytes)"]="White Blood Cell Count (Leukocytes)"
nm.lab.print[nm.lab.print=="white blood cell count (Leukocytes)"]="White Blood Cell Count"


if (y.scale == 'original'){
  nm.lab.print[nm.lab.print=="D-dimer"]="D-dimer"
}
if (y.scale == 'log'){
  nm.lab.print[nm.lab.print=="D-dimer"]="log(D-dimer)"
}
if (nm.lab == 'D-dimer'){
  unit.lab <- 'ng/mL'
}


tmpci.rma.sum = data.frame("Lab"=nm.lab, "days_since_positive" = tt0, tmpci.rma.list[["sum"]])
file.nm = paste0("output/Figure_2f_lab_plot_meta",
                 "_daymax", daymax,ifelse(is.children,"_children", ""), 
                 ifelse(mater, "", "_without_maternity"), "_", method, 
                 "_", nm.lab, "_overall", ".jpg")
jpeg(file.nm, width = 500, height = 500)
add_lines <- function(dat.frame, col){
  lines(dat.frame[,"days_since_positive"], dat.frame[,"mean"], 
        col= col, lwd=2, pch=15, cex=1.5)
  piece_lst <- vector("list", 1)
  piece_lst[[1]] <- c(1000)
  for (j in 1:31){
    if (! is.na(dat.frame[j,"mean"])){
      piece_lst[[length(piece_lst)]] <- c(piece_lst[[length(piece_lst)]], j)
    }else{
      piece_lst[[length(piece_lst) + 1]] <- c(1000)
    }
  }
  for (k in 1:length(piece_lst)){
    if (length(piece_lst[[k]]) >= 3){
      set_use <- setdiff(piece_lst[[k]], 1000)
      polygon(c(dat.frame[set_use,"days_since_positive"],
                rev(dat.frame[set_use,"days_since_positive"])),
              c(dat.frame[set_use,"ci_95L"], rev(dat.frame[set_use,"ci_95U"])),
              col=t_col(col, percent = 80), border=NA)
    }
  }
}

plot(tmpci.rma.sum[,"days_since_positive"], tmpci.rma.sum[,"mean"], type="n", 
     ylim = ylim, ylab = unit.lab, xlab="Days Since Admission",
     main=paste0(nm.lab.print, ' (', nmin.all, ',', nmax.all, ')'),
     cex = 1.75, cex.lab = 1.75, cex.axis = 1.75, cex.main = 2)
add_lines(tmpci.rma.sum, "orchid4")
dev.off()

return(c(n.collect, nmax.all))

}






