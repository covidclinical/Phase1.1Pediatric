table_lab_meta=function(dat.clean.lab, day=0, nm.lab.all, y.scale.all, setting = 'diff',
                        method = 'random effect'){
  tab=NULL
    for(y.scale in y.scale.all){
      for(nm.lab in nm.lab.all){
        tmp=table_lab_meta_kern(dat=dat.clean.lab, nm.lab, day, y.scale=y.scale, 
                                setting = setting, method = method)
        tmp=c(nm.lab, y.scale, tmp)
        
        if (method == 'pooled'){
          names(tmp)=c("labname","scale", "est", "se","lb", "ub","pvalue", 'n0')
        }else{
          names(tmp)=c("labname","scale", "est", "se","lb", "ub","pvalue")
        }
        tab=rbind(tab, tmp)
      }
  }
  tab
  }


