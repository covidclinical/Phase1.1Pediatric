plot_lab_meta_final=function(dat.clean.lab, lab.range, lab.qc, daymax, 
                             nm.lab.all, y.scale.all, method = 'random effect'){
  vec <- c(0,0)
  for(nm.lab in nm.lab.all){
    for(y.scale in y.scale.all){
      count.vec <- plot_lab_meta_kern_final(dat=dat.clean.lab, nm.lab, lab.range, lab.qc, 
                                            daymax, y.scale=y.scale, method = method)
      vec[1] <- vec[1] + count.vec[1]
      vec[2] <- max(vec[2], count.vec[2], na.rm = T)
    }
  }
  return(vec)
}

