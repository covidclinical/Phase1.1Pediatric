data_clean=function(dir.dat, type, code.dict=NULL, site.register, site.ped, site.rm, is.children){
  site.read=site_kern(dir.dat, site.register, site.ped, type=type, site.rm, is.children)
  tab.obs <- read.table("~/111 Dropbox/Molei Liu/Phase1.1_Analysis_Data/from_Molei/input/pediatric_obfuscation.txt",
                        header = T)
  if(type=="Medications"){dat.clean=data_clean_med_kern(site.read, tab.obs)}
  if(type=="Labs"){
  lab.original=data_clean_lab_kern(site.read, code.dict, y.scale="original", tab.obs)
  lab.log=data_clean_lab_kern(site.read, code.dict, y.scale="log", tab.obs)
  dat.clean=list("original"=lab.original, "log"=lab.log)}
  dat.clean
}
