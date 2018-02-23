#Load pacakges for GAPIT

library("BiocInstaller")
#biocLite("multtest")
library("multtest")
library("gplots")
library("LDheatmap")
library("genetics")
library("EMMREML")
library("compiler")
library("scatterplot3d")


#Load EMMA & GAPIT Source code


source("http://www.zzlab.net/GAPIT/emma.txt")
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")

random_geno_list<-list()
for(n in seq(132)) {
random_geno_list[[n]] <- sprintf('randomSubset_%s', n)
}
for(file in random_geno_list) {
dir.create((file))
setwd(file.path(getwd(),file))
 myGAPIT<-GAPIT(
  file.GD=sprintf("/home/birdk/data/AARandomGeneration/%s_GD",file),
  file.GM=sprintf("/home/birdk/data/AARandomGeneration/%s_GM",file),
  file.Ext.GD="txt",
  file.Ext.GM = "txt",
  file.from=1,
  file.to=4,
  PCA.total=6
  )
  setwd("..")
}
