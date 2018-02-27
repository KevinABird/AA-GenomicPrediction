library("BiocInstaller")
#biocLite("multtest")
library("multtest")
library("gplots")
library("LDheatmap")
library("genetics")
library("EMMREML")
library("compiler")
library("scatterplot3d")
source("http://www.zzlab.net/GAPIT/emma.txt")
source("C:\\Users/Kevin/Dropbox/GAPITFunctionsEdit.txt")

MyY=read.csv("C:\\Users/Kevin/Documents/Blup-transformedforGWAS.csv",head=T)
myGAPIT<-GAPIT(
  Y=as.data.frame(MyY[,c(1,2)]),
  file.GD="C:\\Users/Kevin/Documents/At-GenomicSelection-Ruthie/AASubsetGeno/AASubsetGenoGD",
  file.GM="C:\\Users/Kevin/Documents/At-GenomicSelection-Ruthie/AASubsetGeno/AASubsetGenoGM",
  file.Ext.GM = "txt",
  file.from=1,
  file.to=4,
  PCA.total=6,
  group.from=324,
  group.to=324
  
  
)
