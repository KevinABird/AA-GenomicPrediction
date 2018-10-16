library("multtest")
library("gplots")
library("LDheatmap")
library("genetics")
library("EMMREML")
library("compiler")
library("scatterplot3d")

source("http://www.zzlab.net/GAPIT/emma.txt")
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")

for(number in c(1:50)){
AAabs <-read.csv("AfterBlupsAbsLevels04222014.csv",head=T)
myKI <-read.csv("/home/birdk/data/AA_Subset/GAPIT.Kin.VanRaden.csv", head=FALSE)



t=100 #total replicates
s=1/3 #sample of inference, e.g. set it to 1/5 for five fold cross validation
Y.raw=as.data.frame(AAabs[,c(1,2)]) #choose a trait
Y.raw=Y.raw[!is.na(Y.raw[,2]),] #Remove missing data
n=nrow(Y.raw)
n.missing=round(n*s)
storage.ref=matrix(NA,t,1)
storage.inf=matrix(NA,t,1)
#Loop on replicates
for(rep in 1:t){
#Set missing data
sample.missing=sample(1:n,n.missing)
if(n.missing>0){ Y0=Y.raw[-sample.missing,]
}else{Y0=Y.raw}
#Prediction
myGAPIT <- GAPIT(
Y=Y0,
KI=myKI,
group.from=324,
group.to=324,
group.by=10,
kinship.cluster=c("average"),
kinship.group=c("Mean")
)
prediction=myGAPIT$Pred
#Separate reference (with phenotype) and inference (without phenotype)
prediction.ref=prediction[prediction[,3]=="1",]
prediction.inf=prediction[prediction[,3]=="2",]
#Merge prediction with original Y
YP.ref <- merge(Y.raw, prediction.ref, by.x = "Taxa", by.y = "Taxa")
YP.inf <- merge(Y.raw, prediction.inf, by.x = "Taxa", by.y = "Taxa")
#Calculate correlation and store them
r.ref=cor(as.numeric(as.vector(YP.ref[,2])),as.numeric(as.vector(YP.ref[,6]) ))
r.inf=cor(as.numeric(as.vector(YP.inf[,2])),as.numeric(as.vector(YP.inf[,6]) ))
storage.ref[rep,1]=r.ref
storage.inf[rep,1]=r.inf
}#End of for (rep in 1:t)
storage=cbind(storage.ref,storage.inf)
colnames(storage)=c("Reference","Inference")
write.table(storage, file=sprintf("AASubset%s_%s_GAPIT.Cross.Validation.txt",number,colnames(AAabs)[2]), quote = FALSE, sep = "\t", row.names = TRUE,col.names = NA)
}
