library("foreach")
library("readr")
AASubsetGM1<-read_delim("AASubsetGenoGM1.txt",delim="\t")
AASubsetGM2<-read_delim("AASubsetGenoGM2.txt",delim="\t")
AASubsetGM3<-read_delim("AASubsetGenoGM3.txt",delim="\t")
AASubsetGM4<-read_delim("AASubsetGenoGM4.txt",delim="\t")
AASubsetGD1<-read_delim("AASubsetGenoGD1.txt",delim="\t")
AASubsetGD2<-read_delim("AASubsetGenoGD2.txt",delim="\t")
AASubsetGD3<-read_delim("AASubsetGenoGD3.txt",delim="\t")
AASubsetGD4<-read_delim("AASubsetGenoGD4.txt",delim="\t")
IMPORTANT_Call_Method_75_GM1<-read_delim("IMPORTANT_Call_Method_75_GM1.txt",delim="\t")
IMPORTANT_Call_Method_75_GM2<-read_delim("IMPORTANT_Call_Method_75_GM2.txt",delim="\t")
IMPORTANT_Call_Method_75_GM3<-read_delim("IMPORTANT_Call_Method_75_GM3.txt",delim="\t")
IMPORTANT_Call_Method_75_GM4<-read_delim("IMPORTANT_Call_Method_75_GM4.txt",delim="\t")
GenoFile1<-read_delim("Call_Method_75_GD1.txt",delim="\t")
GenoFile2<-read_delim("Call_Method_75_GD2.txt",delim="\t")
GenoFile3<-read_delim("Call_Method_75_GD3.txt",delim="\t")
GenoFile4<-read_delim("Call_Method_75_GD4.txt",delim="\t")
GenicListTAIR10_NoAA<-read_delim("GenicListTAIR10_NoAA.txt",delim="\t")

random_geno_list<-c(NA)

n=1
m=1
while (length(random_geno_list) < 51){

  
assign(sprintf('randomSubset_%s',n),rbind(GenicListTAIR10_NoAA[GenicListTAIR10_NoAA$chromosome_name==1,][sample(1:nrow(GenicListTAIR10_NoAA[GenicListTAIR10_NoAA$chromosome_name==1,]),77,replace=F),],GenicListTAIR10_NoAA[GenicListTAIR10_NoAA$chromosome_name==2,][sample(1:nrow(GenicListTAIR10_NoAA[GenicListTAIR10_NoAA$chromosome_name==2,]),49,replace=F),],GenicListTAIR10_NoAA[GenicListTAIR10_NoAA$chromosome_name==3,][sample(1:nrow(GenicListTAIR10_NoAA[GenicListTAIR10_NoAA$chromosome_name==3,]),73,replace=F),],GenicListTAIR10_NoAA[GenicListTAIR10_NoAA$chromosome_name==4,][sample(1:nrow(GenicListTAIR10_NoAA[GenicListTAIR10_NoAA$chromosome_name==4,]),57,replace=F),],GenicListTAIR10_NoAA[GenicListTAIR10_NoAA$chromosome_name==5,][sample(1:nrow(GenicListTAIR10_NoAA[GenicListTAIR10_NoAA$chromosome_name==5,]),79,replace=F),]))


assign(sprintf('randomSubset_%s_%s',n,'GM1'), foreach(i=1:nrow(get(sprintf('randomSubset_%s',n))), .combine='rbind') %do% subset(IMPORTANT_Call_Method_75_GM1, Chromosome == get(sprintf('randomSubset_%s',n))$chromosome_name[i] & Positions >= get(sprintf('randomSubset_%s',n))$start_position[i]-5000 & Positions <= get(sprintf('randomSubset_%s',n))$end_position[i]+5000))

assign(sprintf('randomSubset_%s_%s',n,'GM2'), foreach(i=1:nrow(get(sprintf('randomSubset_%s',n))), .combine='rbind') %do% subset(IMPORTANT_Call_Method_75_GM2, Chromsome == get(sprintf('randomSubset_%s',n))$chromosome_name[i] & Position >= get(sprintf('randomSubset_%s',n))$start_position[i]-5000 & Position <= get(sprintf('randomSubset_%s',n))$end_position[i]+5000))

assign(sprintf('randomSubset_%s_%s',n,'GM3'), foreach(i=1:nrow(get(sprintf('randomSubset_%s',n))), .combine='rbind') %do% subset(IMPORTANT_Call_Method_75_GM3, Chromsome == get(sprintf('randomSubset_%s',n))$chromosome_name[i] & Position >= get(sprintf('randomSubset_%s',n))$start_position[i]-5000 & Position <= get(sprintf('randomSubset_%s',n))$end_position[i]+5000))

assign(sprintf('randomSubset_%s_%s',n,'GM4'), foreach(i=1:nrow(get(sprintf('randomSubset_%s',n))), .combine='rbind') %do% subset(IMPORTANT_Call_Method_75_GM4, Chromsome == get(sprintf('randomSubset_%s',n))$chromosome_name[i] & Position >= get(sprintf('randomSubset_%s',n))$start_position[i]-5000 & Position <= get(sprintf('randomSubset_%s',n))$end_position[i]+5000))


assign(sprintf('randomSubset_%s_%s',n,'GM1'),unique(get(sprintf('randomSubset_%s_%s',n,'GM1'))))
assign(sprintf('randomSubset_%s_%s',n,'GM2'),unique(get(sprintf('randomSubset_%s_%s',n,'GM2'))))
assign(sprintf('randomSubset_%s_%s',n,'GM3'),unique(get(sprintf('randomSubset_%s_%s',n,'GM3'))))
assign(sprintf('randomSubset_%s_%s',n,'GM4'),unique(get(sprintf('randomSubset_%s_%s',n,'GM4'))))
if(
  nrow(get(sprintf('randomSubset_%s_%s',n,'GM1'))) >= nrow(AASubsetGM1) &&
  nrow(get(sprintf('randomSubset_%s_%s',n,'GM2'))) >= nrow(AASubsetGM2) &&
  nrow(get(sprintf('randomSubset_%s_%s',n,'GM3'))) >= nrow(AASubsetGM3) &&
  nrow(get(sprintf('randomSubset_%s_%s',n,'GM4'))) >= nrow(AASubsetGM4)
) {
  
  random_geno_list[m]<-paste0(sprintf('randomSubset_%s',n))

  
  assign(sprintf('randomSubset_%s_%s',n,'GM1'),get(sprintf('randomSubset_%s_%s',n,'GM1'))[sample(1:nrow(get(sprintf('randomSubset_%s_%s',n,'GM1'))),nrow(AASubsetGM1), replace=FALSE),])

assign(sprintf('randomSubset_%s_%s',n,'GM2'),get(sprintf('randomSubset_%s_%s',n,'GM2'))[sample(1:nrow(get(sprintf('randomSubset_%s_%s',n,'GM2'))),nrow(AASubsetGM2), replace=FALSE),])

assign(sprintf('randomSubset_%s_%s',n,'GM3'),get(sprintf('randomSubset_%s_%s',n,'GM3'))[sample(1:nrow(get(sprintf('randomSubset_%s_%s',n,'GM3'))),nrow(AASubsetGM3), replace=FALSE),])

assign(sprintf('randomSubset_%s_%s',n,'GM4'),get(sprintf('randomSubset_%s_%s',n,'GM4'))[sample(1:nrow(get(sprintf('randomSubset_%s_%s',n,'GM4'))),nrow(AASubsetGM4), replace=FALSE),])
  
assign(sprintf('randomSubset_%s_%s',n,'GD1'),  GenoFile1[,c(get(sprintf('randomSubset_%s_%s',n,'GM1'))$SNP_ID)])
assign(sprintf('randomSubset_%s_%s',n,'GD2'), GenoFile2[,c(get(sprintf('randomSubset_%s_%s',n,'GM2'))$Name)])
assign(sprintf('randomSubset_%s_%s',n,'GD3'), GenoFile3[,c(get(sprintf('randomSubset_%s_%s',n,'GM3'))$Name)])
assign(sprintf('randomSubset_%s_%s',n,'GD4'), GenoFile4[,c(get(sprintf('randomSubset_%s_%s',n,'GM4'))$Name)])

assign(sprintf('randomSubset_%s_%s',n,'GD1'),cbind(Ecotpe_ID=GenoFile1$Ecotype_ID,get(sprintf('randomSubset_%s_%s',n,'GD1'))))

assign(sprintf('randomSubset_%s_%s',n,'GD2'),cbind(Ecotpe_ID=GenoFile1$Ecotype_ID,get(sprintf('randomSubset_%s_%s',n,'GD2'))))
assign(sprintf('randomSubset_%s_%s',n,'GD3'),cbind(Ecotpe_ID=GenoFile1$Ecotype_ID,get(sprintf('randomSubset_%s_%s',n,'GD3'))))
assign(sprintf('randomSubset_%s_%s',n,'GD4'),cbind(Ecotpe_ID=GenoFile1$Ecotype_ID,get(sprintf('randomSubset_%s_%s',n,'GD4'))))

write_delim(get(sprintf('randomSubset_%s_%s',n,'GD1')),sprintf('randomSubset_%s_%s.txt',n,'GD1'),delim = "\t")
write_delim(get(sprintf('randomSubset_%s_%s',n,'GD2')),sprintf('randomSubset_%s_%s.txt',n,'GD2'),delim = "\t")
write_delim(get(sprintf('randomSubset_%s_%s',n,'GD3')),sprintf('randomSubset_%s_%s.txt',n,'GD3'),delim = "\t")
write_delim(get(sprintf('randomSubset_%s_%s',n,'GD4')),sprintf('randomSubset_%s_%s.txt',n,'GD4'),delim = "\t")

write_delim(get(sprintf('randomSubset_%s_%s',n,'GM1')),sprintf('randomSubset_%s_%s.txt',n,'GM1'),delim = "\t")
write_delim(get(sprintf('randomSubset_%s_%s',n,'GM2')),sprintf('randomSubset_%s_%s.txt',n,'GM2'),delim = "\t")
write_delim(get(sprintf('randomSubset_%s_%s',n,'GM3')),sprintf('randomSubset_%s_%s.txt',n,'GM3'),delim = "\t")
write_delim(get(sprintf('randomSubset_%s_%s',n,'GM4')),sprintf('randomSubset_%s_%s.txt',n,'GM4'),delim = "\t")

n<-n+1
m<-m+1

}
}
