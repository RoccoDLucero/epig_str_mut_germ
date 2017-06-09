
ks.test(x = RARE_AFR_CNV.SV.bed.sampled$V1,y = RARE_EUR_CNV.SV.bed.sampled$V1)
ks.test(x = RARE_AFR_CNV.SV.bed.sampled$V1,y = RARE_EAS_CNV.SV.bed.sampled$V1)
ks.test(x = RARE_AFR_CNV.SV.bed.sampled$V1,y = RARE_AMR_CNV.SV.bed.sampled$V1)
ks.test(x = RARE_AFR_CNV.SV.bed.sampled$V1,y = RARE_SAS_CNV.SV.bed.sampled$V1)

ks.test(x = RARE_AMR_CNV.SV.bed.sampled$V1,y = RARE_EAS_CNV.SV.bed.sampled$V1)


RARE_AFR_CNV.SV.bed <- read.csv(
  "C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/ROI_Data/BED_Format/Output/tmpdist/RARE_AFR_CNV.SV.bed.dist",
  sep="", stringsAsFactors=FALSE
  )
RAR_CNV_List = list.files(pattern = "RAR.*AFR.*CNV.SV.bed.sampled*")
plot(ecdf(RARE_AFR_CNV.SV.bed[,1]))
for(file in RAR_CNV_List ){
  lines(ecdf(read.table(file)[,1]),col = 'purple')
}
RAR_CNV_Control_List = list.files(pattern = "RAR.*AFR.*CNV.SV.bed.sampled.*toshuff")
for(file in RAR_CNV_Control_List ){
  read.table(file)[,1]
  lines(ecdf(read.table(file)[,1]),col = 'black')
}

COM_CNV_List = list.files(pattern = "COM.*CNV.SV.bed.sampled*")
lines(ecdf(RARE_AFR_CNV.SV.bed[,1]))
for(file in COM_CNV_List ){
  lines(ecdf(read.table(file)[,1]),col = 'aquamarine')
}

DUP_List = list.files(pattern = "*DUP.SV.bed.sampled*")
lines(ecdf(RARE_AFR_DUP.SV.bed.shuffleChIP),col ='hotpink2')
for(file in DUP_List ){
  lines(ecdf(read.table(file)[,1]),col = 'red')
}


DEL_List = list.files(pattern = "*DEL.SV.bed.sampled*")
lines(ecdf(RARE_AFR_DEL.SV.bed.shuffleChIP[,1]),col ='goldenrod')
for(file in DEL_List ){
  #later add code to plot line color by population 
  lines(ecdf(read.table(file)[,1]),col = 'chartreuse')  ##rgb(runif(5),runif(5),runif(5)))
}

