
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

RARE_EAS_CNV.SV.bed <- read.csv(
  "C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/ROI_Data/BED_Format/Output/tmpdist/RARE_EAS_CNV.SV.bed.dist",
  sep="", stringsAsFactors=FALSE
)
RAR_CNV_List = list.files(pattern = "RAR.*EAS.*CNV.SV.bed.sampled..dist")
plot(ecdf(RARE_EAS_CNV.SV.bed[,1]), xlim = c(0,2000000))
for(file in RAR_CNV_List ){
  a = read.table(file)[,1]
  lines(ecdf(read.table(file)[,1]),col = 'purple')
}
RAR_CNV_Control_List = list.files(pattern = "RAR.*EAS.*CNV.SV.bed.sampled.*toshuff")
for(file in RAR_CNV_Control_List ){
  b = read.table(file)[,1]
  lines(ecdf(read.table(file)[,1]),col = 'black')
}
