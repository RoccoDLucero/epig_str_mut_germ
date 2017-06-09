#MethylKit for PGC Methylation level in 100kb and 1kb windows
setwd("H:/Dropbox/BRL/Methylation_Struct_Mut/Processed_Data")
dev.new()
library(methylKit)
library(graphics)
file.list = list(
    "./GSE63818_Heart_5W_embryo1_methylation_calling.CpG.bed.sorted",
    "./GSE63818_Soma_7W_embryo1_M_methylation_calling.CpG.bed.sorted",
    "./GSE63818_Soma_10W_embryo1_M_methylation_calling.CpG.bed.sorted",
    "./GSE63818_Soma_11W_embryo1_M_methylation_calling.CpG.bed.sorted",
    "./GSE63818_Soma_17W_embryo1_F_methylation_calling.CpG.bed.sorted",
    "./GSE63818_Soma_19W_embryo1_M_methylation_calling.CpG.bed.sorted"
)
sample.list = list("Emb1_5w_Heart","Emb1_M_7w_Soma","Emb1_M_10w_Soma","Emb1_M_11w_Soma",
                   "Emb1_F_17w_soma","Emb1_M_19w_Soma")
#sample.list = list("Emb1_5w_Heart","Emb1_10w_Soma")
treatment = c(1,1,1,1,1,1)
#treatment = c(1,1)
headerMap = list(fraction=T, chr.col=1,start.col=2,end.col=2,
                 coverage.col=5,strand.col=4,freqC.col=8)


myobj = read(location = file.list, sample.id = sample.list, assembly = "hg19",
             context = "CpG",treatment = treatment,header = F,
             pipeline= headerMap)

meth = unite(object = myobj, destrand = F,min.per.group = 1L)
#methAllSmpCvr = unite(object = myobj, destrand = F)

methTiles100kb = tileMethylCounts(object = myobj, win.size =100000,step.size = 100000 )
allTiles100kb = unite(object = methTiles100kb, destrand = F)
dev.new()
#getCoverageStats(object = myobj[[1]],plot = F)
getMethylationStats(object = myobj[[1]],plot = T)
#getCorrelation(allTiles100kb,plot = T)
perc.Meth100kb =percMethylation(allTiles100kb)
perc.Meth100kb = cbind(getData(allTiles100kb)[,1:3],perc.Meth100kb)
write.table(x = perc.Meth100kb,
            file = "./GSE63818_Somatic_100kb_hg19_methylation.txt",
            row.names = F,col.names = F,quote = F,sep = "\t")


methTiles1kb = tileMethylCounts(object = myobj, win.size =1000,step.size = 1000 )
allTiles1kb = unite(object = methTiles1kb, destrand = F)
dev.new()
getCorrelation(allTiles1kb,plot = T)
perc.Meth1kb =percMethylation(allTiles1kb)
perc.Meth1kb = cbind(getData(allTiles1kb)[,1:3],perc.Meth1kb)
write.table(x = perc.Meth1kb,
            file = "./GSE63818_Somatic_1kb_hg19_methylation.txt",
            row.names = F,col.names = T,quote = F,sep = "\t")


################################################################################
################################################################################
#This is not a robust block of code
#fiveWkGetDeltaMeth should compute change in methylation levels for each 
# of the four transitions and the total amount of methylation change that
# can be attributed to each tile/or base
#fiveWkGetDeltaMeth = function(Matrx){
#    for (i in 4:7){
#        Matrx[,(ncol(Matrx)+1)] =
#            as.numeric(Matrx[,(i+1)]) - as.numeric(Matrx[,i])
#    }
#    Matrx$totalDiff = rowSums(abs(Matrx[,9:12]))
#    Matrx
#}
#xx = fiveWkGetDeltaMeth(perc.Meth100kb)
#yy = fiveWkGetDeltaMeth(perc.Meth1kb)
#write.table(x = xx,
#            file = "./GSE63818_Soma_7-19W_embryo1_M_100kb_hg19_meth_plus_del.txt",
#            row.names = F,col.names = T,quote = F,sep = "\t")
#write.table(x = yy,
#            file = "./GSE63818_Soma_7-19W_embryo1_M_1kb_hg19_meth_plus_del.txt",
#            row.names = F,col.names = T,quote = F,sep = "\t")
#rm(xx,yy)

##Performed IntersectBeds at command line for all 1KgPhase3 breakpoints
##Vs Tiled windows.

#GSE63818_PGC_7.19W_embryo1_M_100kb_hg19_meth_V_AllPhase3 <- read.delim("H:/Dropbox/BRL/Methylation_Struct_Mut/Processed_Data/GSE63818_PGC_7-19W_embryo1_M_100kb_hg19_meth_V_AllPhase3.txt", header=FALSE)
#head(GSE63818_PGC_7.19W_embryo1_M_100kb_hg19_meth_V_AllPhase3)

#summary(GSE63818_PGC_7.19W_embryo1_M_100kb_hg19_meth_V_AllPhase3)
#ff =GSE63818_PGC_7.19W_embryo1_M_100kb_hg19_meth_V_AllPhase3[which(GSE63818_PGC_7.19W_embryo1_M_100kb_hg19_meth_V_AllPhase3$V9>=2),]
#dim(ff)
#summary(ff)
#ee = ff[which(ff$V9>8),]
#summary(ee)
#dim(ee)
#ee
#tt = ff[which(ff$V9<=8),]
#summary(tt)


GSE30340_human_sperm_CpG_coverage_hg18 = read.delim("H:/Dropbox/BRL/Methylation_Struct_Mut/Raw_Data/GSE30340_human_sperm_CpG_coverage_hg18.wig", header=FALSE,row.names = F)
GSE30340_human_sperm_CpG_methylation_hg18 = read.delim("H:/Dropbox/BRL/Methylation_Struct_Mut/Raw_Data/GSE30340_human_sperm_CpG_methylation_hg18.wig", header=FALSE,row.names = F)
GSE30340_human_sperm_CpG_coverage_hg18_5cov = GSE30340_human_sperm_CpG_coverage_hg18[which(GSE30340_human_sperm_CpG_coverage_hg18$V4>=5),]
MolaroCpGplusCoverage5_hg18 = merge(GSE30340_human_sperm_CpG_methylation_hg18,GSE30340_human_sperm_CpG_coverage_hg18_5cov,by.x = c(1:3),by.y = c(1:3))
write.table(x = MolaroCpGplusCoverage5_hg18,
            file = "./MolaroCpGplusCoverage5_hg18.txt",
            row.names = F,col.names = F,quote = F,sep = "\t")
########################################################################################
############################################################################################################
MOLARO SPERM DATA and HARRIS CNV DATA
############################################################################################################

library(plyr)
Confounders_Data = read.delim("H:/Dropbox/BRL/Methylation_Struct_Mut/Raw_Data/Harris_Confounders_Data/Table_S1/Confounders_Data.txt")
head(arrange(df = Confounders_Data, Chromosome,Start)[,1:6],400)
mixedorder(Confounders_Data$Chromosome)
#file.list = list("H:/Dropbox/BRL/Methylation_Struct_Mut/Processed_Data/GSE30340_human_sperm_CpG_met_cov5_hg19_molaro.wig")
file.list = list("H:/Dropbox/BRL/Methylation_Struct_Mut/Processed_Data/tmp.wig",
                 "H:/Dropbox/BRL/Methylation_Struct_Mut/Processed_Data/tmp.wig")

#sample.list = list("Molaro_Sperm_CpgMe")
sample.list = list("Molaro_Sperm_CpgMe","Molaro_Sperm_CpgMe")
headerMap = list(fraction=T, chr.col=1,start.col=2,end.col=3,freqC.col=4, coverage.col=5,strand.col = 6)
#treatment = c(0,0)
treatment = c(0,0)


#myobjMolaro = read(location = file.list[[1]], sample.id = sample.list[[1]], assembly = "hg19",
#             context = "CpG",treatment = treatment, header = F,
#             pipeline= headerMap)
#myobjMolaro1 = read.bed(location = file.list[[1]],remove.unsual = T)

myobjMolaro2 = read(location = file.list, sample.id = sample.list, assembly = "hg19",
                    context = "CpG",treatment = treatment, header = F,
                    pipeline= headerMap)
#meth = unite(object = myobjMolaro2, destrand = F,min.per.group = 1L)

methTiles1kb = tileMethylCounts(object = myobjMolaro2, win.size =1000,step.size = 1000 )
allTiles1kb = unite(object = methTiles1kb)
perc.Meth1kb =percMethylation(allTiles1kb)
perc.Meth1kb = cbind(getData(allTiles1kb)[,1:3],perc.Meth1kb[,1])

methTiles100kb = tileMethylCounts(object = myobjMolaro2, win.size =100000,step.size = 100000 )
allTiles100kb = unite(object = methTiles100kb, destrand = F)
perc.Meth100kb =percMethylation(allTiles100kb)
perc.Meth100kb = cbind(getData(allTiles100kb)[,1:3],perc.Meth100kb[,1])

methTiles100bp = tileMethylCounts(object = myobjMolaro2, win.size =100,step.size = 100 )
allTiles100bp = unite(object = methTiles100bp, destrand = F)
perc.Meth100bp =percMethylation(allTiles100bp)
perc.Meth100bp = cbind(getData(allTiles100bp)[,1:3],perc.Meth100bp[,1])


write.table(x = perc.Meth1kb,
            file = "./GSE30340_human_sperm_CpG_1kbTile_hg19_molaro.wig",
            row.names = F,col.names = F,quote = F,sep = "\t")

write.table(x = perc.Meth100kb,
            file = "./GSE30340_human_sperm_CpG_100kbTile_hg19_molaro.wig",
            row.names = F,col.names = F,quote = F,sep = "\t")

write.table(x = perc.Meth100bp,
            file = "./GSE30340_human_sperm_CpG_100bpTile_hg19_molaro.wig",
            row.names = F,col.names = F,quote = F,sep = "\t")


setwd("H:/Dropbox/BRL/Methylation_Struct_Mut/Processed_Data")
Phase3_all_bkpts <- read.delim("./1KG_phase3_all_bkpts.v5_0.1Mb.hg19.bed", header=FALSE)

Phase3SVs = `1KG_phase3_all_bkpts.v5_0.1Mb.hg19`
rm(`1KG_phase3_all_bkpts.v5_0.1Mb.hg19`)
Phase3SVs$row.id = paste(Phase3SVs$V1,Phase3SVs$V2,Phase3SVs$V3,sep = ".")

Phase3SVs = Phase3SVs[,c(5,1:4)]

Rearrangement_Hum_Spec <- read.delim("./Rearrangement_Human_Specific.0.1Mb.hg19.bed", header=FALSE)
Rearrangement_Human_Specific.0.1Mb.hg19$row.id = paste(Rearrangement_Human_Specific.0.1Mb.hg19$V1,
                                                       Rearrangement_Human_Specific.0.1Mb.hg19$V2,
                                                       Rearrangement_Human_Specific.0.1Mb.hg19$V3,sep = ".")
Rearrangement_Human_Specific.0.1Mb.hg19 =Rearrangement_Human_Specific.0.1Mb.hg19[,c(5,1:4)] 

aa = merge(x = Phase3SVs, y = Rearrangement_Human_Specific.0.1Mb.hg19)
aa = aa[,c(1,5,6)]

PGC_7.19W_embryo1_methylation <- read.delim("./GSE63818_PGC_7-19W_embryo1_M_100kb_hg19_methylation.txt", header=FALSE)

GSE63818_PGC_7.19W_embryo1_M_100kb_hg19_methylation$row.id = paste(
    GSE63818_PGC_7.19W_embryo1_M_100kb_hg19_methylation$V1,
    GSE63818_PGC_7.19W_embryo1_M_100kb_hg19_methylation$V2,
    GSE63818_PGC_7.19W_embryo1_M_100kb_hg19_methylation$V3,sep = '.'
)
GSE63818_PGC_7.19W_embryo1_M_100kb_hg19_methylation = GSE63818_PGC_7.19W_embryo1_M_100kb_hg19_methylation[,c(9,4:8)]
colnames(GSE63818_PGC_7.19W_embryo1_M_100kb_hg19_methylation)[2:6] =
    c('mPGC.WK7.me','mPGC.WK10.me','mPGC.WK11.me','mPGC.WK13.me','mPGC.WK19.me')

hg19_0.1Mb_windowCoords.sorted <- read.delim("H:/Dropbox/BRL/Methylation_Struct_Mut/Processed_Data/hg19_0.1Mb_windowCoords.sorted.bed", header=FALSE)
colnames(hg19_0.1Mb_windowCoords.sorted) = c("chr","start","stop")
hg19_0.1Mb_windowCoords.sorted$row.id = paste(hg19_0.1Mb_windowCoords.sorted$chr,
                                              hg19_0.1Mb_windowCoords.sorted$start,
                                              hg19_0.1Mb_windowCoords.sorted$stop,sep = "."
)
bb = merge.data.frame(x = hg19_0.1Mb_windowCoords.sorted,
                      y = GSE63818_PGC_7.19W_embryo1_M_100kb_hg19_methylation,
                      all.x = T,by = "row.id"
)

cc = merge.data.frame(x = bb,y = aa,all.x = T)

sperm_CpG_100kbTile_hg19_molaro <- read.delim("H:/Dropbox/BRL/Methylation_Struct_Mut/Processed_Data/GSE30340_human_sperm_CpG_100kbTile_hg19_molaro.sorted.wig", header=FALSE, stringsAsFactors=FALSE)
sperm_CpG_100kbTile_hg19_molaro$row.id = paste(sperm_CpG_100kbTile_hg19_molaro$V1,
                                               sperm_CpG_100kbTile_hg19_molaro$V2,
                                               sperm_CpG_100kbTile_hg19_molaro$V3,sep = ".")

sperm_CpG_100kbTile_hg19_molaro = sperm_CpG_100kbTile_hg19_molaro[,c(4,5)]
colnames(sperm_CpG_100kbTile_hg19_molaro)[1] = c("hum.sperm.Me")
dd = merge.data.frame(x = cc,y = sperm_CpG_100kbTile_hg19_molaro,all.x = T)
head(dd)
dd = dd[,c(1:9,12,10,11)]

write.table(x = dd, file = "./hg19_PGC-Sperm-humRe-Phs3AllSv.txt",row.names = F)

##Import the rest of the Structural Variant Data Sets from  Li/Harris and 1KG
PGC_7.19W_embryo1_methylation <- read.delim("./GSE63818_PGC_7-19W_embryo1_M_100kb_hg19_methylation.txt", header=FALSE)
Rearrangement_Hum_Spec <- read.delim("./Rearrangement_Human_Specific.0.1Mb.hg19.bed", header=FALSE)
Phase3_all_bkpts <- read.delim("./1KG_phase3_all_bkpts.v5_0.1Mb.hg19.bed", header=FALSE)

CNV_Bipolar_Case <- read.delim("./singletonCNV_Bipolar_Case.0.1Mb.hg19.bed", header=FALSE)
CNV_Bipolar_Control <- read.delim("./singletonCNV_Bipolar_Control.0.1Mb.hg19.bed", header=FALSE)
CNV_Autism_Case <- read.delim("./rareCNV_Autism_Case.0.1Mb.hg19.bed", header=FALSE)
CNV_Autism_Control <- read.delim("./rareCNV_Autism_Control.0.1Mb.hg19.bed", header=FALSE)
CNV_DevDelay_Case <- read.delim("./rareCNV_DevDelay_Case.0.1Mb.hg19.bed", header=FALSE)
CNV_DevDelay_Control <- read.delim("./rareCNV_DevDelay_Control.0.1Mb.hg19.bed", header=FALSE)
CNV_Schiz_Case  <- read.delim("./rareCNV_Schizophrenia_Case.0.1Mb.hg19.bed", header=FALSE)
CNV_Schiz_Control <- read.delim("./rareCNV_Schizophrenia_Control.0.1Mb.hg19.bed", header=FALSE)
CNV_270_HapMap <- read.delim("./CNV_270HapMap.0.1Mb.hg19.bed", header=FALSE)
CNV_450_HapMap <- read.delim("./CNV_450HapMap.0.1Mb.hg19.bed", header=FALSE)
CNV_400_MGL <- read.delim("./CNV_400MGL.0.1Mb.hg19.bed", header=FALSE)
CNV_WTCCC <- read.delim("./CNV_WTCCC.0.1Mb.hg19.bed", header=FALSE)

HarrisCNV = list(PGC_7.19W_embryo1_methylation, Rearrangement_Hum_Spec, Phase3_all_bkpts
                 CNV_Bipolar_Case,CNV_Bipolar_Control,CNV_Autism_Case,CNV_Autism_Control,
                 CNV_DevDelay_Case,CNV_DevDelay_Control,CNV_Schiz_Case,CNV_Schiz_Control,
                 CNV_WTCCC,CNV_270_HapMap,CNV_450_HapMap,CNV_400_MGL)
for (df in HarrisCNV){
    df$row.id = paste(df$V1,df$V2,df$V3,sep = ".")
    df = df[,c(5,1:4)]
    dd = merge(x = dd, y = df)
}








#Get the entire table ordered by increasing methylation level in sperm.
ee = dd[complete.cases(dd),]
ee = ee[order(ee$hum.sperm.Me),]
#ASK some questions:
# 1. Given that we have moved to hg19 coordinates is the lowest 5% in sperm
#    still enriched for Human Specific Rearrangements
#Write the lowest 5% and the other 95% of tiles
#as a table and do intersect beds on each with the sets of HR and all SVS
0.05*length(ee$hum.sperm.Me[complete.cases(ee)])
0.01*length(ee$hum.sperm.Me[complete.cases(ee)]) 
write.table(x = ee[1:1430,c(2:4,10)], file = "./hg19_Sperm_Lo5.txt",
            row.names = F,quote = F,sep = "\t",col.names = F)
write.table(x = ee[1431:length(ee$hum.sperm.Me),c(2:4,10)],
            file = "./hg19_Sperm_Hi95.txt",row.names = F,quote = F,
            sep = "\t",col.names = F)
#
#
ee = ee[order(ee$hum.sperm.Me),]
pct = 5
totalTiles = length(ee$hum.sperm.Me[complete.cases(ee)])
testTiles = round(pct/100*totalTiles,digits = 0)
testEnrich = sum(ee$HumanSpRe[1:testTiles])/testTiles 
otherEnrich = sum(ee$HumanSpRe[(testTiles+1):totalTiles])/(totalTiles - testTiles) 
Sperm5pctHumanSpecifcRe = testEnrich/otherEnrich

pct = 1
totalTiles = length(ee$hum.sperm.Me[complete.cases(ee)])
testTiles = round(pct/100*totalTiles,digits = 0)
testEnrich = sum(ee$HumanSpRe[1:testTiles])/testTiles 
otherEnrich = sum(ee$HumanSpRe[(testTiles+1):totalTiles])/(totalTiles - testTiles) 
Sperm1pctHumanSpecifcRe = testEnrich/otherEnrich

################################################################################
################################################################################
#
#
#
#
#
#
#
#
gg = ee[order(ee$Phs3AllSVs,decreasing = T),]
head(gg,400)
dev.new()
plot(ecdf(gg$hum.sperm.Me))
lines(ecdf(gg$mPGC.WK7))
lines(ecdf(gg$mPGC.WK10))
lines(ecdf(gg$mPGC.WK11))
lines(ecdf(gg$mPGC.WK13))
lines(ecdf(gg$mPGC.WK19))

