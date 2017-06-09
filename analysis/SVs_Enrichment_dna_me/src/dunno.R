df = featureTilesDF[order(featureTilesDF[,foi[1]],decreasing = rev),foi]
source('H:/Dropbox/BRL/Methylation_Struct_Mut/MyScripts/MergeDFsFromBEDs.2.2.2016.R', echo=F)
z = MergeDFsFromBEDs()
#zz = as.matrix(z[,5:ncol(z)])
#class(zz[,5])
colnames(z)
my.foi= colnames(z)[c(5,13)]
getTileEnrichment(featureTilesDF = z,foi = my.foi,pct = 1,precision = 1000 )

#Are the sets the same for all 1% hypomethylated regions:
#identify all of the 1% regions and show when they turn on and off over time...
#Where a section is uniquely hypomethylated look for a subset of SVs 
#that fit there and then test to see if this si enriched vs ither genomic regions and time points
#start with...
length(intersect(x = z[order(z$mPGC.WK7.me),"row.id"][1:1500],y = z[order(z$sperm_CpGMe),"row.id"][1:1500]))
length(intersect(x = z[order(z$mPGC.WK11.me),"row.id"][1:1500],y = z[order(z$mPGC.WK13.me),"row.id"][1:1500]))

#SUZ12 has extraordiary enreachment near
#developmental germline 1% hypomethylatedsites
################################################################################

################################################################################
################################################################################
################################################################################
library(stats)
library(useful)
library(ggplot2)
library(gplots)
library(RColorBrewer)
TF_names = c('Cmyc.ENCODE','Ctcf.ENCODE','Foxa1.ENCODE','Jund.ENCODE',
             'Mef2a.ENCODE','Mef2c.ENCODE','Pou5f1.ENCODE','Prdm1.ENCODE',
             'PRDM14_ES_PEAKS','Rad21.ENCODE','Suz12.ENCODE','Yy1.ENCODE',
             'Znf143.ENCODE','Znf263.ENCODE','Znf274.ENCODE','HAVANA.genes')

SV_names = c('Rearrangement_Hum_Spec','Phase3_all_bkpts',
             'Phase3_mCNV','Phase3_DUP','Phase3_DEL','Phase3_INV',
             'CNV_Bipolar_Case','CNV_Bipolar_Control',
             'CNV_Autism_Case','CNV_Autism_Control',
             'CNV_DevDelay_Case','CNV_DevDelay_Control',
             'CNV_Schiz_Case','CNV_Schiz_Control',
             'CNV_WTCCC','CNV_270_HapMap',
             'CNV_450_HapMap','CNV_400_MGL')

Methylome_names = c('PGC_wk7', 'PGC_wk10','PGC_wk11','PGC_wk13','PGC_wk19', 'Sperm',
                    ' Heart_Wk5','Soma_wk10')

##############################################################################################
##############################################################################################
#Produce a heatmap for 1% lowestmethylation vs. Structrual Variants
#Sperm methylation deserts show strongest enrichment
#MakeHeatmap.hi.lo.Me.CNV = function(MethDf,CNVDf,pct=1,prec = 10000,rev =F){
    MethDf = z[,5:12]
    CNVDf = z[,13:29]
    pct = 1
    rev = F
    prec = 10000
    All.Enrichments1 = data.frame(matrix(data = NA,nrow =500,ncol = 3))
    k =1 #initialize counter
    df =cbind(MethDf,CNVDf)
    for(meth in 1:ncol(MethDf)){
        for(sv in 1:ncol(CNVDf)){
            my2features = c(meth,sv)
            tmp.enr = getTileEnrichment(featureTilesDF = df, foi = my2features,
                                        pct = pct, precision = prec, rev = rev)
            rownames(All.Enrichments1)[k] = tmp.enr[[1]]
            All.Enrichments1[k,] = c(tmp.enr[[1]],tmp.enr[[2]],tmp.enr[[3]])
            k = k+1
        }
    }
    colnames(All.Enrichments1) = names(tmp.enr)
    All.Enrichments1 = All.Enrichments1[complete.cases(All.Enrichments1),]
    #write.table(x = All.Enrichments1,file = "EnrichmentComparisons1.txt",
    #            col.names = T,row.names = F, sep = "\t",quote = F)
    dimsr = list(SV_names,Methylome_names)
    r = matrix((as.numeric(All.Enrichments1[,2])),nrow = length(SV_names), dimnames = dimsr)
    r[r==0] =NA
    r
    #brk.max = max(abs(log2(r)))
    brk.max = 3.2
    dev.new()
    #par(bg = "transparent")
    par(bg = "grey85")
    numBreaks = 256
    pairs.breaks = seq(from=-brk.max ,to=brk.max,length.out=257)
    heatmap.2(log(r,base = 2),na.rm = T,trace = 'none',margins = c(10,18),
              breaks = pairs.breaks,
              symbreaks = any(breaks<0,na.rm = T),
              density.info = 'none', key.xlab = 'Log2 Fold Enrichment', key.title = NA,
              cexCol = 2,cexRow = 1.5,
              col = bluered(numBreaks),Colv = F,Rowv = F,
              main = "Enrichment of Structural Variants in 100kb Tiles\n with HI 1% CpG Methylation\nin human male developmental germline" )
}
colnames(z)

MakeHeatmap.hi.lo.Me.CNV(MethDf = z[,5:12],CNVDf = z[,14:29],pct = 1,rev = F)
##############################################################################################
##############################################################################################
#KS-Test for distrubution of CNV in %methylated regions
source('H:/Dropbox/BRL/Methylation_Struct_Mut/MyScripts/getTileKS.R', echo=F)
#x = getTile.KS(z,c(6,13))
#x 

colnames(z)
All.KS = data.frame(matrix(data = NA,nrow =500,ncol = 3))
k =1 #initialize counter
for(meth in 5:12){
    for(sv in 13:30){
        my.foi = c(meth,sv)
        tmp.enr = getTile.KS(featureTilesDF = z,my.foi)
        #rownames(All.KS)[k] = paste(colnames(z)[my.foi[1]],colnames(z)[my.foi[1]],sep = " ")
        All.KS[k,] = c(tmp.enr[[1]],tmp.enr[[2]],tmp.enr[[3]])
        k = k+1
    }
}
colnames(All.KS) = c('D.Statistic','p.val','test-type')
All.KS = All.KS[complete.cases(All.KS),]
All.KS.sig = All.KS
p.val = 0.5
Bonf = nrow(All.KS)
All.KS.sig[as.numeric(All.KS.sig[,2])>(p.val/Bonf),1] = NA

#write.table(x = All.Enrichments1,file = "EnrichmentComparisons1.txt",
#            col.names = T,row.names = F, sep = "\t",quote = F)
dimsr = list(SV_names,Methylome_names)
#sigKS = All.KS[,2]>(.05/144)

r = matrix((as.numeric(All.KS.sig[,1])),nrow = length(SV_names), dimnames = dimsr)
p = matrix((as.numeric(All.KS.sig[,2])),nrow = length(SV_names), dimnames = dimsr)
r
p
dev.new()
par(bg = "transparent")
par(bg = "white")
numBreaks = 256
#pairs.breaks = seq(from=-bb,to=bb,length.out=257)
redramp = colorRampPalette(c('gray','red'))
heatmap.2(r,na.rm = T,trace = 'none',margins = c(10,18),
          #breaks = pairs.breaks,
          #symbreaks = any(breaks<0,na.rm = T),
          density.info = 'none', key.xlab = 'KS-Test D Statistic', key.title = NA,
          cexCol = 2,cexRow = 1.5, sepcolor = "black", colsep = 0:(ncol(r)+2),rowsep = 1:nrow(r),
          col = redramp(numBreaks),Colv = F,Rowv = F,
          main = "KS TEST D-Statistic for %CpG methylation difference\nin 100kb Tiles with Structural Vs. Rest of genome\n Bonferroni-Corrected p < .01" )




##############################################################################################
##############################################################################################
###TRANSCRIPTION FACTORS AND ANNOTATIONS Vs. DNAme  
#
pct = 1
prec =1000
Me_TF.Enrichments1 = data.frame(matrix(data = NA,nrow =500,ncol = 3))
k =1
for(meth in 5:10){
    for(sv in 29:ncol(z)){
        my.foi = c(meth,sv)
        tmp.enr = getTileEnrichment(featureTilesDF = z,my.foi,pct=pct,precision =  prec )
        rownames( Me_TF.Enrichments1 )[k] = tmp.enr[[1]]
        Me_TF.Enrichments1[k,] = c(tmp.enr[[1]],tmp.enr[[2]],tmp.enr[[3]])
        k = k+1
    }
}
colnames( Me_TF.Enrichments1 ) = names(tmp.enr)
Me_TF.Enrichments1 = Me_TF.Enrichments1[complete.cases(Me_TF.Enrichments1),]
#write.table(x = Me_TF.Enrichments1,file = "TF_vs_DNAMe_EnrichmentComparisons1.txt",
#            col.names = T,row.names = F, sep = "\t",quote = F)
dimsr = list(TF_names,Methylome_names)
r = matrix((as.numeric(Me_TF.Enrichments1[,2])),nrow = length(TF_names), dimnames = dimsr)
dev.new()

heatmap.2(log(r,base = 2),na.rm = T,trace = 'none',margins = c(10,12),
          col = bluered(40),Colv = F,Rowv = T,main = "Enrichment of TF binding (ENCODE)\n in 100kb Tiles with lowest 1% CpG Methylation\nin human male developmental germline" )
################################################################################
################################################################################
###Trying to reset color scale
###TRANSCRIPTION FACTORS AND ANNOTATIONS Vs. DNAme###  
#
pct = 5
prec =1000
Me_TF.Enrichments5 = data.frame(matrix(data = NA,nrow =500,ncol = 3))
k =1
for(meth in 5:10){
    for(sv in 29:ncol(z)){
        my.foi = c(meth,sv)
        tmp.enr = getTileEnrichment(featureTilesDF = z,my.foi,pct=pct,precision =  prec )
        rownames( Me_TF.Enrichments1 )[k] = tmp.enr[[1]]
        Me_TF.Enrichments5[k,] = c(tmp.enr[[1]],tmp.enr[[2]],tmp.enr[[3]])
        k = k+1
    }
}
colnames( Me_TF.Enrichments5 ) = names(tmp.enr)
Me_TF.Enrichments5 = Me_TF.Enrichments5[complete.cases(Me_TF.Enrichments5),]
#write.table(x = Me_TF.Enrichments1,file = "TF_vs_DNAMe_EnrichmentComparisons1.txt",
#            col.names = T,row.names = F, sep = "\t",quote = F)
dimsr = list(TF_names,Methylome_names)
r = matrix((as.numeric(Me_TF.Enrichments5[,2])),nrow = length(TF_names), dimnames = dimsr)
dev.new()
heatmap.2(log(r,base = 2),na.rm = T,trace = 'none',margins = c(10,12),breaks = c(seq(-.5,5,length=41)),
          col = bluered(40),Colv = F,Rowv = T,main = "Enrichment of TF binding (ENCODE)\n in 100kb Tiles with lowest 5% CpG Methylation\nin human male developmental germline" )

