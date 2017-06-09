source('H:/Dropbox/BRL/Methylation_Struct_Mut/MyScripts/MergeDFsFromBEDs.2.2.2016.R', echo=F)
source('H:/Dropbox/BRL/Methylation_Struct_Mut/MyScripts/MergeDFsFromBEDs2_2.18.2016.R', echo=F)
z = MergeDFsFromBEDs2()
colnames(z)
zz = z
z = zz

##for germline deltas d and dd
x = z[,12:16]; y = z[,11:15]; w = (x-y)
##for Somatic deltas
#q = z[,6:10];r = z[,5:9]; w = q-r;

ww = (w[,2:5]-w[1:4])

colnames(w) =  c("d7.10","d10.11","d11.13","d13.19","d19.sprm")
colnames(ww) = c("dd7.11",'dd10.13','dd11.19','dd13.sprm')
#w = cbind(z[,1:16],w,z[,17:ncol(z)])
ww = cbind(z[,1:16],w,ww,z[,17:ncol(z)])
z =ww
#zz = as.matrix(z[,5:ncol(z)])
#class(zz[,5])
colnames(z)
#my.foi= colnames(z)[c(21,22)] #17-21
#jnk = getTileEnrichment2(featureTilesDF = z,foi = my.foi,pct = 1,precision = 1000 )
#jnk[[2]]
#Are the sets the same for all 1% hypomethylated regions:
#identify all of the 1% regions and show when they turn on and off over time...
#Where a section is uniquely hypomethylated look for a subset of SVs 
#that fit there and then test to see if this si enriched vs ither genomic regions and time points
#start with...
#length(intersect(x = z[order(z$mPGC.WK7.me),"row.id"][1:1500],y = z[order(z$sperm_CpGMe),"row.id"][1:1500]))
#length(intersect(x = z[order(z$mPGC.WK11.me),"row.id"][1:1500],y = z[order(z$mPGC.WK13.me),"row.id"][1:1500]))

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
TF_names = colnames(z)[44:146] 

SV_names = c('Rearrangement_Hum_Spec','Phase3_all_bkpts',
             'Phase3_mCNV','Phase3_DUP','Phase3_DEL','Phase3_INV',
             'CNV_Bipolar_Case','CNV_Bipolar_Control',
             'CNV_Autism_Case','CNV_Autism_Control',
             'CNV_DevDelay_Case','CNV_DevDelay_Control',
             'CNV_Schiz_Case','CNV_Schiz_Control',
             'CNV_WTCCC','CNV_270_HapMap',
             'CNV_450_HapMap','CNV_400_MGL')

Methylome_names = c('Heart_Wk5','Soma_wk10','Soma_wk11','Soma_wk13','Soma_wk17','Soma_wk19',
                    'PGC_wk7', 'PGC_wk10','PGC_wk11','PGC_wk13','PGC_wk19', 'Sperm',
                    "d7.10","d10.11","d11.13","d13.19","d19.sprm",
                    "dd7.11",'dd10.13','dd11.19','dd13.sprm')

##############################################################################################
##############################################################################################
##Heatmaps for methylation-level-tiles enrichment with CNV
source('H:/Dropbox/BRL/Methylation_Struct_Mut/MyScripts/getTileEnrichment.R', echo=TRUE)
source('H:/Dropbox/BRL/Methylation_Struct_Mut/MyScripts/getTileEnrichment2.R', echo=TRUE)
source('H:/Dropbox/BRL/Methylation_Struct_Mut/MyScripts/featureTilesHeatmap.R', echo=TRUE)
source('H:/Dropbox/BRL/Methylation_Struct_Mut/MyScripts/featureTilesHeatmap2.R', echo=TRUE)

featureTilesHeatmap(featureTilesDF = z,firstSet = 5:25 ,secondSet = 26:43,pct = 1,rev = F)
featureTilesHeatmap(featureTilesDF = z,firstSet = 5:25 ,secondSet = 26:43,pct = 1,rev = T)
mya = featureTilesHeatmap2(featureTilesDF = z,firstSet = 5:25 ,secondSet = 26:43,pct = 1,rev = F)
myb = featureTilesHeatmap2(featureTilesDF = z,firstSet = 5:25 ,secondSet = 26:43,pct = 1,rev = T)


mmyyaa = mya[seq(from = 1, to =378,by=18)]
mmyybb = myb[seq(from = 1, to =378,by=18)]
SharedTiles = c()
for(set1 in 1:length(mmyyaa)){
    for(set2 in 1:length(mmyyaa)){
        SharedTiles = c(SharedTiles,length(intersect(mmyyaa[[set1]],mmyyaa[[set2]]))                )
    }
}
SharedTilesMatrix = matrix(data = SharedTiles,
                           nrow = length(mmyyaa),
                           dimnames = list(Methylome_names,Methylome_names))
SharedTilesMatrix = 100*SharedTilesMatrix/(max(SharedTilesMatrix))
mycols = colorRampPalette(c("white","red"))
dev.new()
heatmap.2(  SharedTilesMatrix,
            trace = 'none',
            Colv = F,Rowv = F,
            col = mycols(50),
            margins = c(8,8),
            main = paste("Percent intersect for 100kb tiles with lowest 1% methylation,\n",
                         "most negative first and second deltas in germline development")
            
        )

SharedTiles2 = c()
for(set1 in 1:length(mmyyaa)){
    for(set2 in (1:length(mmyybb))){
        SharedTiles2 = c(SharedTiles2,length(intersect(mmyyaa[[set1]],mmyybb[[set2]]))                )
    }
}

SharedTilesMatrix2 = matrix(data = SharedTiles2,
                           nrow = length(mmyyaa),
                           byrow = F,
                           dimnames = list(Methylome_names,(Methylome_names)))

SharedTilesMatrix2 = 100*SharedTilesMatrix2/(max(SharedTilesMatrix2))
mycols = colorRampPalette(c("white","purple"))
dev.new()
heatmap.2(  SharedTilesMatrix2,
            trace = 'none',
            Colv = F,Rowv = F,
            col = mycols(50),
            margins = c(8,8),
            xlab = "Low %Me or negative delta ",
            ylab = "High %Me or positive delta",
            main = paste("Percent intersect for 100kb tiles with lowest and highest 1% methylation,\n",
                         "most negative and most positive first and second deltas in germline development")
)
SV_col = z$CNV_Schiz_Case
hasCNV = z[which(SV_col>0),1:43]
hasnoCNV =  z[which(SV_col==0),1:43]
hasCNV = hasCNV[complete.cases(hasCNV),]
smp =400
hasCNVMethOnly = hasCNV[sample(nrow(hasCNV),smp),11:25]
noCNVMethOnly = hasnoCNV[sample(nrow(hasnoCNV),smp),11:25]

mmee = rbind(colMeans(hasCNV[5:43],na.rm = T),colMeans(hasnoCNV[5:43],na.rm = T))
library(miscTools)
mmeedd = rbind(colMedians(hasCNV[5:43],na.rm = T),colMedians(hasnoCNV[5:43],na.rm = T))
## Organizing the 100kb tiles into those with and without CNV for HSR and Phase 3 mCNV
## Seems to reveal subtle differneces in methyaltion but more dramatic associations with other
## CNV classes
## Apparently where one CNV type is called CNVs of other types will be enriched, That is when
## organizing the data by presence vs. absence of a particular CNV class, the tiles which
## harbor the specific CNV class also carry more CNV from other classes.. (What are the confounders here??)


dev.new()
#min(noCNVMethOnly,na.rm = T)
heatmap.2(as.matrix(noCNVMethOnly),
          na.rm = T,
          trace = 'none',
          Colv = F,Rowv = T,
          col = bluered(30),
          margins = c(8,8),
          xlab = "Cell Type, or stages",
          ylab = "100Kb Tile "#,
          #main = paste("Percent intersect for 100kb tiles with lowest and highest 1% methylation,\n",
          #             "most negative and most positive first and second deltas in germline development")
)


##############################################################################################
##############################################################################################
#KS-Test for distrubution of CNV in %methylated regions
source('H:/Dropbox/BRL/Methylation_Struct_Mut/MyScripts/getTileKS.R', echo=F)
#x = getTile.KS(z,c(6,13))
#x 

colnames(z)
All.KS = data.frame(matrix(data = NA,nrow =500,ncol = 3))
k =1 #initialize counter
for(meth in 5:16){
    for(sv in 17:34){
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
p.val = 0.05
Bonf = nrow(All.KS)
All.KS.sig[as.numeric(All.KS.sig[,2])>(p.val/Bonf),1] = NA

#write.table(x = All.Enrichments1,file = "EnrichmentComparisons1.txt",
#            col.names = T,row.names = F, sep = "\t",quote = F)
dimsr = list(SV_names,Methylome_names)
#sigKS = All.KS[,2]>(.05/144)

r = matrix((as.numeric(All.KS.sig[,1])),nrow = length(SV_names), dimnames = dimsr)
p = matrix((as.numeric(All.KS.sig[,2])),nrow = length(SV_names), dimnames = dimsr)
r[r==0] =NA
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
          cexCol = 1.8,cexRow = 1.5, sepcolor = "black", colsep = 0:(ncol(r)+2),rowsep = 1:nrow(r),
          col = redramp(numBreaks),Colv = F,Rowv = F,
          main = paste("KS TEST D-Statistic for % CpG-methylation\n",
                       "in 100kb Tiles with structural variants vs. rest of genome\n",
                        "Bonferroni-Corrected p < ",p.val,sep = "")
          )



featureTilesHeatmap(featureTilesDF = z,firstSet = 5:25 ,secondSet = 26:43) 

##############################################################################################
##############################################################################################
###TRANSCRIPTION FACTORS AND ANNOTATIONS Vs. DNAme  
#
colnames(z)

for(pct in c(1,5)){
    for(rev in c(T,F)){
    prec =1000
    HILO = NULL
    if(rev == T){HILO = "highest"}
    if(rev== F){HILO = "lowest"}
        Me_TF.Enrichments1 = data.frame(matrix(data = NA,nrow =5000,ncol = 3))
        k =1
        for(meth in 5:25){
            for(tf in 44:146){
                my.foi = c(meth,tf)
                #print(my.foi)}}
                tmp.enr = getTileEnrichment(featureTilesDF = z,my.foi,pct=pct,precision =  prec,rev = rev )
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
        r = matrix((as.numeric(Me_TF.Enrichments1$`fold enrichment`)),nrow = length(TF_names), dimnames = dimsr)
        r[r==0] =1
        r
        dev.new()
        par(bg = "transparent")
        #par(bg = "grey85")
        numBreaks = 256
        pairs.breaks = seq(from=-3.322,to=3.322,length.out=257)
        heatmap.2(log(r,base = 2),na.rm = T,trace = 'none',margins = c(10,18),
                  breaks = pairs.breaks,
                  symbreaks = any(breaks<0,na.rm = T),
                  density.info = 'none', key.xlab = 'Log2 Fold Enrichment', key.title = NA,
                  cexCol = 1.9,cexRow = .5,
                  col = bluered(numBreaks),Colv = F, Rowv = T,
                  main = paste("Enrichment of Transcription Factor Binding\n",
                               "in 100kb Tiles with ",HILO," ",pct,"% methylation",sep ="")
                  )
    }    
}
##########################################################################################
##########################################################################################
colnames(z)

for(pct in c(1)){
    for(rev in c(T,F)){
        prec =1000
        HILO = NULL
        if(rev == T){HILO = "highest"}
        if(rev== F){HILO = "lowest"}
        SV_TF.Enrichments1 = data.frame(matrix(data = NA,nrow =5000,ncol = 3))
        k =1
        for(sv in 26:43){
            for(tf in 44:66){#:146
                my.foi = c(sv,tf)
                #tmp.enr = getTileEnrichment(featureTilesDF = z,my.foi,pct=pct,precision =  prec,rev = rev )
                tmp.enr = getTileEnrichment2.pres(featureTilesDF = z,my.foi,precision =  prec,mode = "pct")
                rownames( SV_TF.Enrichments1 )[k] = tmp.enr[[1]]
                SV_TF.Enrichments1[k,] = c(tmp.enr[[1]],tmp.enr[[2]],tmp.enr[[3]])
                k = k+1
            }
        }
        colnames( SV_TF.Enrichments1 ) = names(tmp.enr)
        SV_TF.Enrichments1 = SV_TF.Enrichments1[complete.cases(SV_TF.Enrichments1),]
        #write.table(x = SV_TF.Enrichments1,file = "TF_vs_SV_EnrichmentComparisons1.txt",
        #            col.names = T,row.names = F, sep = "\t",quote = F)
        dimsr = list(TF_names[1:23],SV_names)
        r = matrix((as.numeric(SV_TF.Enrichments1$`fold enrichment`)),nrow = length(TF_names[1:23]), dimnames = dimsr)
        r[r==0] =1
        r
        dev.new()
        par(bg = "transparent")
        #par(bg = "grey85")
        numBreaks = 256
        pairs.breaks = seq(from=-3.322,to=3.322,length.out=257)
        heatmap.2(log(r,base = 2),na.rm = T,trace = 'none',margins = c(10,18),
                  breaks = pairs.breaks,
                  symbreaks = any(breaks<0,na.rm = T),
                  density.info = 'none', key.xlab = 'Log2 Fold Enrichment', key.title = NA,
                  cexCol = 1.9,cexRow = .5,
                  col = bluered(numBreaks),Colv = F, Rowv = T,
                  main = paste("Enrichment of Transcription Factor Binding in\n",pct,"% ",HILO," CNV enriched ", 
                               "100kb Tiles",sep ="")
        )
    }    
}


#########################################################
#########################################################
#TF by presence/abs
colnames(z)
        SV_TF.Enrichments1 = data.frame(matrix(data = NA,nrow =5000,ncol = 3))
        k =1
        for(sv in 26:43){
            for(tf in 44:146){#:146
                my.foi = c(sv,tf)
                tmp.enr = getTileEnrichment2.pres(featureTilesDF = z,my.foi,precision =  prec,mode = "pct")
                #rownames( SV_TF.Enrichments1 )[k] = tmp.enr[[1]]
                SV_TF.Enrichments1[k,] = c(tmp.enr[[1]])#,tmp.enr[[2]],tmp.enr[[3]])
                k = k+1
            }
        }
        #colnames( SV_TF.Enrichments1 ) = names(tmp.enr)
        SV_TF.Enrichments1 = SV_TF.Enrichments1[complete.cases(SV_TF.Enrichments1),]
        #write.table(x = SV_TF.Enrichments1,file = "TF_vs_SV_EnrichmentComparisons1.txt",
        #            col.names = T,row.names = F, sep = "\t",quote = F)
        dimsr = list(TF_names,SV_names)
        r = matrix((as.numeric(SV_TF.Enrichments1[,1])),nrow = length(TF_names), dimnames = dimsr)
        r[r==0] =1
        r
        dev.new()
        par(bg = "transparent")
        #par(bg = "grey85")
        numBreaks = 256
        pairs.breaks = seq(from=-3.322,to=3.322,length.out=257)
        heatmap.2(log(r,base = 2),na.rm = T,trace = 'none',margins = c(10,18),
                  breaks = pairs.breaks,
                  symbreaks = any(breaks<0,na.rm = T),
                  density.info = 'none', key.xlab = 'Log2 Fold Enrichment', key.title = NA,
                  cexCol = 1.9,cexRow = .5,
                  col = bluered(numBreaks),Colv = F, Rowv = T,
                  main = paste("Enrichment of Transcription Factor Binding in\n", "CNV enriched ", 
                               "100kb Tiles",sep ="")
        )


colnames(z)
SV_TF.Enrichments1 = data.frame(matrix(data = NA,nrow =5000,ncol = 3))
k =1
for(sv in 26:43){
    for(tf in 44:146){#:146
        my.foi = c(tf,sv)
        tmp.enr = getTileEnrichment2.pres(featureTilesDF = z,my.foi,precision =  prec,mode = "pct")
        #rownames( SV_TF.Enrichments1 )[k] = tmp.enr[[1]]
        SV_TF.Enrichments1[k,] = c(tmp.enr[[1]])#,tmp.enr[[2]],tmp.enr[[3]])
        k = k+1
    }
}
#colnames( SV_TF.Enrichments1 ) = names(tmp.enr)
SV_TF.Enrichments1 = SV_TF.Enrichments1[complete.cases(SV_TF.Enrichments1),]
#write.table(x = SV_TF.Enrichments1,file = "TF_vs_SV_EnrichmentComparisons1.txt",
#            col.names = T,row.names = F, sep = "\t",quote = F)
dimsr = list(TF_names,SV_names)
r = matrix((as.numeric(SV_TF.Enrichments1[,1])),nrow = length(TF_names), dimnames = dimsr)
r[r==0] =1

dev.new()
par(bg = "transparent")
#par(bg = "grey85")
numBreaks = 256
pairs.breaks = seq(from=-8,to=8,length.out=257)
heatmap.2(log(r,base = 2),na.rm = T,trace = 'none',margins = c(10,18),
          breaks = pairs.breaks,
          symbreaks = any(breaks<0,na.rm = T),
          density.info = 'none', key.xlab = 'Log2 Fold Enrichment', key.title = NA,
          cexCol = 1.9,cexRow = .5,
          col = bluered(numBreaks),Colv = T, Rowv = T,
          main = paste("Enrichment of CNV in Tiles with Transcription Factor Binding\n","in 100kb Tiles",sep ="")
)

##########################################################################################
##########################################################################################
library(RColorBrewer)
#colnames(z)
#colSums(z[,11:43],na.rm = T)
#Make CNVs unique for case control colunmns 
z[which((z$CNV_Bipolar_Case !=0) & (z$CNV_Bipolar_Control !=0)),c("CNV_Bipolar_Case","CNV_Bipolar_Control")] = 0
z[which((z$CNV_Autism_Case !=0) & (z$CNV_Autism_Control !=0)),c("CNV_Autism_Case","CNV_Autism_Control") ] = 0
z[which((z$CNV_Schiz_Case !=0) & (z$CNV_Schiz_Control !=0)),c("CNV_Schiz_Case","CNV_Schiz_Control") ] = 0
z[which((z$CNV_DevDelay_Case !=0) & (z$CNV_DevDelay_Control !=0)),c("CNV_DevDelay_Case","CNV_DevDelay_Control") ] = 0

rebl = colorRampPalette(c('red',"black"))
grbl    =  colorRampPalette(c('green',"blue"))
ow = colorRampPalette(c('orange',"purple"))
cols = c(rebl(6),grbl(6),ow(9))
pie(rep(1,21),col=cols)
mmm = c()
strt = 0
stp = 100
n = 1
i = n
for(cnvty in rev(c(26:43))){
    dev.new()
    plot(NULL,type = "l",ylim = c(-4,4),xlim = c(0,101/n),ylab = "log2 Fold Enrichment",main = colnames(z)[cnvty] )
    abline(h = c(1,-1),lty = 2)
for (cellty in 17:25){
for(i in seq(strt,(stp-n),n)){
    j = i+n
    #print(c(i,j))
    x = getTileEnrichment2.slice(featureTilesDF = z,foi = c(cellty,cnvty),lo.pct = i,hi.pct = j,precision = 100)
    if (x[[2]][1]==0){ mmm[j] = NA}
    else {mmm[j] =(x[[2]][1])}
    
    mycol = cols[(cellty-4)]
}

sss = matrix(data = mmm,nrow = 10,byrow = T)
sss
ss = as.vector(t(sss))
ss = unlist(ss)
#dev.new()
#plot(ss,type = "l",ylim = c(0,11))

points(log2(ss),type = "l",col =mycol,lwd =2,ylim = c(-10,10))
#points(ss,type = "l",ylim = c(0,11),col ="blue",lwd = 2)
}
}



for(cnvty in rev(c(26,29,32:39))){
    dev.new()
    plot(NULL,type = "l",ylim = c(0,1000),xlim = c(0,101/n),ylab = "CNV Count in methylome quantile",main = colnames(z)[cnvty] )
    abline(h = c(1,-1),lty = 2)
    for (cellty in 17:25){
        for(i in seq(strt,(stp-n),n)){
            j = i+n
            #print(c(i,j))
            x = getTileEnrichment2.sliceSum(featureTilesDF = z,foi = c(cellty,cnvty),lo.pct = i,hi.pct = j,precision = 100)
            #if (x[[2]][1]==0){ mmm[j] = NA}
            mmm[j] =(x[[2]][1])
            
            mycol = cols[(cellty-4)]
        }
        
        sss = matrix(data = mmm,nrow = 10,byrow = T)
        sss
        ss = as.vector(t(sss))
        ss = unlist(ss)
        #dev.new()
        #plot(ss,type = "l",ylim = c(0,11))
        
        points((ss),type = "l",col =mycol,lwd =2)
        #points(ss,type = "l",ylim = c(0,11),col ="blue",lwd = 2)
    }
}
