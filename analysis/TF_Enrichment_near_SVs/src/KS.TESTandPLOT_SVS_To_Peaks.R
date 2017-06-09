#REMEMBER TO USE SETWD  TO THE FOLDER WITH ALL THE .DIST Files
setwd("C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/ROI_Data/BED_Format/Output/SV_samp500/SVToPEAK")
pops = c('AFR','AMR','EUR', 'EAS','SAS')
SVs = c('CNV','DEL','DUP')# , 'INV')
popSVs = sapply(pops,paste,SVs,sep='_')

SV.shuffStatus = c('+', '(shuff.)+')
PRDM14.shuffStatus = c('TruePRDM14', 'ShuffPRDM14')
MyFilePatterns = sapply(SV.shuffStatus,paste,popSVs,sep='')
MyFilePatterns = unlist(sapply(MyFilePatterns,paste,PRDM14.shuffStatus,sep='.*'),recursive = T)
MyFilePatterns = rbind(MyFilePatterns[1:2,1:length(popSVs)],MyFilePatterns[1:2,(1+length(popSVs)):(2*length(popSVs))])

#For every column in the MyFilePatterns matrix create a  
#read the corresponding files #into a 6968x(4*samplingdepth) ; samplingdepth = 10 for now
All_Dist.Pop.SV = vector(mode = "list", length = length(popSVs))
for (popSV in 1:ncol(MyFilePatterns)){
  Dist.Pop.SV = matrix(nrow = 6968)
  colnames(Dist.Pop.SV) = "junk"
  for(MyPattern in MyFilePatterns[,popSV] ){
    tmpList = list.files(pattern = MyPattern) #get all files
    for(distFile in tmpList){                 #For each file stick its data columnwise to the rest of the data for files matching the popSV  
      tmpDat = sort((read.table(distFile)[,1]))
      Dist.Pop.SV = cbind(Dist.Pop.SV,tmpDat[1:6968])           #At the end Dist.Pop.SV should be a 6968*40 matrix
      colnames(Dist.Pop.SV)[ncol(Dist.Pop.SV)] = distFile
      rm(tmpDat)
    }
    Dist.Pop.SV = data.frame(Dist.Pop.SV[,1:ncol(Dist.Pop.SV)])
    All_Dist.Pop.SV[[popSV]] = Dist.Pop.SV
  }
}
names(All_Dist.Pop.SV) = popSVs
All_Dist.Pop.SV = lapply(All_Dist.Pop.SV,function(dat){dat[,2:ncol(dat)]}) #Remove Null Column
#Set boundaries on PEAK to SV distance ...probably apply this at the stage of plotting rather than here...
upperLim = 1000000
lowerLim = 000000
xmin = 0
xmax = 40000000

pdf(file = "testPlots.pdf" ,onefile = T)
dev.new()
lapply(All_Dist.Pop.SV, function(df){
    plotTitle = substr(colnames(df)[1],start = 1,stop = 12)
    plot(x = c() ,xlab= "PEAK to SV distance [bp]" ,ylab = "Cumulative proportion",
         main = paste("Distribution of distances\nfrom PRDM14 peak to",plotTitle,"\nfor variants within 2Mb",sep = " "),
         sub =  "\n[10 samples of 500 SVs per comparison]",
         xlim = c(xmin,(upperLim+500000)),ylim = c(0,1))
    
    lapply(colnames(df)[1:10],function(x){
        lines(ecdf(df[(df>=lowerLim)&(df<upperLim),x]),col = 'black',pch ='.')})  #TrueSV   TruePRDM14
    #lapply(colnames(df)[11:20],function(x){
    #    lines(ecdf(df[(df>=lowerLim)&(df<upperLim),x]),col = 'blue',pch ='.')})  #ShuffSV  TruePRDM14
    lapply(colnames(df)[21:30],function(x){
        lines(ecdf(df[(df>=lowerLim)&(df<upperLim),x]),col = 'red3',pch ='.')})  #TrueSV   ShuffPRDM14 
    lapply(colnames(df)[31:40],function(x){
        lines(ecdf(df[(df>=lowerLim)&(df<upperLim),x]),col = 'green3',pch ='.')})#ShuffSV  ShuffPRDM14
    
    legend('bottomright',legend = c("TrueSV-TruePRDM14","TrueSV-ShuffPRDM14","ShuffSV-ShuffPRDM14"),
           col = c('black','red3','green3'),lty = c(1,1,1),cex = .8,bg = 'white')
    }
) 
dev.off()

KSTESTS = lapply(All_Dist.Pop.SV, function(df){ 
    KS_List_True.v.PeakShuff = vector(mode = "list", length = 100)
    for(n in seq(1:10)){
        for(m in seq(1:10)){
            KS_List_True.v.PeakShuff[[((10*(n-1))+m)]] = ks.test(x = df[,n], y = df[,(m+20)])
        }
        
    }
   return(KS_List_True.v.PeakShuff) 
})
#Need to get an estimate of pval and D for all tests.... 
aa = lapply(KSTESTS,function(df){
    a = c()
    b = c()
    for (m in 1:100){
        a = c(a,KSTESTS[[df]][[m]]$p.value)
        b = c(b,KSTESTS[[df]][[m]]$statistic)
    }
    return(c(mean(a)*10, mean(b)))
})

#The goal here is to plot SV-Peak distances for each SVtype with all pops on one plot
pdf(file = "PoPPlots.pdf" ,onefile = T)
#for CNVs
dev.new()
dev.new()
plotTitle = substr(colnames(df)[1],start = 1,stop = 12)
#pdf(file = "testPlots.pdf" ,onefile = T)
plot(x = c() ,xlab= "PEAK to SV distance [bp]" ,ylab = "Cumulative proportion",
     main = paste("Distribution of distances\nfrom PRDM14 peak to",plotTitle,"\nfor variants within 200Mb",sep = " "),
     sub =  "\n[10 samples of 500 SVs per comparison]",
     xlim = c(xmin,(upperLim+500000)),ylim = c(0,1))
lapply(All_Dist.Pop.SV[grep("DUP",names(All_Dist.Pop.SV))], function(df){
   lapply(colnames(df)[1:10],function(x){
        lines(ecdf(df[(df>=lowerLim)&(df<upperLim),x]),col = 'black',pch ='.')})  #TrueSV   TruePRDM14
    #lapply(colnames(df)[11:20],function(x){
    #    lines(ecdf(df[(df>=lowerLim)&(df<upperLim),x]),col = 'blue',pch ='.')})  #ShuffSV  TruePRDM14
    lapply(colnames(df)[21:30],function(x){
        lines(ecdf(df[(df>=lowerLim)&(df<upperLim),x]),col = 'firebrick',pch ='.')})  #TrueSV   ShuffPRDM14 
    lapply(colnames(df)[31:40],function(x){
        lines(ecdf(df[(df>=lowerLim)&(df<upperLim),x]),col = 'darkorchid',pch ='.')})#ShuffSV  ShuffPRDM14
}
)
legend('bottomright',legend = c("TrueSV-TruePRDM14","TrueSV-ShuffPRDM14","ShuffSV-ShuffPRDM14"),
       col = c('black','red3','green3'),lty = c(1,1,1),cex = .8,bg = 'white')

dev.off()

