
#############################################################################################################
############################################################################################################
#Function Definitions#

MyMakeList =  function(df, cols){
    outList = list()
    for(col in cols){outList[[col]] =df[,col]}
    return(outList)
}

myTrimList = function(list,uppr,lowr = 0){
    outList = list()
    for (i in a){outList[[1+length(outList)]] = i[which(i>=lowr&i<=uppr)]}
    return(outList)
} 

############################################################################################################
############################################################################################################
setwd("C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/ROI_Data/BED_Format/Output/SV_samp500/PeakToSV")
pops = c('AFR','AMR','EUR', 'EAS','SAS')
SVs = c('CNV','DEL','DUP')#, 'INV') #Inversions discluded due to reduced frequency.Sampling
# method from previous (Linux) processing must be modified to n~100
popSVs = sapply(pops,paste,SVs,sep='_')

SV.shuffStatus = c('+', '(shuff.)+')
PRDM14.shuffStatus = c('TruePRDM14', 'ShuffPRDM14')
MyFilePatterns = sapply(SV.shuffStatus,paste,popSVs,sep='')
MyFilePatterns = unlist(sapply(MyFilePatterns,paste,PRDM14.shuffStatus,sep='.*'),recursive = T)
MyFilePatterns = rbind(MyFilePatterns[1:2,1:length(popSVs)],MyFilePatterns[1:2,(1+length(popSVs)):(2*length(popSVs))])

#For every column in the MyFilePatterns matrix create a read the corresponding files
#into a 6968x(4*samplingdepth) matrix; samplingdepth = 10 for now
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
#Set boundaries on PEAK to SV distance 
MBnumber = 5L
upperLim =1000000*MBnumber
lowerLim = 1
xmin = 0




##############################################################################################################
##############################################################################################################
library(Rmisc)
#Plot the distribution of Peak-SV distances for the mean of (here n=10) all samples for each SV-POP
#Later add 95% CI as a pale polygon with appropriate transparency to show all comparisons
MyGetPlots = function(MBnumber){
    options(scipen = 0) 
    #MBnumber= 10L# as.integer(MBnumber)
    upperLim = 1000000*MBnumber
    lowerLim = 0
    xmin = 0
    pdf(file = paste(MBnumber,"MbSV_Peak_POP-SV_Plots_meanCI.pdf",sep = '') ,onefile = T)
    lapply(All_Dist.Pop.SV, function(df){  #For each DataFrame plot the appropriate statistics
        plotTitle = substr(colnames(df)[1],start = 1,stop = 12)
        dev.new() #for testing
        df = All_Dist.Pop.SV[[1]] #for testing
        plot(x = c() ,xlab= "PEAK to SV distance [bp]" ,ylab = "Cumulative proportion",
             main = paste("Distribution of distances\nfrom PRDM14 peak to",plotTitle,
                          "\n(for variants within ",MBnumber," Mb of a peak)",sep = " "),
             sub =  "\n[mean and 95%CI for 10 samples of 500 SVs]",
             xlim = c(xmin,(upperLim+500000)),ylim = c(0,1))
        lowerLim = 0
        upperLim = 20000000
        for(rng in c(1:10,11:20,21:30.31:40)){
            xxo = MyMakeList(df,cols = rng)
            xxo = myTrimList(list = xxo,uppr = upperLim,lowr = lowerLim)
            xxo =unlist(xxo)
            lines(ecdf(xxo))
        }
        #foo = list()
        #for (i in xxo){foo[1+length(foo)] = CI(i,1:10]),ci = .95)
        #}
        
        
        xxx = df
        #Select the rows for which all values are
        
        #xxo = xxx[complete.cases(xxx[(xxx >= lowerLim)&(xxx<=upperLim),1:10]),1:10]
        xxp = xxx[complete.cases(xxx[(xxx >= lowerLim)&(xxx<=upperLim),21:30]),21:30]
        xxq = xxx[complete.cases(xxx[(xxx >= lowerLim)&(xxx<=upperLim),31:40]),31:40]
        #xxr = xxx[complete.cases(xxx[(xxx >= lowerLim)&(xxx<=upperLim),11:20]),11:20]
        
        #Write a function to draw the mean and conf int polygon from the code below....
        confFram = as.data.frame(matrix(ncol = 3,nrow = nrow(xxo),dimnames = list(c(),c('upper','mean','lower'))))
        for (i in 1:nrow(confFram)){
            confFram[i,] =CI(x = unlist(xxo[i,1:10]),ci = .95)
        }
        fun.ecdf <- ecdf(confFram$lower)
        lower.ecdf = fun.ecdf(sort(confFram$lower))
        fun.ecdf <- ecdf(confFram$upper)
        upper.ecdf = fun.ecdf(sort(confFram$upper))
        polygon(c(confFram$lower,rev(confFram$upper)),
                c(lower.ecdf,rev(upper.ecdf)),
                col = 'grey',border = F)
        lines(ecdf(confFram[,2]))
        truSVtruPK = confFram$mean
        
        confFram = as.data.frame(matrix(ncol = 3,nrow = nrow(xxp),dimnames = list(c(),c('upper','mean','lower'))))
        for (i in 1:nrow(confFram)){
            confFram[i,] =CI(x = unlist(xxp[i,1:10]),ci = .95)
        }
        fun.ecdf <- ecdf(confFram$lower)
        lower.ecdf = fun.ecdf(sort(confFram$lower))
        fun.ecdf <- ecdf(confFram$upper)
        upper.ecdf = fun.ecdf(sort(confFram$upper))
        polygon(c(confFram$lower,rev(confFram$upper)),
                c(lower.ecdf,rev(upper.ecdf)),
                col = (adjustcolor('firebrick1',alpha.f = 0.5)),border = F)
        lines(ecdf(confFram[,2]),col = 'firebrick4')
        truSVshfPK = confFram$mean
        
        confFram = as.data.frame(matrix(ncol = 3,nrow = nrow(xxq),dimnames = list(c(),c('upper','mean','lower'))))
        for (i in 1:nrow(confFram)){
            confFram[i,] =CI(x = unlist(xxq[i,1:10]),ci = .95)
        }
        fun.ecdf <- ecdf(confFram$lower)
        lower.ecdf = fun.ecdf(sort(confFram$lower))
        fun.ecdf <- ecdf(confFram$upper)
        upper.ecdf = fun.ecdf(sort(confFram$upper))
        polygon(c(confFram$lower,rev(confFram$upper)),
                c(lower.ecdf,rev(upper.ecdf)),
                col = (adjustcolor('darkorchid1',alpha.f = 0.5)),border = F)
        lines(ecdf(confFram[,2]),col = 'darkorchid4')
        shfSVshfPK = confFram$mean
        
        
        #Perform KS-test for mean distribution functions:
        options("scipen"=-100, "digits"=4)
        a = ks.test(x = truSVtruPK,y = truSVshfPK)
        b = ks.test(x = truSVtruPK,y = shfSVshfPK)
        text(x = (.8*(upperLim+500000)),y = .25,labels = paste("KS-test:\nred-black:",
                                                               signif(a$p.value,3),"\npurple-black:",signif(b$p.value,4)))
        
        
        
        legend('bottomright',legend = c("TrueSV-TruePRDM14","TrueSV-ShuffPRDM14","ShuffSV-ShuffPRDM14"),
               col = c('black','firebrick','darkorchid'),lty = c(1,1,1),lwd=2,cex = 1,bg = 'white')
    }
    ) 
    dev.off()
}
################################################################################################################
################################################################################################################
MyGetPlots(30)

MyMakeList =  function(df, cols){
   outList = list()
    for(col in cols){outList[[col]] =df[,col]}
    return(outList)
}

myTrimList = function(list,uppr,lowr = 0){
    outList = list()
    print(length(outList))
    for (i in a){outList[[1+length(outList)]] = i[which(i>=lowr&i<=uppr)]}
    return(outList)
} 


