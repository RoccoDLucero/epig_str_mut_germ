#REMEMBER TO USE SETWD  TO THE FOLDER WITH ALL THE .DIST Files
setwd("C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/ROI_Data/BED_Format/Output/SV_samp500/PeakToSV")
#options("scipen"=-100, "digits"=4)
pops = c('AFR','AMR','EUR', 'EAS','SAS')
SVs = c('CNV','DEL','DUP')# , 'INV')
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
MBnumber =2L
##############################################################################################################
##############################################################################################################
MyGetCombinedPlots = function(MBnumber = 1){ 
library(Rmisc)
library(lattice)
library(gridBase)
require(ggplot2)
require(grid)
options(scipen = 0) 
MBnumber = as.integer(MBnumber)
#Set boundaries on PEAK to SV distance
upperLim = as.integer(1000000*MBnumber) #Sets max distance to consider
lowerLim = 1
xmin = 0       
#Plot the distribution of Peak-SV distances for the mean of (here n=10) all samples for each SV-POP
#Later add 95% CI as a pale polygon with appropriate transparency to show all comparisons
pdf(file = paste(MBnumber,"Mb_range_SV_Peak_SV_Plots_AllPOPS_CI.pdf",sep = "") ,onefile = T)
popColors = sapply(c('orange2','red2','green4','navy','darkorchid4'),function(myCol){adjustcolor(myCol,alpha.f = 0.5)})


for (SVtype in SVs){
    #dev.new() #for a quick look
    plotTitle = SVtype
    plot(x = c() ,xlab= "PEAK to SV distance [bp]" ,ylab = "Cumulative proportion",
        main = paste("Population-specific distribution of ranges\nfrom PRDM14 peak to rare ",plotTitle,"s\n(for variants within ",MBnumber,"Mb of a peak)",sep = ""),
        sub =  "\n[95%CI for 10 samples of 500 SVs]",
        xlim = c(xmin,(upperLim+500000)),ylim = c(0,1))
    
    DistribMeans = data.frame(matrix(data = NA,nrow = nrow(All_Dist.Pop.SV[[1]]),ncol = length(pops)))
    myDFpos = grep(SVtype,names(All_Dist.Pop.SV))
    for (pos in 1:length(myDFpos)){
        tmp = All_Dist.Pop.SV[[myDFpos[pos]]][,1:10]
        tmp = tmp[(tmp >= lowerLim)&(tmp<=upperLim),1:10]
        tmp = tmp[complete.cases(tmp),]
        confFram = as.data.frame(matrix(ncol = 3,nrow = nrow(tmp),
                          dimnames = list(c(),c('upper','mean','lower'))))
        for (i in 1:nrow(confFram)){
            confFram[i,] =CI(x = unlist(tmp[i,]),ci = .95)
        }
        fun.ecdf <- ecdf(confFram$lower)
        lower.ecdf = fun.ecdf(sort(confFram$lower))
        fun.ecdf <- ecdf(confFram$upper)
        upper.ecdf = fun.ecdf(sort(confFram$upper))
        polygon(c(confFram$lower,rev(confFram$upper)),
                c(lower.ecdf,rev(upper.ecdf)),
                col = popColors[pos],border = F)
        #lines(ecdf(confFram[,2]))
        fun.ecdf <- ecdf(confFram$mean)
        mean.ecdf = fun.ecdf(sort(confFram$mean))
        DistribMeans[1:length(mean.ecdf),pos] = mean.ecdf
        legend('bottomright',legend = pops,
               col = popColors,lty = c(1),lwd = 3,cex = 1.5,bg = 'white')
    }    
    #before restarting the loop we need to make the ks-test heatmap
    pthresh = (.05/4)
    sigMatrix = matrix(nrow = ncol(DistribMeans),ncol = ncol(DistribMeans))
    DMatrix = matrix(nrow = ncol(DistribMeans),ncol = ncol(DistribMeans))
    for (i in 1:ncol(DistribMeans)){
        for (j in 1:ncol(DistribMeans)){
            tmp1 = ks.test(DistribMeans[,i],DistribMeans[,j])
            sigMatrix[i,j] = (tmp1$p.value<pthresh)
            DMatrix[i,j] = tmp1$statistic
        }
    }
    #qp <- qplot(mpg, wt, data=mtcars)
    #print(qp, vp=viewport(.8, .75, .2, .2))
    #dev.new()
    #plot.new()
    #vp =viewport(x = 500000,y = .2,width = 60000,height = .3)
    #pushViewport(vp)#iewport(x = 500000,y = .2,width = 60000,height = .3))
    #aa = image(z = sigMatrix, col=c("red","blue")[sigMatrix = T])
    #bb = grid(nx=5, ny=5, lty=1)
    #print(aa, vp=vp)
    
    #pushViewport(viewport())
    #xvars <- rnorm(25)
    #yvars <- rnorm(25)
    #xyplot(yvars~xvars)
    #pushViewport(viewport(x=.6,y=.8,width=.25,height=.25,just=c("left","top")))
    #grid.rect()
    #par(plt = gridPLT(), new=TRUE)
    #plot(xvars,yvars)
    #popViewport(2)
    
    print(sigMatrix)
}

dev.off()
}
################################################################################################################
################################################################################################################