#For all PRDM14 ChIP SEQ peaks: take chromStart +/- 10kb as window:
#Compare start and end of SV to see if it is in the window
require('gtools')
require('hash')

##EVENTUALLY WRITE THESE IN FUNCTION FORM  FOR REUSE:

 #Turn this into an import statement for all of the populations and all of the SV classes.
trial_no = 1
#IMPLEMENT TESTS AT SEVERAL WINDOW SIZES AND OUTPUT THE COUNTS IN A VECTOR
#Import PRDM14 ChIP-Seq PEAKS and make various ranges:
PRDM14_PEAKS = read.delim("C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/ROI_Data/PEAKS/hg38_PRDM14_PEAKS.bed", stringsAsFactors=FALSE,header = F)
colnames(PRDM14_PEAKS) = c("chrom","chromStart","chromEnd")
ChIP_WindowSizeVec = c(10000)
PRDM14PeakWindows = as.data.frame(cbind(PRDM14_PEAKS$chromStart-(ChIP_WindowSizeVec[1]/2),PRDM14_PEAKS$chromStart+(ChIP_WindowSizeVec[1]/2)))
PRDM14PeakWindows = as.data.frame(cbind(PRDM14_PEAKS$chrom,PRDM14PeakWindows),stringsAsFactors = F)
rm(PRDM14_PEAKS)
PRDM14PeakWindows[,1] = as.integer(gsub("chr",'',x = PRDM14PeakWindows[,1]))
PRDM14PeakWindows = PRDM14PeakWindows[complete.cases(PRDM14PeakWindows),]
samp_portion = 1.0
randPeaks = sample(rownames(PRDM14PeakWindows),size = (nrow(PRDM14PeakWindows)*samp_portion))

PopulationNamesVec = c('AFR','AMR','EAS','EUR','SAS') 
SVTYPE_Vec = c('CNV','DEL','DUP','INV')
CHIP_PRDM14_POP_SC_DF = data.frame(POP_SV = character(),HIT_COUNT = integer(), BinSize = integer(), stringsAsFactors = F)
for(i in 1:length(PopulationNamesVec)){
  for(j in 1:length(SVTYPE_Vec)){
    POP_SV_NAME = paste(PopulationNamesVec[i],SVTYPE_Vec[j],sep = '')
    POP_SVTYPE = read.delim(paste('C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/ROI_Data/',POP_SV_NAME,".txt",sep = ''),stringsAsFactors=FALSE)
    POP_SVTYPE = cbind(POP_SVTYPE$chrom,POP_SVTYPE$chromStart,POP_SVTYPE$chromEND)
    POP_SVTYPE = cbind(as.integer(POP_SVTYPE[,1]),as.integer(POP_SVTYPE[,2]),as.integer(POP_SVTYPE[,3]))
    POP_SVTYPE = as.data.frame(POP_SVTYPE,stringsAsFactors = F)
    
    POP_SVTYPE =(POP_SVTYPE[with(POP_SVTYPE,order(POP_SVTYPE[1]),),])
    POP_SVTYPE = POP_SVTYPE[complete.cases(POP_SVTYPE),]
    
    count = 0
    binSize= ChIP_WindowSizeVec[1] ##Later add the for loop to process multiple window sizes
    for(peak in 1:nrow(PRDM14PeakWindows[randPeaks,])){   #For each ChipSeq peak-window
      for(SV in 1:nrow(POP_SVTYPE[,])){                     #Test every SV at both ends for inclusion in the window
        if(PRDM14PeakWindows[peak,1] == POP_SVTYPE[SV,1]){     #For events occurring on the same chromosome
          if(PRDM14PeakWindows[peak,2] <= POP_SVTYPE[SV,2]){   #Check that the left window coordinate is less than BP coordinate
            if(PRDM14PeakWindows[peak,3] >= POP_SVTYPE[SV,2]){ #And that the right window coordinate is greater than the BP coordinate 
              #print(".")
              count = count+1
            }
          }  
          if(PRDM14PeakWindows[peak,2] <= POP_SVTYPE[SV,3]){   #Check that the left window coordinate is less than BP coordinate
            if(PRDM14PeakWindows[peak,3] >= POP_SVTYPE[SV,3]){ #And that the right window coordinate is greater than the BP coordinate 
              #print(".")
              count = count+1    
            }                   
          }    
        }  
      }
    }#teachers[nrow(teachers) + 1, ] <- c( "ted", 50)
    CHIP_PRDM14_POP_SC_DF[nrow(CHIP_PRDM14_POP_SC_DF)+1,] = c(POP_SV_NAME,count,binSize)
  }
  
}
write.table(CHIP_PRDM14_POP_SC_DF, file = paste('C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/ROI_Data/',
                                                'RARE_SVs_PRDM14_ChIP_HITS_',as.character(binSize),'_Proportion_',
                                                as.character(samp_portion),'_trial_',as.character(trial_no),".txt",sep = ''),sep = "\t",quote = F,row.names = F)


#IMPLEMENT TESTS AT SEVERAL WINDOW SIZES AND OUTPUT THE COUNTS IN A VECTOR

#For now use hand typed table of SV counts:

Rare_SV_total_counts = matrix(data = c(1043,13235,1827,200,729,8942,1066,156,823,9149,1471,187,897,8675,1608,165,763,8804,1614,173),nrow = 5,byrow = T)
colnames(Rare_SV_total_counts) = c('CNV','DEL','DUP','INV')
rownames(Rare_SV_total_counts) = c('AFR','AMR','EAS','EUR','SAS')

RARE_SVs_PRDM14_ChIP_HITS_10000kb = matrix(data = c(37,616,96,8,23,425,47,4,39,441,73,11,38,371,61,5,22,349,96,11),nrow = 5,byrow = T)
colnames(RARE_SVs_PRDM14_ChIP_HITS_10000kb) = c('CNV','DEL','DUP','INV')
rownames(RARE_SVs_PRDM14_ChIP_HITS_10000kb) = c('AFR','AMR','EAS','EUR','SAS')

Rare_SV_total_counts = Rare_SV_total_counts-RARE_SVs_PRDM14_ChIP_HITS_10000kb

#Perform all populationwise comparisons within CNV_Types
#POP_SV_Fisher_Result = matrix(data = NA,nrow = 40,ncol = 3)
POP_SV_Fisher_Result = data.frame("Comparison" = character(), "Odds.Ratio" = numeric(),"p.Value" =numeric(),stringsAsFactors = F)
currentRow =1
for(cnv in 1:ncol(RARE_SVs_PRDM14_ChIP_HITS_10000kb)){
  SV_TYPE = colnames(RARE_SVs_PRDM14_ChIP_HITS_10000kb)[cnv]
  for(pop1 in 1:(nrow(RARE_SVs_PRDM14_ChIP_HITS_10000kb)-1)){
    for(pop2 in (pop1+1):(nrow(RARE_SVs_PRDM14_ChIP_HITS_10000kb))){
    #a = c(i,j)
    a = matrix(c(RARE_SVs_PRDM14_ChIP_HITS_10000kb[pop1,cnv],RARE_SVs_PRDM14_ChIP_HITS_10000kb[pop2,cnv],Rare_SV_total_counts[pop1,cnv],Rare_SV_total_counts[pop2,cnv]),nrow = 2,byrow = T)
    tested = fisher.test(a)
    comparison = (paste(rownames(RARE_SVs_PRDM14_ChIP_HITS_10000kb)[pop1],rownames(RARE_SVs_PRDM14_ChIP_HITS_10000kb)[pop2],SV_TYPE,sep = '.'))
    POP_SV_Fisher_Result[currentRow,1]= comparison
    POP_SV_Fisher_Result[currentRow,2] = tested$estimate[[1]]
    POP_SV_Fisher_Result[currentRow,3] = tested$p.value
    currentRow = currentRow+1
    }
  }
}
rownames(x = POP_SV_Fisher_Result) = POP_SV_Fisher_Result$Comparison
POP_SV_Fisher_Result = POP_SV_Fisher_Result[,2:3]
POP_SV_Fisher_Result[which(POP_SV_Fisher_Result[,2]<=.05),]
POP_SV_Fisher_Result

POP_SV_Chi_Result = data.frame("Comparison" = character(), "ChiSq Statistic" = numeric(),"p.Value" =numeric(),stringsAsFactors = F)
currentRow =1
for(cnv in 1:ncol(RARE_SVs_PRDM14_ChIP_HITS_10000kb)){
  SV_TYPE = colnames(RARE_SVs_PRDM14_ChIP_HITS_10000kb)[cnv]
  for(pop1 in 1:(nrow(RARE_SVs_PRDM14_ChIP_HITS_10000kb)-1)){
    for(pop2 in (pop1+1):(nrow(RARE_SVs_PRDM14_ChIP_HITS_10000kb))){
      #a = c(i,j)
      a = matrix(c(RARE_SVs_PRDM14_ChIP_HITS_10000kb[pop1,cnv],RARE_SVs_PRDM14_ChIP_HITS_10000kb[pop2,cnv],Rare_SV_total_counts[pop1,cnv],Rare_SV_total_counts[pop2,cnv]),nrow = 2,byrow = T)
      tested = chisq.test(a)
      comparison = (paste(rownames(RARE_SVs_PRDM14_ChIP_HITS_10000kb)[pop1],rownames(RARE_SVs_PRDM14_ChIP_HITS_10000kb)[pop2],SV_TYPE,sep = '.'))
      POP_SV_Chi_Result[currentRow,1]= comparison
      POP_SV_Chi_Result[currentRow,2] = tested$statistic[[1]]
      POP_SV_Chi_Result[currentRow,3] = tested$p.value
      currentRow = currentRow+1
    }
  }
}
POP_SV_Chi_Result[which(POP_SV_Chi_Result[,3]<=.05),]

##MAke sure chromosome lengths are not exceeded...
##ENSURE that windows around peaks are not overlapping...
##FIND A WAY TO CORRECT FOR TWO BP HITS IN THE SAME PEAK INTERVAL

