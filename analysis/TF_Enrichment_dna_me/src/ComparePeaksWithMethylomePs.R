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
#Center the window on the chip-seq peak
PRDM14PeakWindows = as.data.frame(cbind(PRDM14_PEAKS$chrom,PRDM14PeakWindows),stringsAsFactors = F)
rm(PRDM14_PEAKS)
PRDM14PeakWindows[,1] = as.integer(gsub("chr",'',x = PRDM14PeakWindows[,1]))
PRDM14PeakWindows = PRDM14PeakWindows[complete.cases(PRDM14PeakWindows),]
samp_portion = 1.0
randPeaks = sample(rownames(PRDM14PeakWindows),size = (nrow(PRDM14PeakWindows)*samp_portion))

#hg38_Molaro_Lowest1 <- read.delim("C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/ROI_Data/Methylation/hg38_Molaro_Lowest1.bed", header=FALSE)
#hg38_Molaro_Lowest5 <- read.delim("C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/ROI_Data/Methylation/hg38_Molaro_Lowest5.bed", header=FALSE)
hg38_BCM_Lowest1 <- read.delim("C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/ROI_Data/Methylation/hg38_BCM_Lowest1.bed", header=FALSE)
hg38_BCM_Lowest5  <- read.delim("C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/ROI_Data/Methylation/hg38_BCM_Lowest5.bed", header=FALSE)
CurrentMethylome = hg38_BCM_Lowest5 


DNAMe_PRDM14_DF = data.frame(Methylome = character(),HIT_COUNT = integer(), stringsAsFactors = F)


count = 0
binSize= ChIP_WindowSizeVec[1]
for(peak in 1:nrow(PRDM14PeakWindows[randPeaks,])){   #For each ChipSeq peak-window
  for(MeBin in 1:nrow(CurrentMethylome[,])){                     #Test every Methyaltion window at both ends for inclusion in the ChIP-seq window
    if(PRDM14PeakWindows[peak,1] == CurrentMethylome[MeBin,1]){     #For events occurring on the same chromosome
      if(PRDM14PeakWindows[peak,2] <= CurrentMethylome[MeBin,2]){   #Check that the left window coordinate is less than left methylome coordinate
        if(PRDM14PeakWindows[peak,3] >= CurrentMethylome[MeBin,2]){ #And that the right window coordinate is greater than the left methylome coordinate 
          #print(".")
          count = count+1
        }##Need to add logic that prevents double counting...perhaps introduce new boolean variable 'overlaping' that once True increases "count" and moves on to next iteration 
      }  
      if(PRDM14PeakWindows[peak,2] <= CurrentMethylome[MeBin,3]){   #Check that the left window coordinate is less than BP coordinate
        if(PRDM14PeakWindows[peak,3] >= CurrentMethylome[MeBin,3]){ #And that the right window coordinate is greater than the BP coordinate 
          #print(".")
          count = count+1    
        }                   
      }    
    }  
  }
}

DNAMe_PRDM14_DF[nrow(DNAMe_PRDM14_DF)+1,] = c("BCM_lowest_5",count)



#UPDATE write operation
#write.table(CHIP_PRDM14_POP_SC_DF, file = paste('C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/ROI_Data/',
#                                                'RARE_SVs_PRDM14_ChIP_HITS_',as.character(binSize),'_Proportion_',
#                                                as.character(samp_portion),'_trial_',as.character(trial_no),".txt",sep = ''),sep = "\t",quote = F,row.names = F)
