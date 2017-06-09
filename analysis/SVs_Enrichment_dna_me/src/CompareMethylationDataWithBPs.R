
#Preprocess the methyation BED files to make compatible with this comparison script... 
AcceptableChrNames = c(as.character((1:23)),paste("chr",as.character((1:23)),sep = ''))

hg38_BCM_Lowest1 <- read.delim("C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/ROI_Data/Methylation/hg38_BCM_Lowest1.bed", header=FALSE, stringsAsFactors=FALSE)
hg38_BCM_Lowest1 = hg38_BCM_Lowest1[which(is.element(el = hg38_BCM_Lowest1[,1],set = AcceptableChrNames)),1:3]
hg38_BCM_Lowest1[,1] = c(gsub(pattern = "chr",replacement = "",x = hg38_BCM_Lowest1[,1]))
row.names(hg38_BCM_Lowest1) = NULL

hg38_BCM_nonMD1 <- read.delim("C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/ROI_Data/Methylation/hg38_BCM_nonMD1.bed", header=FALSE, stringsAsFactors=FALSE)
hg38_BCM_nonMD1 = hg38_BCM_nonMD1[which(is.element(el = hg38_BCM_nonMD1[,1],set = AcceptableChrNames)),1:3]
hg38_BCM_nonMD1[,1] = c(gsub(pattern = "chr",replacement = "",x = hg38_BCM_nonMD1[,1]))
row.names(hg38_BCM_nonMD1) = NULL
##Write the above as a function "preprocess Me BED"
#function


CurrentWindows = hg38_BCM_nonMD1
#CurrentWindows = hg38_BCM_Lowest1

PopulationNamesVec = c('AFR','AMR','EAS','EUR','SAS') 
SVTYPE_Vec = c('CNV','DEL','DUP','INV')
DNAMe_POP_SV_DF = data.frame(POP_SV = character(),HIT_COUNT = integer(), stringsAsFactors = F)
for(i in 1:length(PopulationNamesVec)){
  for(j in 1:length(SVTYPE_Vec)){
    POP_SV_NAME = paste(PopulationNamesVec[i],SVTYPE_Vec[j],sep = '_')
    POP_SVTYPE = read.delim(paste('C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/ROI_Data/','RARE_',POP_SV_NAME,".txt",sep = ''),stringsAsFactors=FALSE)
    POP_SVTYPE = cbind(POP_SVTYPE$chrom,POP_SVTYPE$chromStart,POP_SVTYPE$chromEND)
    POP_SVTYPE = cbind(as.integer(POP_SVTYPE[,1]),as.integer(POP_SVTYPE[,2]),as.integer(POP_SVTYPE[,3]))
    POP_SVTYPE = as.data.frame(POP_SVTYPE,stringsAsFactors = F)
    
    POP_SVTYPE =(POP_SVTYPE[with(POP_SVTYPE,order(POP_SVTYPE[1]),),])
    POP_SVTYPE = POP_SVTYPE[complete.cases(POP_SVTYPE),]
    #These POP_SVTYPE Lists should be generated once and written to a file......    
    
    count = 0
    for(DNAMeWin in 1:nrow(CurrentWindows)){   #For each methylation window (sizes Vary)
      for(SV in 1:nrow(POP_SVTYPE[,])){                     #Test every SV at both ends for inclusion in the window
        if(CurrentWindows[DNAMeWin,1] == POP_SVTYPE[SV,1]){     #For events occurring on the same chromosome
          if(CurrentWindows[DNAMeWin,2] <= POP_SVTYPE[SV,2]){   #Check that the left window coordinate is less than BP coordinate
            if(CurrentWindows[DNAMeWin,3] >= POP_SVTYPE[SV,2]){ #And that the right window coordinate is greater than the BP coordinate 
              count = count+1
            }
          }  
          if(CurrentWindows[DNAMeWin,2] <= POP_SVTYPE[SV,3]){   #Check that the left window coordinate is less than BP coordinate
            if(CurrentWindows[DNAMeWin,3] >= POP_SVTYPE[SV,3]){ #And that the right window coordinate is greater than the BP coordinate 
              count = count+1    
            }                   
          }    
        }  
      }
    }
    DNAMe_POP_SV_DF[nrow(DNAMe_POP_SV_DF)+1,] = c(POP_SV_NAME,count)
  }
  
}
write.table(DNAMe_POP_SV_DF, file = paste('C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/ROI_Data/',
                             'RARE_SVs_NormalMe_HITS',".txt",sep = ''),sep = "\t",quote = F,row.names = F)

#Also need to normalize for number od individuals sampled per population

#IMPLEMENT TESTS AT SEVERAL WINDOW SIZES AND OUTPUT THE COUNTS IN A VECTOR

RARE_SVs_nonMD_HITS_First_Try <- read.delim("C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/ROI_Data/RARE_SVs_NormalMe_HITS_First_Try.txt", stringsAsFactors=FALSE)
RARE_SVs_HypoMe_HITS <- read.delim("C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/ROI_Data/RARE_SVs_HypoMe_HITS.txt", stringsAsFactors=FALSE)

Rare_SV_total_counts_nonMD = matrix(data = RARE_SVs_nonMD_HITS_First_Try[1:12,2],nrow = 3,byrow = T)
colnames(Rare_SV_total_counts_nonMD) = c('CNV','DEL','DUP','INV')
rownames(Rare_SV_total_counts_nonMD) = c('AFR','AMR','EAS')#,'EUR','SAS')

Rare_SV_total_counts_hypoMe = matrix(data = RARE_SVs_HypoMe_HITS[,2],nrow = 5,byrow = T)
colnames(Rare_SV_total_counts_hypoMe) = c('CNV','DEL','DUP','INV')
rownames(Rare_SV_total_counts_hypoMe) = c('AFR','AMR','EAS','EUR','SAS')
Rare_SV_total_counts_hypoMe = Rare_SV_total_counts_hypoMe[1:3,]#Remove later

#Perform all populationwise comparisons within CNV_Types
#POP_SV_Fisher_Result = matrix(data = NA,nrow = 40,ncol = 3)
POP_SV_Fisher_Result = data.frame("Comparison" = character(), "Odds.Ratio" = numeric(),"p.Value" =numeric(),stringsAsFactors = F)
currentRow = 1
for(cnv in 1:ncol(Rare_SV_total_counts_hypoMe)){
  SV_TYPE = colnames(Rare_SV_total_counts_hypoMe)[cnv]
  for(pop in 1:(nrow(Rare_SV_total_counts_hypoMe)-1)){
    for(pop2 in (pop1+1):(nrow(Rare_SV_total_counts_hypoMe))){
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


##MAke sure chromosome lengths are not exceeded...
##ENSURE that windows around peaks are not overlapping...
##FIND A WAY TO CORRECT FOR TWO BP HITS IN THE SAME PEAK INTERVAL