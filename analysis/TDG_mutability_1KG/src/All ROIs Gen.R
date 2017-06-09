#Sample of INFO Column Data"SVTYPE=CNV;END=755966;CS=DUP_gs;AC=3,206;AF=0.00059904,0.04113419;NS=2504;
#   AN=5008;EAS_AF=0.001,0.0615;EUR_AF=0.001,0.0417;AFR_AF=0.0,0.0303;AMR_AF=0.0014,0.0259;SAS_AF=0.0,0.045;SITEPOST=1.0"

require(hash)
  #Import the 1000Genomes Phase3 SV calls for all individuals as 'trimmed_phase3_sv_map' (The VCF with Metatdata Header lines removed)
trimmed_phase3_sv_map = read.delim("C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf/trimmed_phase3_sv_map.vcf", stringsAsFactors=FALSE)
  #Later set locations in filesystem as variables...

popID_FilesVec = c("AfroPop.txt", "AMR_POP.txt","AsianPop.txt","EURO_POP.txt", "SAS_POP.txt")
PopulationNamesVec = c('AFR','AMR','EAS','EUR','SAS') 
popHash= hash(popID_FilesVec,PopulationNamesVec) #Confirm that Appropriate Samples and population will be paired
SVTYPE_Vec = c('DUP','CNV','DEL','INV')

#Separate the data set by population:
for(popFile in popID_FilesVec){
  #Load the list of SampleIDs pertaining to each population in the VCF File 
  currentPopIDs = read.delim(paste("C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/Populations/",popFile,sep ='' ), header=FALSE, stringsAsFactors=FALSE)
  validPopIDs = intersect(currentPopIDs[,1],colnames(trimmed_phase3_sv_map))
  rm(currentPopIDs)
  
  currentPopAllSVs =  cbind(trimmed_phase3_sv_map[,1:8],trimmed_phase3_sv_map[,validPopIDs]) 
  #Now we have a list of all SVS per population, trim down by type 
  #Set up data frames to handle each SV type, populate the DF and then export to a text file...
  for(SVTYPE in SVTYPE_Vec){
    #For The current population identify the Rare SVs in the current population of each SVTYPE in SVTYPE_VEC
    SV_Subset = currentPopAllSVs[grep(paste("SVTYPE=", SVTYPE,sep = ""), currentPopAllSVs$INFO), ] #Use the INFO column to get all SVs of a given type
    #Now find the population specifc allele frequency for each allele...
    currentAF_Var = paste(popHash[[popFile]],'_AF',sep = '')
    Rare_SV = data.frame('Phase3_SV_Locus'= numeric(nrow(SV_Subset))) #Empty DF to hold the data #currentAF_Var
    rownames(Rare_SV) = rownames(SV_Subset) #preserve the rownames of the original 1000Genomes VCF
    for(i in 1:nrow(SV_Subset)){
      a = unlist(strsplit(SV_Subset[i,8],split = c(";"))) #Cannot assume order of elements of fixed number in info column must read strings
      a = a[grep(currentAF_Var,a)]
      #print(a)
      a = unlist(strsplit(a,"="))
      a = unlist(strsplit(a,","))
      #print(a)
      
      b = as.numeric(a[3])
      c = as.numeric(a[4])
      d = as.numeric(a[5])
      e = as.numeric(a[6])
      a = as.numeric(a[2])
      Rare_SV[i,1] = a
      Rare_SV[i,2] = b
      Rare_SV[i,3] = c
      Rare_SV[i,4] = d
      Rare_SV[i,5] = e
    }
    #Next Select all of the SV indices for rows with at least one non-zero value: These are the EAS_DUPS
    Rare_SV[Rare_SV==0] = NA #Convert all zeros to NA so we can remove SVs not-present in the population
    Rare_SV = Rare_SV[apply(!is.na(Rare_SV),1,any),] #Select all of the SV indices for rows with at least one non-NA value:
    Rare_SV = (Rare_SV[rowSums(is.na(Rare_SV)) != 5 ,])
    #Next Select all of the SV indices where the AF is below 1%
    #Rare_EAS_DUPS_1_PCT = (Rare_EAS_DUPS[apply(which(Rare_EAS_DUPS<=0.01),1,all),])
    Rare_SV_1_PCT = Rare_SV[which(Rare_SV<=0.01),]
    Rare_SV_1_PCT = (Rare_SV_1_PCT[rowSums(is.na(Rare_SV_1_PCT)) != 5 ,])
    Rare_SV_1_PCT = format(Rare_SV_1_PCT, scientific = F)
    #This is a list of the SV row positions based on the PHASE3 call set data frame
    #This information will allow is to Characterize rare variants by population
    Rare_SV_1_PCT = currentPopAllSVs[rownames(Rare_SV_1_PCT),c(1,2,8)]
    for(i in 1:nrow(Rare_SV_1_PCT)){
      a = unlist(strsplit(Rare_SV_1_PCT[i,"INFO"],split = c(";"))) #Cannot assume order of elements of fixed number in info column must read strings
      a = a[grep("^END=",a)] #END=
      #print(a)
      a = unlist(strsplit(a,"="))
      a = unlist(strsplit(a,","))
      #print(a)
      
      a = as.numeric(a[2])
      Rare_SV_1_PCT[i,4] = a
      }
    colnames(Rare_SV_1_PCT) = c('chrom','chromStart',"INFO" ,'chromEND')
    class = rep("Structural Variant", nrow(Rare_SV_1_PCT))
    name = rep(popHash[[popFile]],nrow(Rare_SV_1_PCT))
    type = rep("1000 Gen Phase 3", nrow(Rare_SV_1_PCT))
    subtype = rep(SVTYPE, nrow(Rare_SV_1_PCT))
    strand = rep("+", nrow(Rare_SV_1_PCT))
    phase = rep(".", nrow(Rare_SV_1_PCT))
    score = rep(1.0, nrow(Rare_SV_1_PCT))
    Rare_SV_1_PCT = cbind(class,name,type,subtype,Rare_SV_1_PCT[,c(1,2,4)],strand,phase,score)
    Rare_SV_1_PCT = format(Rare_SV_1_PCT, scientific = F)
    #Convert to LFF Format for GENBOREE
    
    #Export to a file or variable...
    write.table(Rare_SV_1_PCT, file = paste('C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/ROI_Data/','RARE_',popHash[[popFile]],'_',SVTYPE,".txt",sep = ''),sep = "\t",quote = F,row.names = F)
    #write.table(Rare_SV_1_PCT, file = paste('C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/ROI_Data/',popHash[[popFile]],SVTYPE,"_with_rownames.txt"),sep = "\t",quote = F,row.names = T)
}
  
  #ADD code to create a table for each population and RARE SV type total counts to use for comparison with other count
  #information such as the PRDM14 Chip-SEQ enrichment
  
}

