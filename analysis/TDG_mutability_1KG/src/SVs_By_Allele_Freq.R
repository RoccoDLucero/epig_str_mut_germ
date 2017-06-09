#This script is written to generate ROIs based on 1000Genomes Phase3 Structural Variants
#mapped to the modified GRCh38 genome.
#Methylation data will need to be mapped to this modified GRCh38 in a separate script



#Import the 1000Genomes Phase3 SV calls for all individuals as 'trimmed_phase3_sv_map' (The VCF with Metatdata Header lines removed) 
#SV_Subset_DEL <- trimmed_phase3_sv_map[grep("SVTYPE=DEL", trimmed_phase3_sv_map$INFO), ]
#SV_Subset_DUP <- trimmed_phase3_sv_map[grep("SVTYPE=DUP", trimmed_phase3_sv_map$INFO), ]
#SV_Subset_INV <- trimmed_phase3_sv_map[grep("SVTYPE=INV", trimmed_phase3_sv_map$INFO), ]
#SV_Subset_CNV <- trimmed_phase3_sv_map[grep("SVTYPE=CNV", trimmed_phase3_sv_map$INFO), ]

#Later define a function that takes as input an SV_TYPE, a population ID, and a frequency threshold to return
#The rare, or unique structural variants as list of positions and individuals with each SV. 

Allele_Freqs_DUP =  as.data.frame(SV_Subset_DUP[,8],row.names = rownames(SV_Subset_DUP),stringsAsFactors = F) #Column 8 is INFO
#INFO column contains allele frequency information by population
#Rare East Asian DUPS
Rare_EAS_DUPS = data.frame('EAS_AF' = numeric(nrow(SV_Subset_DUP)))
rownames(Rare_EAS_DUPS) = rownames(Allele_Freqs_DUP)

for(i in 1:nrow(Allele_Freqs_DUP)){
  a = unlist(strsplit(Allele_Freqs_DUP[i,1],split = c(";"))) #Cannot assume order of elements of fixed number in info column must read strings
  a = a[grep("EAS_AF=",a)]
  #print(a)
  a = unlist(strsplit(a,"="))
  a = unlist(strsplit(a,","))
  #print(a)
  
  b = as.numeric(a[3])
  c = as.numeric(a[4])
  d = as.numeric(a[5])
  e = as.numeric(a[6])
  a = as.numeric(a[2])
  Rare_EAS_DUPS[i,1] = a
  Rare_EAS_DUPS[i,2] = b
  Rare_EAS_DUPS[i,3] = c
  Rare_EAS_DUPS[i,4] = d
  Rare_EAS_DUPS[i,5] = e
    }
#Rare_EAS_DUPS[1:500,]

#Next Select all of the SV indices for rows with at least one non-zero value: These are the EAS_DUPS
Rare_EAS_DUPS[Rare_EAS_DUPS==0] = NA #Convert all zeros to NA so we can remove SVs not-present in the population
Rare_EAS_DUPS = Rare_EAS_DUPS[apply(!is.na(Rare_EAS_DUPS),1,any),] #Select all of the SV indices for rows with at least one non-NA value:
#Rare_EAS_DUPS[is.na(Rare_EAS_DUPS)] = 99 

#Next Select all of the SV indices where the AF is below 1%

#Rare_EAS_DUPS_1_PCT = (Rare_EAS_DUPS[apply(which(Rare_EAS_DUPS<=0.01),1,all),])
Rare_EAS_DUPS_1_PCT = Rare_EAS_DUPS[which(Rare_EAS_DUPS<=0.01),]
Rare_EAS_DUPS_1_PCT = (Rare_EAS_DUPS_1_PCT[apply(!is.na(Rare_EAS_DUPS_1_PCT),1,any),]) #For some reason certain rows are lost...and become all NA. Mysterious

#This is a list of the SV row positions based on the PHASE3 call set data frame
#This information will allow is to Characterize rare variants by population