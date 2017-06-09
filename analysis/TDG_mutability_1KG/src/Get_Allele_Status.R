#Should be a generic function for getting the allele status from the 1000Genomes or Hapmap SNP calls

#Make a list of carriers and non-carriers of some allele of interest from the phase 3 1000 Genomes data:
#Use vcf exported from 1000 Genomes Browser at given locus, e.g. at rs4135094, to generate TDG_Calls.
#by looking for all individuals with genotype not "0" or "0|0" as "Carriers"
#by looking for all individuals with genotype "0" or "0|0" as "Non-Carriers"

#At some later point we need to ensure that all of the individuals included each group
#have e.g. functional alleles or meet other important criteria that must be computed for these individuals
#by incorporating more of the available data.


#Make lists of Carriers and Non-Carriers of TDGG199S in phase3 1000Genomes dataset.
#Note that only called individuals are included. In the larger data set of SVs some individuals will be included that lacked data for rs4135094
#This must be dealt with later when pulling from the total SV set in the context of TDG status....
#Use vcf exported from 1000 Genomes Browser at TDG locus at rs4135094.
#Carrieres will be het and possibly homo for G199S Allele

##CurrentLocusVCF = #Get from ARGS or other script...should be a VCF filename under the proper directory
  CurrentLocusVCF = "TDG_Variants_phase3.vcf"  
  #Later this should be expanded to take any input/output path... 
VariantsPhase3 = read.delim(paste("C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/VariantCalls/",CurrentLocusVCF,sep = ''))
VarCalls = VariantsPhase3[10:length(VariantsPhase3)]
Carriers = VarCalls[which(VarCalls != '0|0')]
NonCarriers = VarCalls[which(VarCalls == '0|0')]
HomoCarriers = Carriers[which(Carriers == '1|1')]

AlleleName = "TDG_G199S"
write.table(colnames(NonCarriers),col.names = F, file = paste('C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/AlleleStatus/',AlleleName,
                                                           '_NonCarriers',".txt",sep = ''),sep = "\t",quote = F,row.names = F)

write.table(colnames(Carriers),col.names = F, file = paste('C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/AlleleStatus/',AlleleName,
                                                           '_Carriers',".txt",sep = ''),sep = "\t",quote = F,row.names = F)

write.table(colnames(HomoCarriers),col.names = F, file = paste('C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/AlleleStatus/',AlleleName,
                                                           '_HomoCarriers',".txt",sep = ''),sep = "\t",quote = F,row.names = F)
