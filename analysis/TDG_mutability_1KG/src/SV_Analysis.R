#Prelinary Processing...
#Import the 1000Genomes Phase3 SV calls for all individuals as 'trimmed_phase3_sv_map' (The VCF with Metatdata Header lines removed) 

#Make a list of carriers and non-carriers of TDG G-199-S in the phase 3 1000 Genomes data:
  #Use vcf exported from 1000 Genomes Browser at TDG locus at rs4135094 To generate TDG_Calls.
  #by looking for all individuals with genotype not "0" or "0|0"
  #by looking for all individuals with genotype "0" or "0|0"
  #by identifying individuals without data at rs4135094
  
  #At some later point we need to ensure that all of the individuals included in the Carrier group have functional TDG alleles
  #This is less important for the Non-Carrier group since we are presupposing a gain of function in the carriers.

#Superpopulations: all SVs from 1000Genomes Phase 3
#For HG38
#trimmed_phase3_sv_map = read.delim("C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf/trimmed_phase3_sv_map.vcf", stringsAsFactors=FALSE)
#For HG19
trimmed_phase3_sv_map = read.delim("../Raw_Data/ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf", stringsAsFactors=FALSE)

#Trim down the SV calls to only include Deletions Duplications Inversions and mCNVs:
#Remove the row if SVTYPE== {"ALU","LINE1","SVA","INS","NUMT"} in the 'INFO' column.
#Select Rows by Mutatiion typ DEL DUP INV CNV
SV_Subset_DEL <- trimmed_phase3_sv_map[grep("SVTYPE=DEL", trimmed_phase3_sv_map$INFO), ]
SV_Subset_DUP <- trimmed_phase3_sv_map[grep("SVTYPE=DUP", trimmed_phase3_sv_map$INFO), ]
SV_Subset_INV <- trimmed_phase3_sv_map[grep("SVTYPE=INV", trimmed_phase3_sv_map$INFO), ]
SV_Subset_CNV <- trimmed_phase3_sv_map[grep("SVTYPE=CNV", trimmed_phase3_sv_map$INFO), ]


#Make lists of Carriers and Non-Carriers of TDGG199S in phase3 1000Genomes dataset.
#Note that only called individuals are included. In the larger data set of SVs some individuals will be included that lacked data for rs4135094
#This must be dealt with later when pulling from the total SV set in the context of TDG status....
  #Use vcf exported from 1000 Genomes Browser at TDG locus at rs4135094.
  #Carrieres will be het and possibly homo for G199S Allele
TDG_Variants_phase3 = read.delim("C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/TDG_Variants_phase3.vcf")
TDG_Calls = TDG_Variants_phase3[10:length(TDG_Variants_phase3)]
Carriers = TDG_Calls[which(TDG_Calls != '0|0')]
NonCarriers = TDG_Calls[which(TDG_Calls == '0|0')]
HomoCarriers = Carriers[which(Carriers == '1|1')]

write.table(colnames(Carriers),col.names = F, file = paste('C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/AlleleStatus/',
                                   'TDG_G199S_Carriers',".txt",sep = ''),sep = "\t",quote = F,row.names = F)

#Compare number of duplications called per individual based on TDG allele status:
##CarrierDups = data.frame(SV_Subset_DUP[,intersect(colnames(Carriers),colnames(trimmed_phase3_sv_map))],stringsAsFactors = T)
CarrierDups = SV_Subset_DUP[,intersect(colnames(Carriers),colnames(trimmed_phase3_sv_map))]
NonCarrierDups = SV_Subset_DUP[,intersect(colnames(NonCarriers),colnames(trimmed_phase3_sv_map))]

Samples =  unlist(colnames(CarrierDups))
SV_Counts = c()
for(Sample in Samples){
  count = 0
  Loci = unlist(CarrierDups[Sample])
  for(locus in Loci){
    if(locus != "0|0"){
      if(locus != "0"){
        count = count+1
      }
    }
  }
  SV_Counts = c(SV_Counts,count)
}
CarrierDupCounts = SV_Counts

Samples =  unlist(colnames(NonCarrierDups))
SV_Counts = c()
for(Sample in Samples){
  count = 0
  Loci = unlist(NonCarrierDups[Sample])
  for(locus in Loci){
    if(locus != "0|0"){
      if(locus != "0"){
        count = count+1
      }
    }
  }
  SV_Counts = c(SV_Counts,count)
}
NonCarrierDupCounts = SV_Counts

plot.new()
#hist(sample(NonCarrierDupCounts,563), xlim=c(0,90), ylim = c(0,1000), col=rgb(.1,.1,.1,.5),breaks = 40)
hist(NonCarrierDupCounts, xlim=c(0,90), ylim = c(0,400), col=rgb(.1,.1,.1,.5),breaks = 40)
hist(CarrierDupCounts, col =rgb(.8,0,1,.5),breaks = 40 , add = T)
t.test(NonCarrierDupCounts,CarrierDupCounts)

#Compare number of deletions called per individual based on TDG allele status:
CarrierDels = SV_Subset_DEL[,intersect(colnames(Carriers),colnames(trimmed_phase3_sv_map))]
NonCarrierDels = SV_Subset_DEL[,intersect(colnames(NonCarriers),colnames(trimmed_phase3_sv_map))]

Samples =  unlist(colnames(CarrierDels))
SV_Counts = c()
for(Sample in Samples){
  count = 0
  Loci = unlist(CarrierDels[Sample])
  for(locus in Loci){
    if(locus != "0|0"){
      if(locus != "0"){
        count = count+1
      }
    }
  }
  SV_Counts = c(SV_Counts,count)
}
CarrierDelCounts = SV_Counts

Samples =  unlist(colnames(NonCarrierDels))
SV_Counts = c()
for(Sample in Samples){
  count = 0
  Loci = unlist(NonCarrierDels[Sample])
  for(locus in Loci){
    if(locus != "0|0"){
      if(locus != "0"){
        count = count+1
      }
    }
  }
  SV_Counts = c(SV_Counts,count)
}
NonCarrierDelCounts = SV_Counts

plot.new()
#hist(sample(NonCarrierDupCounts,563), xlim=c(0,90), ylim = c(0,1000), col=rgb(.1,.1,.1,.5),breaks = 40)
hist(NonCarrierDelCounts, ylim = c(0,400), col=rgb(.1,.1,.1,.5),breaks = 40)
hist(CarrierDelCounts, col =rgb(.8,0,1,.5),breaks = 40 , add = T)
t.test(NonCarrierDelCounts,CarrierDelCounts) 

#########################################################################################################
#Compare number of deletions called per individual based on TDG allele status:
CarrierCNVs = SV_Subset_CNV[,intersect(colnames(Carriers),colnames(trimmed_phase3_sv_map))]
NonCarrierCNVs = SV_Subset_CNV[,intersect(colnames(NonCarriers),colnames(trimmed_phase3_sv_map))]

Samples =  unlist(colnames(CarrierCNVs))
SV_Counts = c()
for(Sample in Samples){
  count = 0
  Loci = unlist(CarrierCNVs[Sample])
  for(locus in Loci){
    if(locus != "0|0"){
      if(locus != "0"){
        count = count+1
      }
    }
  }
  SV_Counts = c(SV_Counts,count)
}
CarrierCNVCounts = SV_Counts

Samples =  unlist(colnames(NonCarrierCNVs))
SV_Counts = c()
for(Sample in Samples){
  count = 0
  Loci = unlist(NonCarrierCNVs[Sample])
  for(locus in Loci){
    if(locus != "0|0"){
      if(locus != "0"){
        count = count+1
      }
    }
  }
  SV_Counts = c(SV_Counts,count)
}
NonCarrierCNVCounts = SV_Counts

plot.new()
#hist(sample(NonCarrierDupCounts,563), xlim=c(0,90), ylim = c(0,1000), col=rgb(.1,.1,.1,.5),breaks = 40)
hist(NonCarrierCNVCounts, ylim = c(0,150), col=rgb(.1,.1,.1,.5),breaks = 40)
hist(CarrierCNVCounts, col =rgb(.8,0,1,.5),breaks = 40 , add = T)
t.test(NonCarrierCNVCounts,CarrierCNVCounts) 

#Possible Differences could be due to population History...Need to check the trios etc. for de Novo rate...
#Need to subset by superpopulation and correct for population bias...select evenly from each superpop where possible...
