#Prelinary Processing...
#Make a list of carriers and non-carriers of TDG G-199-S in the phase 3 1000 Genomes data:
#Use vcf exported from 1000 Genomes Browser at TDG locus at rs4135094.
#Look for all individuals with genotype not "0" or "0|0"

AfroPop <- read.delim("C:/Users/Rocco/Dropbox/BRL/Methylation_Struct_Mut/AfroPop.txt", stringsAsFactors=FALSE)
#This is for the contintal African population excluding ACB and ASW

#Trim down the SV calls to only include Deletions Duplications Inversions and mCNVs:
#Remove the row if SVTYPE== {"ALU","LINE1","SVA","INS","NUMT"} in the 'INFO' column.
#Select Rows by Mutatiion typ DEL DUP INV CNV
SV_Subset_DEL <- trimmed_phase3_sv_map[grep("SVTYPE=DEL", trimmed_phase3_sv_map$INFO), ]
SV_Subset_DUP <- trimmed_phase3_sv_map[grep("SVTYPE=DUP", trimmed_phase3_sv_map$INFO), ]
SV_Subset_INV <- trimmed_phase3_sv_map[grep("SVTYPE=INV", trimmed_phase3_sv_map$INFO), ]
SV_Subset_CNV <- trimmed_phase3_sv_map[grep("SVTYPE=CNV", trimmed_phase3_sv_map$INFO), ]


#Make lists of Carriers and Non-Carriers of TDGG199S in phase3 1000Genomes dataset.
#Use vcf exported from 1000 Genomes Browser at TDG locus at rs4135094.
#Carrieres will be het and possibly homo for G199S Allele
TDG_Calls = TDG_Variants_phase3[10:length(TDG_Variants_phase3)]
Carriers = TDG_Calls[which(TDG_Calls != '0|0')]
NonCarriers = TDG_Calls[which(TDG_Calls == '0|0')]

AFRO_Carriers = Carriers[intersect(colnames(Carriers),AfroPop[,1])]
AFRO_NonCarriers = NonCarriers[intersect(colnames(NonCarriers),AfroPop[,1])]


#Compare number of duplications called per individual based on TDG allele status:
##CarrierDups = data.frame(SV_Subset_DUP[,intersect(colnames(Carriers),colnames(trimmed_phase3_sv_map))],stringsAsFactors = T)
CarrierDups = SV_Subset_DUP[,intersect(colnames(AFRO_Carriers),colnames(trimmed_phase3_sv_map))]
NonCarrierDups = SV_Subset_DUP[,intersect(colnames(AFRO_NonCarriers),colnames(trimmed_phase3_sv_map))]

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

#plot.new()
#hist(sample(NonCarrierDupCounts,563), xlim=c(0,90), ylim = c(0,1000), col=rgb(.1,.1,.1,.5),breaks = 40)
hist(NonCarrierDupCounts, xlim=c(0,50), ylim = c(0,40), col=rgb(.1,.1,.1,.5),breaks = 40)
hist(CarrierDupCounts, col =rgb(.8,0,1,.5),breaks = 40 , add = T)
t.test(NonCarrierDupCounts,CarrierDupCounts)

#Compare number of deletions called per individual based on TDG allele status:
CarrierDels = SV_Subset_DEL[,intersect(colnames(AFRO_Carriers),colnames(trimmed_phase3_sv_map))]
NonCarrierDels = SV_Subset_DEL[,intersect(colnames(AFRO_NonCarriers),colnames(trimmed_phase3_sv_map))]

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

#plot.new()
#hist(sample(NonCarrierDupCounts,563), xlim=c(0,90), ylim = c(0,1000), col=rgb(.1,.1,.1,.5),breaks = 40)
hist(NonCarrierDelCounts, ylim = c(0,30), col=rgb(.1,.1,.1,.5),breaks = 40)
hist(CarrierDelCounts, col =rgb(.8,0,1,.5),breaks = 40 , add = T)
t.test(NonCarrierDelCounts,CarrierDelCounts) 

#########################################################################################################
#Compare number of deletions called per individual based on TDG allele status:
CarrierCNVs = SV_Subset_CNV[,intersect(colnames(AFRO_Carriers),colnames(trimmed_phase3_sv_map))]
NonCarrierCNVs = SV_Subset_CNV[,intersect(colnames(AFRO_NonCarriers),colnames(trimmed_phase3_sv_map))]

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

#plot.new()
#hist(sample(NonCarrierDupCounts,563), xlim=c(0,90), ylim = c(0,1000), col=rgb(.1,.1,.1,.5),breaks = 40)
hist(NonCarrierCNVCounts, ylim = c(0,30), col=rgb(.1,.1,.1,.5),breaks = 40)
hist(CarrierCNVCounts, col =rgb(.8,0,1,.5),breaks = 40 , add = T)
t.test(NonCarrierCNVCounts,CarrierCNVCounts) 

#Possible Differences could be due to population History...Need to check the trios etc. for de Novo rate...
#Need to subset by superpopulation and correct for population bias...select evenly from each superpop where possible...
