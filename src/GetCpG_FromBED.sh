#!/bin/bash
##########################################################################################
#This script should take zipped BED files for methylation calls and produce new files
#containing only the CpG positions.

#CpG_Positions are in column10
#Overview: For all files in folder if the file is tab delimited and has fields matching
#FieldsString, Then take all of the entries that correspond to CpG positions with
#Read depth of at least 5 reads.



#Example command that worked for a single file
#zcat GSE63818_PGC_10W_embryo1_M_methylation_calling.bed.gz| awk '$5>4 && $10=="CpG"{print $0}' >../Processed_Data/GSE63818_PGC_10W_embryo1_M_methylation_calling.CpG.bed 
#Now write a loop function for this
awk '$5>4,$10=="CpG"{print $0 }' $inFile > $destFolder CpG_$infile 
