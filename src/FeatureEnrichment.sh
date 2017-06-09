
#Methylation Bedfiles are *.Me.bed
#Structural Variant BedFiles are *.SV.bed
#ChIP-seq peaks are *.ChIP.bed


#Set -a as rhe Me file, -b as the SV file  
for file_a in *Me.bed; do for file_b in *.SV.bed; do printf "%s\t%s\t" "$file_a" "$file_b"; intersectBed -a $file_a  -b $file_b |wc -l; done;done > Me_SV_counts.txt

#Create the set of randomly shuffled files
# or i in $(seq 3:100); do 
for file in *.bed; do shuffleBed -i $file -g ../hg38.chrom.sizes > ./ShuffleBeds/shuff.$file; done
#Run with shuffled methylome windows 
 for file_a in ./ShuffleBeds/s*Me.bed; do for file_b in *.SV.bed; do printf "%s\t%s\t" "$file_a" "$file_b"; intersectBed -a $file_a  -b $file_b |wc -l; done;done > shuffle_Me_SV_counts.txt

#Generate files containing all SVs of a given SV type
for a in *DUP*; do cat $a >> tmpfoo.tmp;done

rm AllSV_tmp.bed
for a in *SV.bed; do cat $a >>AllSV_tmp.bed;done
for file in *Me.bed; do printf "AllSVs\t%s\t" "$file"; intersectBed -a $file -b AllSV_tmp.bed|wc -l;done


#Idea : Run closestBED to identify the SV closest to each PRDM14 peak. Record this in an ordered Hitlist of SVs assigned to each each peak.
#  For each prdm14 remove all of the SVs in its Hitlist and re-run closestBed. Do this for the nearest 20 SVs or all of the SVs within ~10kb.
# This will tell us if certain PRDM14 peaks associate very strongly with SV annotations and what type of SVs tend to be closest to the peaks.
# We can then calculate the average distance from each peak to its nearest 1,2,...n  SVs.

