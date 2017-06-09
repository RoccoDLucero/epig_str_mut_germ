#Run within the ../BED_Format folder
cd /mnt/hgfs/Dropbox_Shared/BRL/Methylation_Struct_Mut/ROI_Data/BED_Format
#For all the 'TRUE'RARE_SV files, sample 500 SVs and then get the distance from every 'TRUE' PRDM14 Peak to the nearest SV in the sampled set 
for ((i=1; i<=10; i++));
do for SV in RARE*.SV.bed;
do shuf -n 500 $SV| closestBed -a hg38_PRDM14_PEAKS_sorted.ChIP.bed -b stdin -d |cut -f7 > ./Output/SV_samp500/PeakToSV/$SV.500sampled$i.TruePRDM14.dist;
done;
done

#For all the 'TRUE'RARE_SV files, sample 500 SVs and then get the distance from every 'SHUFFLED' PRDM14 Peak to the nearest SV in the sampled set 
for ((i=1; i<=10; i++));
do for SV in RARE*.SV.bed;
do shuf -n 500 $SV| closestBed -a ./ShuffleBeds/shuff.hg38_PRDM14_PEAKS_sorted.ChIP.bed -b stdin -d |cut -f7 > ./Output/SV_samp500/PeakToSV/$SV.500sampled$i.ShuffPRDM14.dist;
done;
done

###
###REVERSE THE ORIENTATION OF THE DISTANCE CALCULATION SO NOW WE ARE ORIENTED ON THE SVs

#For all the 'TRUE'RARE_SV files, sample 500 SVs and then get the distance from every SV in the sampled set to the nearest 'TRUE' PRDM14 Peak 
for ((i=1; i<=10; i++));
do for SV in RARE*.SV.bed;
do shuf -n 500 $SV| closestBed -b hg38_PRDM14_PEAKS_sorted.ChIP.bed -a stdin -d |cut -f7 > ./Output/SV_samp500/SVToPEAK/$SV.500sampled$i.TruePRDM14.dist;
done;
done

#For all the 'TRUE'RARE_SV files, sample 500 SVs and then get the distance from every SV in the sampled set to the nearest 'SHUFFLED' PRDM14 Peak
for ((i=1; i<=10; i++));
do for SV in RARE*.SV.bed;
do shuf -n 500 $SV| closestBed -b ./ShuffleBeds/shuff.hg38_PRDM14_PEAKS_sorted.ChIP.bed -a stdin -d |cut -f7 > ./Output/SV_samp500/SVToPEAK/$SV.500sampled$i.ShuffPRDM14.dist;
done;
done
########################################################################
########################################################################
cd ShuffleBeds

#For all the 'SHUFFLED' RARE_SV files, sample 500 SVs and then get the distance from every 'TRUE' PRDM14 Peak to the nearest SV in the sampled set 
for ((i=1; i<=10; i++));
do for SV in shuff.RARE*.SV.bed;
do shuf -n 500 $SV| closestBed -a ../hg38_PRDM14_PEAKS_sorted.ChIP.bed -b stdin -d |cut -f7 > ../Output/SV_samp500/PeakToSV/$SV.500sampled$i.TruePRDM14.dist;
done;
done

#For all the 'SHUFFLED' RARE_SV files, sample 500 SVs and then get the distance from every 'SHUFFLED' PRDM14 Peak to the nearest SV in the sampled set 
for ((i=1; i<=10; i++));
do for SV in shuff.RARE*.SV.bed;
do shuf -n 500 $SV| closestBed -a ../ShuffleBeds/shuff.hg38_PRDM14_PEAKS_sorted.ChIP.bed -b stdin -d |cut -f7 > ../Output/SV_samp500/PeakToSV/$SV.500sampled$i.ShuffPRDM14.dist;
done;
done

###
###REVERSE THE ORIENTATION OF THE DISTANCE CALCULATION SO NOW WE ARE ORIENTED ON THE SVs

#For all the 'SHUFFLED' RARE_SV files, sample 500 SVs and then get the distance from every SV in the sampled set to the nearest 'TRUE' PRDM14 Peak
for ((i=1; i<=10; i++));
do for SV in shuff.RARE*.SV.bed;
do shuf -n 500 $SV| closestBed -b ../hg38_PRDM14_PEAKS_sorted.ChIP.bed -a stdin -d |cut -f7 > ../Output/SV_samp500/SVToPEAK/$SV.500sampled$i.TruePRDM14.dist;
done;
done

#For all the 'SHUFFLED' RARE_SV files, sample 500 SVs and then get the distance from every SV in the sampled set to the nearest 'SHUFFLED' PRDM14 Peak 
for ((i=1; i<=10; i++));
do for SV in shuff.RARE*.SV.bed;
do shuf -n 500 $SV| closestBed -b ../ShuffleBeds/shuff.hg38_PRDM14_PEAKS_sorted.ChIP.bed -a stdin -d |cut -f7 > ../Output/SV_samp500/SVToPEAK/$SV.500sampled$i.ShuffPRDM14.dist;
done;
done