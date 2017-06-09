for file in *.bed; do wc $file  >> SVs_Total.txt;done
for file in *.bed; do printf $file printf "\t" >> SVs_PRDM14Peaks_10000bpWindow_Overlaps1.txt ; windowBed -a hg38_PRDM14_PEAKS_num_sorted.bed -w 10000 -b $file |wc >> SVs_PRDM14Peaks_10000bpWindow_Overlaps1.txt ;done
