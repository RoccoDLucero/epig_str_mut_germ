

plot(ecdf(foo[,"RARE_AFR_CNV.SV.bed.dist"]),col = 'black',lty = 3)
lines(ecdf(RARE_AFR_CNV.SV.bed.10KshuffleChIP[1:69970,1]),col = 'red',lty = 3)
lines(ecdf(foo[,"RARE_AMR_CNV.SV.bed.dist"]),col = 'darkred')
lines(ecdf(foo[,"RARE_EUR_CNV.SV.bed.dist"]),col = 'darkblue')
lines(ecdf(foo[,"RARE_EAS_CNV.SV.bed.dist"]),col = 'darkgreen')
lines(ecdf(foo[,"RARE_SAS_CNV.SV.bed.dist"]),col = 'purple')
##For Rare CNVs AFR seems to be closest to random distribution and the other populations
##Show a shift to a greater distance **NOTE:Double check AFR_RareCNVs
lines(ecdf(shuffleChIPfoo[,"RARE_AFR_CNV.SV.bed.shuffleChIP.dist"]),col = 'grey', lty =2)
lines(ecdf(shuffleChIPfoo[,"RARE_AMR_CNV.SV.bed.shuffleChIP.dist"]),col = 'pink', lty =2)
lines(ecdf(shuffleChIPfoo[,"RARE_EUR_CNV.SV.bed.shuffleChIP.dist"]),col = 'lightblue', lty =2)
lines(ecdf(shuffleChIPfoo[,"RARE_EAS_CNV.SV.bed.shuffleChIP.dist"]),col = 'lightgreen', lty =2)
lines(ecdf(shuffleChIPfoo[,"RARE_SAS_CNV.SV.bed.shuffleChIP.dist"]),col = 'magenta', lty =2)


plot(ecdf(foo[,"RARE_AFR_INV.SV.bed.dist"]),col = 'black')
lines(ecdf(foo[,"RARE_AMR_INV.SV.bed.dist"]),col = 'red')
lines(ecdf(foo[,"RARE_EUR_INV.SV.bed.dist"]),col = 'blue')
lines(ecdf(foo[,"RARE_EAS_INV.SV.bed.dist"]),col = 'green')
lines(ecdf(foo[,"RARE_SAS_INV.SV.bed.dist"]),col = 'purple')
lines(ecdf(shuffleChIPfoo[,"RARE_AFR_INV.SV.bed.shuffleChIP.dist"]),col = 'grey', lty =2)
lines(ecdf(shuffleChIPfoo[,"RARE_AMR_INV.SV.bed.shuffleChIP.dist"]),col = 'pink', lty =2)
lines(ecdf(shuffleChIPfoo[,"RARE_EUR_INV.SV.bed.shuffleChIP.dist"]),col = 'lightblue', lty =2)
lines(ecdf(shuffleChIPfoo[,"RARE_EAS_INV.SV.bed.shuffleChIP.dist"]),col = 'lightgreen', lty =2)
lines(ecdf(shuffleChIPfoo[,"RARE_SAS_INV.SV.bed.shuffleChIP.dist"]),col = 'magenta', lty =2)

plot(ecdf(foo[,"RARE_AFR_DUP.SV.bed.dist"]),col = 'black')
lines(ecdf(foo[,"RARE_AMR_DUP.SV.bed.dist"]),col = 'red')
lines(ecdf(foo[,"RARE_EUR_DUP.SV.bed.dist"]),col = 'blue')
lines(ecdf(foo[,"RARE_EAS_DUP.SV.bed.dist"]),col = 'green')
lines(ecdf(foo[,"RARE_SAS_DUP.SV.bed.dist"]),col = 'purple')
lines(ecdf(shuffleChIPfoo[,"RARE_AFR_DUP.SV.bed.shuffleChIP.dist"]),col = 'grey', lty =2)
lines(ecdf(shuffleChIPfoo[,"RARE_AMR_DUP.SV.bed.shuffleChIP.dist"]),col = 'pink', lty =2)
lines(ecdf(shuffleChIPfoo[,"RARE_EUR_DUP.SV.bed.shuffleChIP.dist"]),col = 'lightblue', lty =2)
lines(ecdf(shuffleChIPfoo[,"RARE_EAS_DUP.SV.bed.shuffleChIP.dist"]),col = 'lightgreen', lty =2)
lines(ecdf(shuffleChIPfoo[,"RARE_SAS_DUP.SV.bed.shuffleChIP.dist"]),col = 'magenta', lty =2)

plot(ecdf(foo[,"RARE_AFR_DEL.SV.bed.dist"]),col = 'black')
lines(ecdf(foo[,"RARE_AMR_DEL.SV.bed.dist"]),col = 'red')
lines(ecdf(foo[,"RARE_EUR_DEL.SV.bed.dist"]),col = 'blue')
lines(ecdf(foo[,"RARE_EAS_DEL.SV.bed.dist"]),col = 'green')
lines(ecdf(foo[,"RARE_SAS_DEL.SV.bed.dist"]),col = 'purple')
lines(ecdf(shuffleChIPfoo[,"RARE_AFR_DEL.SV.bed.shuffleChIP.dist"]),col = 'grey', lty =2)
lines(ecdf(shuffleChIPfoo[,"RARE_AMR_DEL.SV.bed.shuffleChIP.dist"]),col = 'pink', lty =2)
lines(ecdf(shuffleChIPfoo[,"RARE_EUR_DEL.SV.bed.shuffleChIP.dist"]),col = 'lightblue', lty =2)
lines(ecdf(shuffleChIPfoo[,"RARE_EAS_DEL.SV.bed.shuffleChIP.dist"]),col = 'lightgreen', lty =2)
lines(ecdf(shuffleChIPfoo[,"RARE_SAS_DEL.SV.bed.shuffleChIP.dist"]),col = 'magenta', lty =2)
###########################################################################
###########################################################################




