

plot(ecdf(foo[,"RARE_AFR_CNV.SV.bed.dist"]),col = 'black',xlim = c(0,5000000),)
lines(ecdf(foo[,"RARE_AMR_CNV.SV.bed.dist"]),col = 'red')
lines(ecdf(foo[,"RARE_EUR_CNV.SV.bed.dist"]),col = 'darkblue')
lines(ecdf(foo[,"RARE_EAS_CNV.SV.bed.dist"]),col = 'darkgreen')
lines(ecdf(foo[,"RARE_SAS_CNV.SV.bed.dist"]),col = 'purple')
lines(ecdf(shuffleSVfoo[,"shuff.RARE_AFR_CNV.SV.bed.shuffleSVs.dist"]),col = 'grey', lty =2)
lines(ecdf(shuffleSVfoo[,"shuff.RARE_AMR_CNV.SV.bed.shuffleSVs.dist"]),col = 'pink', lty =2)
lines(ecdf(shuffleSVfoo[,"shuff.RARE_EUR_CNV.SV.bed.shuffleSVs.dist"]),col = 'lightblue', lty =2)
lines(ecdf(shuffleSVfoo[,"shuff.RARE_EAS_CNV.SV.bed.shuffleSVs.dist"]),col = 'lightgreen', lty =2)
lines(ecdf(shuffleSVfoo[,"shuff.RARE_SAS_CNV.SV.bed.shuffleSVs.dist"]),col = 'magenta', lty =2)


plot(ecdf(foo[,"RARE_AFR_INV.SV.bed.dist"]),col = 'black')
lines(ecdf(foo[,"RARE_AMR_INV.SV.bed.dist"]),col = 'red')
lines(ecdf(foo[,"RARE_EUR_INV.SV.bed.dist"]),col = 'blue')
lines(ecdf(foo[,"RARE_EAS_INV.SV.bed.dist"]),col = 'green')
lines(ecdf(foo[,"RARE_SAS_INV.SV.bed.dist"]),col = 'purple')
lines(ecdf(shuffleSVfoo[,"shuff.RARE_AFR_INV.SV.bed.shuffleSVs.dist"]),col = 'black', lty =2)
lines(ecdf(shuffleSVfoo[,"shuff.RARE_AMR_INV.SV.bed.shuffleSVs.dist"]),col = 'red', lty =2)
lines(ecdf(shuffleSVfoo[,"shuff.RARE_EUR_INV.SV.bed.shuffleSVs.dist"]),col = 'blue', lty =2)
lines(ecdf(shuffleSVfoo[,"shuff.RARE_EAS_INV.SV.bed.shuffleSVs.dist"]),col = 'green', lty =2)
lines(ecdf(shuffleSVfoo[,"shuff.RARE_SAS_INV.SV.bed.shuffleSVs.dist"]),col = 'purple', lty =2)

plot(ecdf(foo[,"RARE_AFR_DUP.SV.bed.dist"]),col = 'black')
lines(ecdf(foo[,"RARE_AMR_DUP.SV.bed.dist"]),col = 'red')
lines(ecdf(foo[,"RARE_EUR_DUP.SV.bed.dist"]),col = 'blue')
lines(ecdf(foo[,"RARE_EAS_DUP.SV.bed.dist"]),col = 'green')
lines(ecdf(foo[,"RARE_SAS_DUP.SV.bed.dist"]),col = 'purple')
lines(ecdf(shuffleSVfoo[,"shuff.RARE_AFR_DUP.SV.bed.shuffleSVs.dist"]),col = 'black', lty =2)
lines(ecdf(shuffleSVfoo[,"shuff.RARE_AMR_DUP.SV.bed.shuffleSVs.dist"]),col = 'red', lty =2)
lines(ecdf(shuffleSVfoo[,"shuff.RARE_EUR_DUP.SV.bed.shuffleSVs.dist"]),col = 'blue', lty =2)
lines(ecdf(shuffleSVfoo[,"shuff.RARE_EAS_DUP.SV.bed.shuffleSVs.dist"]),col = 'green', lty =2)
lines(ecdf(shuffleSVfoo[,"shuff.RARE_SAS_DUP.SV.bed.shuffleSVs.dist"]),col = 'purple', lty =2)

plot(ecdf(foo[,"RARE_AFR_DEL.SV.bed.dist"]),col = 'black')
lines(ecdf(foo[,"RARE_AMR_DEL.SV.bed.dist"]),col = 'red')
lines(ecdf(foo[,"RARE_EUR_DEL.SV.bed.dist"]),col = 'blue')
lines(ecdf(foo[,"RARE_EAS_DEL.SV.bed.dist"]),col = 'green')
lines(ecdf(foo[,"RARE_SAS_DEL.SV.bed.dist"]),col = 'purple')
lines(ecdf(shuffleSVfoo[,"shuff.RARE_AFR_DEL.SV.bed.shuffleSVs.dist"]),col = 'black', lty =2)
lines(ecdf(shuffleSVfoo[,"shuff.RARE_AMR_DEL.SV.bed.shuffleSVs.dist"]),col = 'red', lty =2)
lines(ecdf(shuffleSVfoo[,"shuff.RARE_EUR_DEL.SV.bed.shuffleSVs.dist"]),col = 'blue', lty =2)
lines(ecdf(shuffleSVfoo[,"shuff.RARE_EAS_DEL.SV.bed.shuffleSVs.dist"]),col = 'green', lty =2)
lines(ecdf(shuffleSVfoo[,"shuff.RARE_SAS_DEL.SV.bed.shuffleSVs.dist"]),col = 'purple', lty =2)
###########################################################################
###########################################################################




