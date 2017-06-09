plot(ecdf(foo[,"RARE_AFR_CNV.SV.bed.dist"]),col = 'black')
lines(ecdf(foo[,"RARE_AMR_CNV.SV.bed.dist"]),col = 'red')
lines(ecdf(foo[,"RARE_EUR_CNV.SV.bed.dist"]),col = 'blue')
lines(ecdf(foo[,"RARE_EAS_CNV.SV.bed.dist"]),col = 'green')
lines(ecdf(foo[,"RARE_SAS_CNV.SV.bed.dist"]),col = 'purple')

plot(ecdf(foo[,"RARE_AFR_INV.SV.bed.dist"]),col = 'black')
lines(ecdf(foo[,"RARE_AMR_INV.SV.bed.dist"]),col = 'red')
lines(ecdf(foo[,"RARE_EUR_INV.SV.bed.dist"]),col = 'blue')
lines(ecdf(foo[,"RARE_EAS_INV.SV.bed.dist"]),col = 'green')
lines(ecdf(foo[,"RARE_SAS_INV.SV.bed.dist"]),col = 'purple')

plot(ecdf(foo[,"RARE_AFR_DUP.SV.bed.dist"]),col = 'black')
lines(ecdf(foo[,"RARE_AMR_DUP.SV.bed.dist"]),col = 'red')
lines(ecdf(foo[,"RARE_EUR_DUP.SV.bed.dist"]),col = 'blue')
lines(ecdf(foo[,"RARE_EAS_DUP.SV.bed.dist"]),col = 'green')
lines(ecdf(foo[,"RARE_SAS_DUP.SV.bed.dist"]),col = 'purple')

plot(ecdf(foo[,"RARE_AFR_DEL.SV.bed.dist"]),col = 'black')
lines(ecdf(foo[,"RARE_AMR_DEL.SV.bed.dist"]),col = 'red')
lines(ecdf(foo[,"RARE_EUR_DEL.SV.bed.dist"]),col = 'blue')
lines(ecdf(foo[,"RARE_EAS_DEL.SV.bed.dist"]),col = 'green')
lines(ecdf(foo[,"RARE_SAS_DEL.SV.bed.dist"]),col = 'purple')
###########################################################################
###########################################################################
plot(ecdf(shuffleSVfoo[,"shuff.RARE_AFR_CNV.SV.bed.shuffleSVs.dist"]),col = 'black')
lines(ecdf(shuffleSVfoo[,"shuff.RARE_AMR_CNV.SV.bed.shuffleSVs.dist"]),col = 'red')
lines(ecdf(shuffleSVfoo[,"shuff.RARE_EUR_CNV.SV.bed.shuffleSVs.dist"]),col = 'blue')
lines(ecdf(shuffleSVfoo[,"shuff.RARE_EAS_CNV.SV.bed.shuffleSVs.dist"]),col = 'green')
lines(ecdf(shuffleSVfoo[,"shuff.RARE_SAS_CNV.SV.bed.shuffleSVs.dist"]),col = 'purple')

plot(ecdf(shuffleSVfoo[,"shuff.RARE_AFR_INV.SV.bed.shuffleSVs.dist"]),col = 'black')
lines(ecdf(shuffleSVfoo[,"shuff.RARE_AMR_INV.SV.bed.shuffleSVs.dist"]),col = 'red')
lines(ecdf(shuffleSVfoo[,"shuff.RARE_EUR_INV.SV.bed.shuffleSVs.dist"]),col = 'blue')
lines(ecdf(shuffleSVfoo[,"shuff.RARE_EAS_INV.SV.bed.shuffleSVs.dist"]),col = 'green')
lines(ecdf(shuffleSVfoo[,"shuff.RARE_SAS_INV.SV.bed.shuffleSVs.dist"]),col = 'purple')

plot(ecdf(shuffleSVfoo[,"shuff.RARE_AFR_DUP.SV.bed.shuffleSVs.dist"]),col = 'black')
lines(ecdf(shuffleSVfoo[,"shuff.RARE_AMR_DUP.SV.bed.shuffleSVs.dist"]),col = 'red')
lines(ecdf(shuffleSVfoo[,"shuff.RARE_EUR_DUP.SV.bed.shuffleSVs.dist"]),col = 'blue')
lines(ecdf(shuffleSVfoo[,"shuff.RARE_EAS_DUP.SV.bed.shuffleSVs.dist"]),col = 'green')
lines(ecdf(shuffleSVfoo[,"shuff.RARE_SAS_DUP.SV.bed.shuffleSVs.dist"]),col = 'purple')

plot(ecdf(shuffleSVfoo[,"shuff.RARE_AFR_DEL.SV.bed.shuffleSVs.dist"]),col = 'black')
lines(ecdf(shuffleSVfoo[,"shuff.RARE_AMR_DEL.SV.bed.shuffleSVs.dist"]),col = 'red')
lines(ecdf(shuffleSVfoo[,"shuff.RARE_EUR_DEL.SV.bed.shuffleSVs.dist"]),col = 'blue')
lines(ecdf(shuffleSVfoo[,"shuff.RARE_EAS_DEL.SV.bed.shuffleSVs.dist"]),col = 'green')
lines(ecdf(shuffleSVfoo[,"shuff.RARE_SAS_DEL.SV.bed.shuffleSVs.dist"]),col = 'purple')

