By_Population_SV_PRDM14_Hits <- read.csv("C:/Users/Rocco/Desktop/By_Population_SV_PRDM14_Hits.csv", stringsAsFactors=FALSE)
By_Population_SV_PRDM14_Hits_1 = By_Population_SV_PRDM14_Hits[,c(1,2,3,5,6)]
View(By_Population_SV_PRDM14_Hits_1)

#a = as.matrix(By_Population_SV_PRDM14_Hits_1[,2:5],byrow = T)
a = as.matrix(By_Population_SV_PRDM14_Hits_1[,c(4,5,2,3)],byrow = T)

for (n in 1:nrow(a)){
  b = fisher.test(x = matrix(a[n,],nrow = 2,byrow = T))
  print((By_Population_SV_PRDM14_Hits_1[n,1]))
  print(b$estimate)
  print(b$p.value)
  print(b$conf.int)
  cat("\n")
}

# (RareHits/RareTotal) / (CommonHits/CommonTotal)
c = 1/((By_Population_SV_PRDM14_Hits[c(2,6,10,14,18),2]/By_Population_SV_PRDM14_Hits[c(2,6,10,14,18),4])/(By_Population_SV_PRDM14_Hits[c(2,6,10,14,18),5]/By_Population_SV_PRDM14_Hits[c(2,6,10,14,18),7]))
c
d = 1/((By_Population_SV_PRDM14_Hits[,2]/By_Population_SV_PRDM14_Hits[,4])/(By_Population_SV_PRDM14_Hits[,5]/By_Population_SV_PRDM14_Hits[,7]))
d
class(a[1,])
class(a)
b$p.value

BySV_PRDM14_Hits <- read.csv("C:/Users/Rocco/Desktop/BySV_PRDM14_Hits.csv")
BySV_PRDM14_Hits_1 = BySV_PRDM14_Hits[,c(1,2,3,5,6)]
a2 = as.matrix(BySV_PRDM14_Hits_1[,c(4,5,2,3)],byrow = T)

for (n in 1:nrow(a2)){
  b2 = fisher.test(x = matrix(a2[n,],nrow = 2,byrow = T))
  print((BySV_PRDM14_Hits_1[n,1]))
  print(b2$estimate)
  print(b2$p.value)
  print(b2$conf.int)
  cat("\n")
}
a2
b2

