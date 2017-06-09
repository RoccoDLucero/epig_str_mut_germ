##Reproduce the Analysis from Harris et.al. PLOS 2013 but include the PRDM14 ES-Cell
# ChIP-seq peak counts from Chia et. al.
#The structural variants are described in Jian Li et.al 2012.
#Here the genome is binned in 100kb Windows for the counts:
Confounders_Data <- read.delim("./input/Harris_Confounders_Data/Table_S1/Confounders_Data.txt")
hg18_PRDM14_vs_Molaro_nonMD1 <- read.delim("./input/hg18_PRDM14_vs_Molaro_nonMD1.txt", header=FALSE)
hg18_PRDM14_vs_Molaro_Lowest1 <- read.delim("./input/hg18_PRDM14_vs_Molaro_Lowest1.txt", header=FALSE)
PRDM14ES_PeaksPerWindow = rbind(hg18_PRDM14_vs_Molaro_Lowest1[,c(1:3,7)],hg18_PRDM14_vs_Molaro_nonMD1[,c(1:3,7)])

colnames(PRDM14ES_PeaksPerWindow) = c(colnames(Confounders_Data)[3:5],'PRDM14_ES_PEAKS')
rm(hg18_PRDM14_vs_Molaro_Lowest1,hg18_PRDM14_vs_Molaro_nonMD1)
#Later I need to ensure that merge assigns values properly to each row...
Confounders_Data = merge.data.frame(x = Confounders_Data[,3:17],y = PRDM14ES_PeaksPerWindow,by = c('Chromosome','Start','End'))
Confounders_Data = Confounders_Data[,colnames(Confounders_Data)[c(1:10,16,11:15)]]

#Define all of the leave one out formulae by generating complete and leave one out combinations of factors for each of the response variables:
#Example Formula: (CNV_270HapMap)~(CpG_Island_Coverage)+(LINE_Coverage)+(LTR_Coverage)+(Satellite_Coverage)+(SINE_Coverage)+(Methylation)

#Stored analyses of the form 'Response.Desctiption.Linkfunction."glm"

LOODescriptions = c('complete','noMeth','noCpG','noLINE','noLTR','noSat','noSINE','noPRDM14')
FactorVars = colnames(Confounders_Data)[c(4:5,7:11)]
ResponseVars = colnames(Confounders_Data)[12:16]
LinkFuncs = c("negBin", "pois","linear")             

#Generate Variable Names for output
Model_OutputVars = c()
for (response in ResponseVars){
  for (func in LinkFuncs){
    for (description in LOODescriptions){
      my_model = paste(response,'.',description,'.',func,'.glm',sep = '')
      Model_OutputVars[length(Model_OutputVars)+1] = my_model
    }
  }
}               
Model_OutputVars = matrix(data = Model_OutputVars, nrow = length(LOODescriptions))


# Write all of the formula Tails without response variables
LOO_Combs =combn(x = rev(seq(1:(length(LOODescriptions)-1))),m = length(LOODescriptions)-2) #generate column vectors describing each of the LOO combinations
My_FormulaTails = c()
My_FormulaTails[length(My_FormulaTails)+1] = paste('~(',paste(FactorVars,collapse = ')+('),')',sep = '')
for(comb in 1:ncol(LOO_Combs)){
  TmpString = '~'
  for(i in LOO_Combs[,comb]){
    TmpString = paste(TmpString,'(',FactorVars[i],')+',sep = '')
  }
  FormulaTail =substr(TmpString,1,nchar(TmpString)-1)
  My_FormulaTails[length(My_FormulaTails)+1] = FormulaTail
}           

#Prepend Response variable into formula
Full_Formulas = c()
for(resp in ResponseVars){
  for(tail in My_FormulaTails){
    Full_Formulas[length(Full_Formulas)+1] = paste('(',resp,')',tail,sep = '')
  }  
}
Full_Formulas = matrix(data = Full_Formulas, nrow = length(LOODescriptions))
Model_OutputVars = Model_OutputVars[,c(seq(1,length(LinkFuncs)*length(ResponseVars),
                                           length(LinkFuncs)),seq(2,length(LinkFuncs)*length(ResponseVars),length(LinkFuncs)),
                                       seq(3,length(LinkFuncs)*length(ResponseVars),length(LinkFuncs)))]
#Full_Formulas_Matched = cbind(Full_Formulas,Full_Formulas,Full_Formulas)
#Generate tuples (output storage variable, appropriate formula)

#These packages are required for the Regression analysis:
require(foreign)
library(reshape2)
require(ggplot2)
require(MASS) # Contains the glm.nb functionality

#For Negative Binomaial
(rm(i,j,tail,resp,response,TmpString,func,comb,description))
rm(i,j,SV_Poisson,tmp)
SV_Pois = list()
SV_Neg_Bin = list()
SV_Lin = list()
for(j in 1:ncol(Full_Formulas)){
  for(i in 1:nrow(Full_Formulas)){
    #For Negative Binomaial
    SV_Neg_Bin[[Model_OutputVars[i,j]]] = glm.nb(formula = Full_Formulas[i,j], data =  Confounders_Data)
    #For Poisson
    #SV_Pois[[Model_OutputVars[i,j+5]]] = glm(formula = Full_Formulas[i,j], data =  Confounders_Data,family = "poisson")    
    #For Linear Model
    #SV_Lin[[Model_OutputVars[i,j+10]]] = glm(formula = Full_Formulas[i,j], data =  Confounders_Data,family ="linear")
  }
}

aicDiffs = c()
for(respGroup in seq(1,length(ResponseVars))){
  for(AIC in seq(1,length(LOODescriptions))){
    complete = names(SV_Neg_Bin)[(1+(respGroup-1)*length(LOODescriptions))]
    LOO = names(SV_Neg_Bin)[AIC+(respGroup-1)*length(LOODescriptions)]
    if((1+(respGroup-1)*length(LOODescriptions))!=AIC+(respGroup-1)*length(LOODescriptions)){
      improvement = SV_Neg_Bin[[LOO]]$aic - SV_Neg_Bin[[complete]]$aic
      aicDiffs = c(aicDiffs,improvement)
    }
  }
}
aicDiffs = as.data.frame(matrix(aicDiffs,byrow = T,nrow = length(ResponseVars)))
rownames(aicDiffs) = ResponseVars
colnames(aicDiffs) = LOODescriptions[2:length(LOODescriptions)]  
aicDiffs

barnames = c('methylation','CpG_Island','LINE','LTR','Satellite','SINE','PRDM14 Peaks')
groupnames = c('CNV\n270HapMap\nCount','CNV\n450HapMap\nCount',
               'CNV\nWTCCC\nCount','rare CNV\nSchizophrenia Case\nCount',
               'rare CNV\nSchizophrenia Control\nCount' )
MyColors = c("blue","red","white","light blue",'purple','yellow','green')
pdf(file = "./output/NegBin_Harris_Confounders_and_PRDM14pks.pdf",width = 9)
par(mar=c(6,6,6,6),xpd=T)
barplot(t(as.matrix(aicDiffs)), beside=TRUE,
        ylab = "Negative Binomial Model \n Improvement upon inclusion of factor",
        main = "GLM with Inclusion of PRDM14",
        col = MyColors, axes = F, ylim = c(-25,300), cex.axis = .1,
        cex.names = .6, names.arg = groupnames, xpd = F
)

legend(x= 38,y = 250, barnames,fill = MyColors,cex = .6)
axis(side = 2,at = seq(-20,280,by = 10))
dev.off()

