colnames(z)
#[1] "row.id"                 "chr"                    "start"                  "stop"                   "mPGC.WK7.me"           
#[6] "mPGC.WK10.me"           "mPGC.WK11.me"           "mPGC.WK13.me"           "mPGC.WK19.me"           "sperm_CpGMe"           
#[11] "Rearrangement_Hum_Spec" "Phase3_all_bkpts"       "Phase3_mCNV"            "Phase3_DUP"             "Phase3_DEL"            
#[16] "Phase3_INV"             "CNV_Bipolar_Case"       "CNV_Bipolar_Control"    "CNV_Autism_Case"        "CNV_Autism_Control"    
#[21] "CNV_DevDelay_Case"      "CNV_DevDelay_Control"   "CNV_Schiz_Case"         "CNV_Schiz_Control"      "CNV_WTCCC"             
#[26] "CNV_270_HapMap"         "CNV_450_HapMap"         "CNV_400_MGL"            "Cmyc.ENCODE"            "Ctcf.ENCODE"           
#[31] "Foxa1.ENCODE"           "Jund.ENCODE"            "Mef2a.ENCODE"           "Mef2c.ENCODE"           "Pou5f1.ENCODE"         
#[36] "Prdm1.ENCODE"           "PRDM14_ES_PEAKS"        "Rad21.ENCODE"           "Suz12.ENCODE"           "Yy1.ENCODE"            
#[41] "Znf143.ENCODE"          "Znf263.ENCODE"          "Znf274.ENCODE"          "Havana.gene.gencode" 
mrx1 = matrix(data = rep(NA,((44-29+1)^2)),nrow = (44-29+1),dimnames = list(colnames(z)[29:44],colnames(z)[29:44]))
for(i in (29:44)){
    for(j in (29:44)){
        mrx1[(i-28),(j-28)] = getTileEnrichment.pres(featureTilesDF = z,featurePair = c(i,j),mode =  "pct")
        
    }
}
dev.new()

heatmap.2(mrx1,na.rm = T,trace = 'none',margins = c(10,12),
          col = bluered(40),Colv = T,Rowv = F,main = "TF-TF Enrichment in 100 kb tiles" )


mrx2 = matrix(data = rep(NA,((44-29+1)*(28-12+1))),nrow = (28-12+1),dimnames = list(colnames(z)[12:28],colnames(z)[29:44]))
for(i in (12:28)){
    for(j in (29:44)){
        mrx2[(i-11),(j-28)] = getTileEnrichment.pres(featureTilesDF = z,featurePair = c(i,j),mode =  "pct")
        
    }
}
dev.new()

heatmap.2(mrx2,na.rm = T,trace = 'none',margins = c(10,12),
          col = bluered(40),Colv = T,Rowv = F,main = "SV-TF Enrichment in 100 kb tiles" )


vv = getTileEnrichment.pres(featureTilesDF = z,featurePair = c(44,11),mode =  "enr")
vv
