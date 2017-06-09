######  Rename this or generalize it to allow enrichment by
######  modes: (percent thresholds v count, count v count, pres abs v count
######                    pres abs v pres abs)
######  
######  getTileEnrichment.R
######  By: Rocco Lucero
######  Created: February 2, 2016
######  Last Modified February 10, 2016 by Rocco Lucero
######
######  getTileEnrichment. R is a script that takes as input a data frame of
######  genome feature  counts (or levels) in uniformly-sized genomic
######  ranges (tiles) and computes the enrichment of the overlap
######  between a pair of the available features within the set of tiles
######  defined by a fixed percentage supplied by the user. Fold enrichment
######  is reported.
######
######  FOI IS A VECTOR OF COL INDICES OR COLNAMES

getTileEnrichment.slice = function(featureTilesDF,foi,hi.pct = 1,lo.pct = 0,precision = 100){
    set.seed(20160213)
    #featureTilesDF =z
    #foi = c(10,13)
    #colnames(z)[foi]
    df = featureTilesDF[order(featureTilesDF[,foi[1]]),c(1,foi)]
    df = df[complete.cases(df),]
    print(head(df,30))
    #head(df)
    
    comp.slice.Enrichment = function(df=df,hi.pct = hi.pct,lo.pct = lo.pct){
        totalTiles = nrow(df)
        hi.testTiles = round(hi.pct/100*totalTiles,digits = 0)
        if(lo.pct != 0){
            lo.testTiles = round(lo.pct/100*totalTiles,digits = 0)
            hi.df = df[1:hi.testTiles,]
            lo.df = df[1:lo.testTiles,]
            slice.df = hi.df[!((hi.df$row.id) %in% (lo.df$row.id)),]
            print(rownames(slice.df))
            x = (df$row.id %in% slice.df$row.id)
            other.df = df[!x,]
        }
        else{
            hi.df = df[1:hi.testTiles,]
            slice.df = hi.df
            x = (df$row.id %in% slice.df$row.id)
            other.df = df[!x,]
        }
        
        slice.Enrich = sum(slice.df[,3])/hi.testTiles 
        other.Enrich = sum(other.df[,3])/(totalTiles - hi.testTiles) 
        Target.enrichment = slice.Enrich/other.Enrich
        
    }
    
    perm.P.Val = function(df=df,hi.pct = hi.pct,lo.pct = lo.pct,trgEnr,prec = precision){
        #sample "testTiles" number of tiles randomly 1000 times
        #and compute enrichmnet p value
        pct = hi.pct - lo.pct
        perm.dist = c(rep(NULL,prec))
        for(i in 1:prec){
            testTiles = round(pct/100*nrow(df),digits = 0)
            my.samp = sample(x = rownames(df),size = testTiles,replace = F)
            df = df[my.samp,]
            sampledEnr = comp.slice.Enrichment(df=df,hi.pct = hi.pct,lo.pct = lo.pct)
            perm.dist[i] = sampledEnr
        }
        ####
        if(trgEnr>=1){length(perm.dist[which(perm.dist>trgEnr)])/prec}
        else{length(perm.dist[which(perm.dist<trgEnr)])/prec}
    }
    
    enr = comp.slice.Enrichment(df=df,hi.pct = hi.pct,lo.pct = lo.pct)
    p.val = perm.P.Val(df = df,hi.pct = hi.pct,lo.pct = lo.pct,trgEnr = enr,prec = precision)
    if (p.val == 0){p.val = paste("<",(1/precision),sep = " ")}
    result = list(paste(colnames(df)[1],"vs.",colnames(df)[2],sep=" "),enr,p.val)
    names(result) = c("comparison","fold enrichment","p-value")
    result

}
mm = getTileEnrichment.slice(featureTilesDF = z,foi = c(10,13),hi.pct = 10,lo.pct = 5)
mm
nn = getTileEnrichment(featureTilesDF = z,foi = c(10,13),pct = 5,rev = T)
nn

