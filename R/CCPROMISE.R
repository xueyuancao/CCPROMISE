CCPROMISE<-function(geneSet,                    
                  ESet,                         
                  MSet,                         
                  promise.pattern,              
                  strat.var=NULL,               
                  prlbl=NULL,                   
                  EMlbl=c('Expr', 'Mthyl'),     
                  nbperm=FALSE,                 
                  max.ntail=100,                
                  nperms=10000,                 
                  seed=13)                      
{
  Edat<-exprs(ESet)
  Mdat<-exprs(MSet)
  pheDat<-phenoData(ESet)

  phdat<-pData(pheDat)
  Edat<-Edat[, is.element(dimnames(Edat)[[2]], dimnames(Mdat)[[2]]) & 
               is.element(dimnames(Edat)[[2]], dimnames(phdat)[[1]])]
  Mdat<-Mdat[, is.element(dimnames(Mdat)[[2]], dimnames(Edat)[[2]]) & 
               is.element(dimnames(Mdat)[[2]], dimnames(phdat)[[1]])]
  phdat<-phdat[is.element(dimnames(phdat)[[1]], dimnames(Edat)[[2]])& 
               is.element(dimnames(phdat)[[1]], dimnames(Mdat)[[2]]), ]
  
  #data checking
  if(any(dimnames(Edat)[[2]]!=dimnames(Mdat)[[2]])) 
  {
    message('Data Error!!!!!!!!!!!!!!!!!!!!!!!!!')
    break()
  }  
  if(any(dimnames(Edat)[[2]]!=dimnames(phdat)[[1]] ))
   {
    message('Data Error!!!!!!!!!!!!!!!!!!!!!!!!!')
    break()
  } 
  
  canres<-CANN(geneSet=geneSet, Edat=Edat, Mdat=Mdat, EMlbl=EMlbl, phdat=phdat)
  canndat<-canres$FSTccscore
  canndatf<-apply(canndat[, -1], 2, as.numeric)
  dimnames(canndatf)[[1]]<-canndat[,1]
  
  expr<-as.matrix(canndatf[, 1:(ncol(canndatf)/2)])
  mtyl<-as.matrix(canndatf[, (ncol(canndatf)/2+1):ncol(canndatf)])
  
  exprsetx<-new("ExpressionSet", exprs=expr, phenoData = pheDat) 
  exprsety<-new("ExpressionSet", exprs=mtyl, phenoData = pheDat)      
  
  PRres<-PROMISE2(exprSet=exprsetx, exprSet2=exprsety,
                    geneSet=NULL, 
                    promise.pattern=promise.pattern, 
                    strat.var=strat.var,  
                    nperms=nperms,
                    nbperm=nbperm,
                    max.ntail=max.ntail,
                    seed=seed)$generes
  if (is.null(prlbl ))  
  {prlbl<-c(promise.pattern[, 'endpt.vars'], paste('PR', nrow(promise.pattern), sep=""))}
  if (!nbperm)
  dimnames(PRres)[[2]]<-c('Gene', paste(paste(EMlbl[1], prlbl, sep="_"), 'Stat', sep="."), 
                            paste(paste(EMlbl[2], prlbl, sep="_"), 'Stat', sep="."), 
                            paste(prlbl[length(prlbl)], 'Stat', sep="."),
                            paste(paste(EMlbl[1], prlbl, sep="_"), 'Pval', sep="."), 
                            paste(paste(EMlbl[2], prlbl, sep="_"), 'Pval', sep="."), 
                            paste(prlbl[length(prlbl)], 'Pval', sep="."))
  if (nbperm) 
  dimnames(PRres)[[2]]<-c('Gene', paste(paste(EMlbl[1], prlbl, sep="_"), 'Stat', sep="."), 
                            paste(paste(EMlbl[2], prlbl, sep="_"), 'Stat', sep="."), 
                            paste(prlbl[length(prlbl)], 'Stat', sep="."),
                            paste(paste(EMlbl[1], prlbl, sep="_"), 'Pval', sep="."), 
                            paste(paste(EMlbl[2], prlbl, sep="_"), 'Pval', sep="."), 
                            paste(prlbl[length(prlbl)], 'Pval', sep="."), 'nperm')                                                              
  return(list(PRres=PRres, CCres=canres$CCres, 
              FSTccscore=canres$FSTccscore, CCload=canres$CCload))
} 