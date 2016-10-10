CANN<-function(geneSet,       
               Edat,          
               Mdat,          
               EMlbl=c("Expr", "Methyl"),      
               phdat)         
{
  res<-NULL
  gen<-NULL
  candtl<-NULL
  sbjnam<-dimnames(phdat)[[1]]
  for (i in 1:length(geneSet))
  { 
    thisgene<-geneSet[[i]]
    thisPrb<-geneIds(thisgene)
    geneName<-setName(thisgene)
    if (is.null(Mdat)) thisMthyl2<-NULL
    if (!is.null(Mdat)) thisMthyl2<-datSel(Mdat, thisPrb)   
    if (is.null(Edat)) thisExpr2<-NULL
    if (!is.null(Edat))   thisExpr2<-datSel(Edat, thisPrb)
    
    
    if (is.null(thisMthyl2))
    {
      nrE<-nrow(thisExpr2)
      ncE<-ncol(thisExpr2)
      x=matrix(prcomp(thisExpr2, retx=TRUE, scale=TRUE)$x[,1], nrow=nrE, ncol=1)
      y<-NULL
      this.res<-matrix(c(geneName, ncE, 0, rep(NA, 3)), nrow=1, ncol=6) 
      this.gen<-c(geneName, as.vector(x), rep(NA, nrow(x)))
      this.candtl<-NULL    
    }
    
    if (is.null(thisExpr2))
    {
      nrM<-nrow(thisMthyl2)
      ncM<-ncol(thisMthyl2) 
      y=matrix(prcomp(thisMthyl2, retx=TRUE, scale=TRUE)$x[,1], nrow=nrM, ncol=1)
      x<-NULL
      this.res<-matrix(c(geneName, 0, ncM, rep(NA, 3)), nrow=1, ncol=6) 
      this.gen<-c(geneName, rep(NA, nrow(y)), as.vector(y))
      this.candtl<-NULL 
    }
    
    if (!is.null(thisExpr2) & !is.null(thisMthyl2))
    {
      nrE<-nrow(thisExpr2)
      ncE<-ncol(thisExpr2)
      
      nrM<-nrow(thisMthyl2)
      ncM<-ncol(thisMthyl2)
           
      this.candtl<-NULL
      cc.res <-cancor(thisExpr2, thisMthyl2) 
      Mcnt<-thisMthyl2
      for(j in 1:ncol(Mcnt))  Mcnt[, j]<-Mcnt[, j] - cc.res$ycenter[j]
      Ecnt<-thisExpr2
      for(j in 1:ncol(Ecnt)) Ecnt[, j]<-Ecnt[, j] - cc.res$xcenter[j] 
      cr<-cc.res$cor[1]
      x=as.matrix(Ecnt[, 1:nrow(cc.res$xcoef)])%*%cc.res$xcoef[,1]*sign(cr)
      y=as.matrix(Mcnt[, 1:nrow(cc.res$ycoef)])%*%cc.res$ycoef[,1]
      this.candtl<-c(geneName, paste(names(cc.res$xcoef[,1]), collapse="|"),  
                               paste(cc.res$xcoef[,1], collapse="|"), 
                               paste(names(cc.res$ycoef[,1]), collapse="|"),  
                               paste(cc.res$ycoef[,1], collapse="|"))
      capture.output(ccperm<-p.perm(thisExpr2, thisMthyl2, 
                                    nboot = 999, rhostart = 1, type = "Wilks")$p, file='NUL')
      capture.output(ccasym<-p.asym(cc.res$cor, nrE, ncE, ncM, 
                                    tstat = "Wilks")$p[1], file='NUL')
      
      this.res<-matrix(c(geneName, ncE, ncM, cr, ccperm, ccasym), nrow=1, ncol=6)
      this.gen<-c(geneName, as.vector(x), as.vector(y))    
    }
    dimnames(this.res)[[2]]<-c('Gene', paste('n.', EMlbl, sep=""), 
                               'CanonicalCR', 'WilksPermPval', 'WilksAsymPval')
    res<- rbind(res, this.res)
    gen<- rbind(gen, this.gen)
    candtl<-rbind(candtl, this.candtl)
  }
  dimnames(gen)[[2]]<-c("Gene", rep(sbjnam, 2))
  dimnames(gen)[[1]]<-gen[,1]
  if (!is.null(candtl)) 
  {
    dimnames(candtl)[[2]]<-c("Gene", 'Eprb', 'Eload', 'Mprb', 'Mload')
    dimnames(candtl)[[1]]<-candtl[,1]
  }
  dimnames(res)[[1]]<-res[,1]
  return(list(CCres=res, FSTccscore=gen, CCload=candtl))
}