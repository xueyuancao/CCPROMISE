PrbCor<-function(geneSet,     
               Edat,          
               Mdat,          
               EMlbl=c("Expr", "Methyl"),   
               phdat,         
               pcut=0.05) 
{
  res<-NULL
  gen<-NULL
  sbjnam<-dimnames(phdat)[[1]]
  for (i in 1:length(geneSet))
  { 
    thisgene<-geneSet[[i]]
    thisPrb<-geneIds(thisgene)
    geneName<-setName(thisgene)
    if (is.null(Mdat)) thisMthyl2<-NULL
    if (!is.null(Mdat)) thisMthyl2<-datSel(Mdat, thisPrb, prbcol=FALSE)   
    if (is.null(Edat)) thisExpr2<-NULL
    if (!is.null(Edat))   thisExpr2<-datSel(Edat, thisPrb, prbcol=FALSE)   
    
    if (is.null(thisMthyl2))
    {
      nrE<-nrow(thisExpr2)
      ncE<-ncol(thisExpr2)
      rnmE<-rownames(thisExpr2)
      this.res<-cbind(rep(geneName, nrE), rnmE, rep(NA, nrE), rep(NA, nrE), rep(NA, nrE))
      rownames(this.res)<-rnmE
      colnames(this.res)<-c('Gene', EMlbl, 'Spearman.rstat', 'Spearman.p')
      this.gen<-cbind(rnmE, thisExpr2, matrix(NA, nrow=nrE, ncol=ncE))
      dimnames(this.gen)[[2]]<-c("ProbPair", rep(sbjnam, 2))
    }
    
    if (is.null(thisExpr2))
    {
      nrM<-nrow(thisMthyl2)
      ncM<-ncol(thisMthyl2)
      rnmM<-rownames(thisMthyl2)    
      this.res<-cbind(rep(geneName, nrM), rep(NA, nrM),rnmM,  rep(NA, nrM), rep(NA, nrM))
      rownames(this.res)<-rnmM
      colnames(this.res)<-c('Gene', EMlbl, 'Spearman.rstat', 'Spearman.p')
      this.gen<-cbind(rnmM, matrix(NA, nrow=nrM, ncol=ncM), thisMthyl2)
      dimnames(this.gen)[[2]]<-c("ProbPair", rep(sbjnam, 2))
    }
    
    if (!is.null(thisExpr2) & !is.null(thisMthyl2))
    {
      nrE<-nrow(thisExpr2)
      ncE<-ncol(thisExpr2)
      rnmE<-rownames(thisExpr2)
          
      nrM<-nrow(thisMthyl2)
      ncM<-ncol(thisMthyl2)
      rnmM<-rownames(thisMthyl2)
      
      this.res<-NULL
      this.gen<-NULL
      for (k in 1:nrE)
      {
         spres<-t(apply(thisMthyl2, 1, function(prb, eprb){
                      tt<-cor.test(prb, as.numeric(eprb), method="spearman")
                      return(c(round(tt$estimate, 4), round(tt$p.value, 10)))
                      }, thisExpr2[k, ]))
         tt.res<-cbind(rep(geneName, nrM), rep(rnmE[k], nrM), rownames(spres), spres)
         rownames(tt.res)<-paste(rnmE[k],rownames(spres), sep='*')
         
         # sign of methylation signal is changed to reflect the correlation with expression. 
         msign<-apply(thisMthyl2, 2, function(prb, corsign){
                       return(prb*corsign)}, sign(spres[, 1]))
         if (is.vector(msign)) msign<-matrix(msign, nrow=1, ncol=length(msign))
         thisEk<-as.numeric(thisExpr2[k, ])
         tt.gen<-cbind(matrix(rep(thisEk, nrM), nrow=nrM, ncol=ncE, byrow=TRUE), msign)
         rownames(tt.gen)<-paste(rnmE[k],rownames(spres), sep='*')
         tt.gen2<-cbind(rownames(tt.gen), tt.gen) 
         this.res<-rbind(this.res, tt.res)
         this.gen<-rbind(this.gen, tt.gen2)
      }
      colnames(this.res)<-c('Gene', EMlbl, 'Spearman.rstat', 'Spearman.p')
      dimnames(this.gen)[[2]]<-c("ProbPair", rep(sbjnam, 2))
    }
    res<- rbind(res, this.res)
    gen<- rbind(gen, this.gen[!is.element(this.gen[,2], gen[,1]),])
  }
  dimnames(gen)[[2]]<-c("ProbPair", rep(sbjnam, 2))
  
  keep<-!is.na(res[, 'Spearman.p']) & as.numeric(res[, 'Spearman.p']) < pcut
  genf<-gen[keep, ]
  genf<-genf[!duplicated(genf[, 'ProbPair']), ]
  resf<-res[keep, ]  
  return(list(res=resf, gen=genf))
}