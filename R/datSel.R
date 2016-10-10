datSel<-function(dat,           
                 thisPrb,       
                 prbcol=TRUE)   
{
  ind<-is.element(dimnames(dat)[[1]], thisPrb)
  if (sum(ind)==0) return(NULL)
  if (prbcol)
  {
    if (sum(ind)==1)
    {
       dat2<-data.frame(matrix(dat[ind, ], nrow=ncol(dat), ncol=1))
       dimnames(dat2)[[2]]<-dimnames(dat)[[1]][ind]
    }
    if (sum(ind)>1)
    {
       dat2<-t(dat[ind, ])
    }
  }
  if (!prbcol)
  {
    if (sum(ind)==1)
    {
       dat2<-data.frame(matrix(dat[ind, ],  nrow=1, ncol=ncol(dat)))
       dimnames(dat2)[[1]]<-dimnames(dat)[[1]][ind]
    }
    if (sum(ind)>1)
    {
       dat2<-dat[ind, ]
    }
  }  
  return(dat2)
}