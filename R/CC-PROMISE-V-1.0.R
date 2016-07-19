
#> BEGIN DESCRIPTION
#> Package: CCPROMISE
#> Type: Package
#> Title: PROMISE analysis with Canonical Correlation for Two Forms of High Dimensional Genetic Data
#> Version: 1.0
#> Author: Xueyuan Cao <xueyuan.cao@stjude.org> and Stanley.pounds <stanley.pounds@stjude.org>
#> Maintainer: Xueyuan Cao <xueyuan.cao@stjude.org>
#> Description: Perform Canonical correlation between two forms of high demensional genetic data, and associate the first compoent of each form of data with a specific biologically interesting pattern of associations with multiple endpoints of variables. A probe level analysis is also implemented.
#> Depends: R (>= 2.15.0), CCP, PROMISE, Biobase, GSEABase
#> License: GPL (>= 2)
#> biocViews: Microarray, Bioinformatics, Gene expression
#> LazyLoad: yes
#> END DESCRIPTION

#> BEGIN NAMESPACE
#> exportPattern("^([^.]|\\..+\\.)")
#> export("CCPROMISE", "PrbPROMISE", "CANN")
#> importFrom(CCP, p.perm, p.asym)
#> importFrom(stats, cancor, cor.test, prcomp)
#> importFrom(PROMISE, promise.genestat)
#> importMethodsFrom(Biobase, exprs, pData, phenoData, ExpressionSet)
#> importFrom(GSEABase, setName,  geneIds)
#> END NAMESPACE

#> BEGIN CCPROMISE-package
#> \Rdversion{3.3.0}
#> \docType{package}
#> \title{
#>   PRojection Onto the Most Interesting Statistical Evidence with Canonical Correlation
#> }
#> \description{
#>   a tool to indentify genes that are correlated between two set of genomic variables and are associated with a predefined  
#>   pattern of associations with multiple endpoint variables.
#> }
#> \details{
#>   \tabular{ll}{
#>     Package: \tab CCPROMISE\cr
#>     Type: \tab Package\cr
#>     Version: \tab 1.0.0\cr
#>     Date: \tab 2015-6-8\cr
#>     License: \tab GPL (>=2)\cr
#>     LazyLoad: \tab yes\cr
#>   }
#>   The CCPROMISE (Canonical correlation  with PROMISE analysis) is performed by calling function CCPROMISE.
#>   The two forms of genomic data such as gene expression and methylation are passed through minimal ExpressionSet;
#>   the gene annotation (defining relationship between a gene and the two forms of genomic data), phenotypic data and definition of 
#    of pattern of association are passed through data frames; \emph{PROMISE2} and \emph{CANN} are called internally by CCPROMISE.
#>   R routines for calculating association statistics with individual endpoint variable are same as in \emph{PROMISE} package. 
#>   Please refer to \emph{PROMISE} package for writing user defined routines.
#> }
#> \author{Xueyuan Cao \email{Xueyuan.cao@stjude.org}, Stanley Pounds \email{stanley.pounds@stjude.org}
#>
#>      Maintainer: Xueyuan Cao \email{xueyuan.cao@stjude.org}
#> }
#> \references{
#>  Hotelling H. (1936). Relations between two sets of variables. Biometrika, 28, 321-327
#>  
#>  Pounds S, Cheng C, Cao X, Crews KR, Plunkett W, Gandhi V, Rubnitz J, Ribeiro RC, Downing JR, and Lamba J (2009) 
#>  PROMISE: a tool to identify genomic features with a specific biologically interesting pattern of associations 
#>  with multiple endpoint variables. Bioinformatics 25: 2013-2019
#>  
#>  Wilks, S. S. (1935) On the independence of k sets of normally distributed statistical variables. Econometrica, 3 309-326.
#> }
#> \keyword{ package }
#> \examples{
#> ## load data
#>  data(exmplESet)
#>  data(exmplMSet)
#>  data(exmplGeneSet)
#>  data(exmplPat)
#>  ## Perform CCPROMISE test
#> test<- CCPROMISE(geneSet=exmplGeneSet, 
#>              ESet=exmplESet, 
#>              MSet=exmplMSet, 
#>              promise.pattern=exmplPat,
#>              strat.var=NULL,
#>              prlbl=NULL, 
#>              EMlbl=c("Expr", "Methyl"),
#>              nbperm=TRUE,
#>              max.ntail=10,
#>              nperms=100,
#>              seed=13)    
#> }

#> END CCPROMISE-package

#> BEGIN exmplGeneSet
#> \title{Example of Gene Annotation}
#> \description{An exmple of gene set collection to annotate both form of genomic data to genes. The gene names can be extracted by method of setName() and probe ids can be extracted by method of geneIds()}
#> \usage{data(exmplGeneSet)}
#> \keyword{misc}
#> END exmplGeneSet

#> BEGIN exmplESet
#> \title{Example of Expression Set}
#> \description{an ExpressionSet class contains minimum of exprs (expression matrix) of gene expression and phenoData (AnnotatedDataFrame of end point data).}
#> \usage{data(exmplESet)}
#> \keyword{misc}
#> END exmplESet

#> BEGIN exmplMSet
##> \title{Example of Methylation Set}
#> \description{an ExpressionSet class contains minimum of exprs (matrix) of DNA methylation and phenoData (AnnotatedDataFrame of end point data).}
#> \usage{data(exmplMSet)}
#> \keyword{misc}
#> END exmplMSet

#> BEGIN exmplPat
#> \title{Example of Phenotype Pattern Definition Set}
#> \description{An exmple of phenotype pattern definition set with three columns: stat.coef, stat.func, and endpt.vars; It defines an association pattern for three phenotypes.}
#> \usage{data(exmplPat)}
#> \keyword{misc}
#> END exmplPat


##############################################################################
#> BEGIN datSel
#> \title{Subset of Genetic Data}
#> \description{Subset of genetic data for a gene according to probe names}
                                #> \arguments{
datSel<-function(dat,           #> \item{dat}{a data frame of genetic data with row as probes and column as subjects}
                 thisPrb,       #> \item{thisPrb}{a vector of probe names that are annotated to a gene}
                 prbcol=TRUE)   #> \item{prbcol}{indicator to return selected data: TRUE, columns are probes; FALSE, rows are probes}
                                #> }
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
#> \details{The function extract the genetic data for a gene according to probe names that are annotated the gene. 
#>  The function is internally called by \emph{CANN}.}
#> \value{a data frame}
#> END datSel

##############################################################################
#> BEGIN CANN
#> \title{Canonical Correlation of Two Sets of Genomic Data}
#> \description{Compute canonical correlation between two sets of genomic data.}
                              #> \arguments{
CANN<-function(geneSet,       #> \item{geneSet}{a gene set collection to annotate probes to gene}
               Edat,          #> \item{Edat}{data frame of the first form of genomic data, such as gene expression data with row being probes and column being subjects. The column names should match the row names \emph{phdat}}
               Mdat,          #> \item{Mdat}{data frame of the second form of genomic data, such as methylation data with row being probes and column being subjects. The column names should match the row names \emph{phdat}}
               EMlbl=c("Expr", "Methyl"),      #> \item{EMlbl}{lablel of the genomic data, default=c("Expr", "Methyl") for \emph{Edat} and \emph{Mdat}}
               phdat)         #> \item{phdat}{phenotype data with row being subjects and column being phenotype variables. The row names should match the column names of \emph{Edat} and \emph{Mdat}}
                              #> }
{
  res<-NULL
  gen<-NULL
  candtl<-NULL
  sbjnam<-dimnames(phdat)[[1]]
  for (i in 1:length(geneSet))
  { 
    thisgene<-geneSet[[i]]
    thisPrb<-geneIds(thisgene)
    if (is.null(Mdat)) thisMthyl2<-NULL
    if (!is.null(Mdat)) thisMthyl2<-datSel(Mdat, thisPrb)   
    if (is.null(Edat)) thisExpr2<-NULL
    if (!is.null(Edat))   thisExpr2<-datSel(Edat, thisPrb)
    
    
    if (is.null(thisMthyl2))
    {
      x=matrix(prcomp(thisExpr2, retx=TRUE, scale=TRUE)$x[,1], nrow=nrow(thisExpr2), ncol=1)
      y<-NULL
      this.res<-matrix(c(setName(thisgene), ncol(thisExpr2), 0, rep(NA, 3)), nrow=1, ncol=6) 
      this.gen<-c(setName(thisgene), as.vector(x), rep(NA, nrow(x)))
      this.candtl<-NULL    
    }
    
    if (is.null(thisExpr2))
    {
      y=matrix(prcomp(thisMthyl2, retx=TRUE, scale=TRUE)$x[,1], nrow=nrow(thisMthyl2), ncol=1)
      x<-NULL
      this.res<-matrix(c(setName(thisgene), 0, ncol(thisMthyl2), rep(NA, 3)), nrow=1, ncol=6) 
      this.gen<-c(setName(thisgene), rep(NA, nrow(y)), as.vector(y))
      this.candtl<-NULL 
    }
    
    if (!is.null(thisExpr2) & !is.null(thisMthyl2))
    {
      this.candtl<-NULL
      cc.res <-cancor(thisExpr2, thisMthyl2) 
      Mcnt<-thisMthyl2
      for(j in 1:ncol(Mcnt))  Mcnt[, j]<-Mcnt[, j] - cc.res$ycenter[j]
      Ecnt<-thisExpr2
      for(j in 1:ncol(Ecnt)) Ecnt[, j]<-Ecnt[, j] - cc.res$xcenter[j] 
      cr<-cc.res$cor[1]
      x=as.matrix(Ecnt[, 1:nrow(cc.res$xcoef)])%*%cc.res$xcoef[,1]*sign(cr)
      y=as.matrix(Mcnt[, 1:nrow(cc.res$ycoef)])%*%cc.res$ycoef[,1]
      this.candtl<-c(setName(thisgene), paste(names(cc.res$xcoef[,1]), collapse="|"),  paste(cc.res$xcoef[,1], collapse="|"), 
                     paste(names(cc.res$ycoef[,1]), collapse="|"),  paste(cc.res$ycoef[,1], collapse="|"))
      ccperm<-p.perm(thisExpr2, thisMthyl2, nboot = 999, rhostart = 1, type = "Wilks")$p
      ccasym<-p.asym(cc.res$cor, nrow(thisExpr2), ncol(thisExpr2), ncol(thisMthyl2), tstat = "Wilks")$p[1]
      
      this.res<-matrix(c(setName(thisgene), ncol(thisExpr2), ncol(thisMthyl2), cr, ccperm, ccasym), nrow=1, ncol=6)
      this.gen<-c(setName(thisgene), as.vector(x), as.vector(y))    
    }
    dimnames(this.res)[[2]]<-c('Gene', paste('n.', EMlbl, sep=""), 'CanonicalCR', 'WilksPermPval', 'WilksAsymPval')
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
  return(list(res=res, gen=gen, candtl=candtl))
}
#> \details{The function performs Canonical correlation between two forms genomic data for each gene (Edat and Mdat) defined by \emph{gann}. If a gene only has one form of genomic data, the first 
#>           prrcinpal component is used; If one form of data has numberof probesets exceeding the number of subjects, the first number of subjects probesets are used. The function return a list 
#>           of three components. See \emph{value} for details.  
#> }
#> \value{
#>       \item{$res}{canonical correlation result: a data frame with row for each each gene and six columns (Gene: gene names; n.EMlbl[1]: number of probes of first form genomic data;
#>           n.EMlbl[2]: number of probes of second form genomic data; CanonicalCR: Canonical correlation of first components; WilksPermPval: permuatation p value of Wilks' Lambda;
#>           WilksAsymPval: p value of F-approximations of Wilks' Lambda).}  
#>       \item{$gen}{the first component of canonical correlation: a data frame with row for each gene, first half of columns for first component of first form genomic 
#>           data and second half of columns for first component of second form genomic data.}
#>       \item{$candtl}{a data frame of loading (each row is for a gene, first column is gene names, 
#>           second column is the probeset ids of first form genomic data seperated by '|', third column is the load for each probeset in first form genomic data seperated by '|',  
#>           fourth column is the probeset ids of second form genomic data seperated by '|', fifth column is the load for each probeset in second form genomic data seperated by '|')}
#> }
#> \author{Xueyuan Cao \email{Xueyuan.cao@stjude.org}, Stanley Pounds \email{stanley.pounds@stjude.org}}
#> \references{Hotelling H. (1936). Relations between two sets of variables. Biometrika, 28, 321-327} 
#> \seealso{\code{\link{CCPROMISE}}}
#> \examples{
#>  ## load  exmplEdat exmplMdat
#>  data(exmplESet)
#>  data(exmplMSet)
#>  data(exmplGeneSet)
#>  ## Perform canonical correlation test
#> test1<- CANN(geneSet=exmplGeneSet, 
#>              Edat=exprs(exmplESet), 
#>              Mdat=exprs(exmplMSet), 
#>              EMlbl=c("Expr", "Methyl"), 
#>              phdat=pData(exmplESet))   
#> }
#> END CANN


##############################################################################
#> BEGIN PROMISE2
#> \title{PROMISE Analysis of Two Genomic Sets}
#> \description{PROMISE analysis of two genomic sets with multiple phenotypes.}
                                         #> \arguments{ 
PROMISE2<- function (exprSet,            #> \item{exprSet}{expression set of first genomic data}
                     exprSet2,           #> \item{exprSet2}{expression set of second genomic data}
                     geneSet = NULL,     #> \item{geneSet}{geneSet should be NULL.}
                     promise.pattern,    #> \item{promise.pattern}{PROMISE pattern}
                     strat.var = NULL,   #> \item{strat.var}{stratum variable}
                     nbperm=FALSE,       #> \item{nbperm}{indicator of fast permuation using negative binomial strategy, taking two valid values: FALSE or TRUE. The default is FALSE.}
                     max.ntail=100,      #> \item{max.ntail}{number of sucess if nbperm = T. Further permutation will not be performed for gene(s) or gene set(s) which max.ntail permutated statistics are greater or equal to the observed statistics, The default is 100.}
                     nperms = 10000,     #> \item{nperms}{number of permutation, default = 10,000}
                     seed = 13)          #> \item{seed}{random seed, default = 13}
                                         #> }
{
    ph.pattern <- promise.pattern
    endpt.data <- pData(phenoData(exprSet))
    array.data <- exprs(exprSet)
    array.data2 <- exprs(exprSet2)
    if (is.data.frame(array.data)) 
        array.data <- as.matrix(array.data)
    if (!is.numeric(array.data)) {
        print("ERROR: The array data is not numeric")
        return()
    }
    array.names <- dimnames(array.data)[[2]]
    endpt.names <- dimnames(endpt.data)[[1]]
    array.inc <- is.element(array.names, endpt.names)
    endpt.inc <- is.element(endpt.names, array.names)
    nsamp.inc <- sum(array.inc)
    message(paste("No. of Matching array ids in endpt data and array data:", 
        nsamp.inc))
    if (nsamp.inc <= 4) {
        message(paste("Array IDs in endpt data:"))
        message(endpt.names)
        message("Array IDs in array data:")
        message(array.names)
        stop("Error: Unable to match array identifiers in endpt data and array data.")
    }
    Y <- array.data[, array.inc]
    Y2 <- array.data2[, array.inc]
    X <- endpt.data[endpt.inc, ]
    array.names <- array.names[array.inc]
    endpt.names <- endpt.names[endpt.inc]
    array.ord <- order(array.names)
    endpt.ord <- order(endpt.names)
    Y <- Y[, array.ord]
    Y2 <- Y2[, array.ord]
    X <- X[endpt.ord, ]
    Y <- as.matrix(Y)
    Y2 <- as.matrix(Y2)
    
    array.names <- array.names[array.ord]
    endpt.names <- endpt.names[endpt.ord]
    if (any(array.names != endpt.names)) {
        name.mtx <- cbind(endpt.arrayid = endpt.names, array.arrayid = array.names)
        message("----------------- IDs ------------------")
        message(name.mtx)
        stop("Error: Unable to align IDs in array data and endpt data")
    }
    endpt.endptcol <- is.element(dimnames(X)[[2]], unique(unlist(strsplit(c(as.character(ph.pattern$endpt.vars), 
        strat.var), ","))))
    phtype <- X[, endpt.endptcol]
    names(phtype) <- dimnames(endpt.data)[[2]][endpt.endptcol]
    if (!is.numeric(phtype[, dimnames(phtype)[[2]] != strat.var]) & 
        !is.data.frame(phtype[, dimnames(phtype)[[2]] != strat.var])) {
        print("ERROR: The clinical data is not numeric or data frame")
        return()
    }
    if (is.vector(phtype)) 
        phtype <- matrix(phtype, length(phtype), 1)
    probes <- dimnames(array.data)[[1]]
    nprobes <- length(probes)
    if (!is.null(geneSet)) {
       print("ERROR: geneSet should be null")
        return()
    }
    
    if(is.null(geneSet)) {
        message(paste("Computing Stats for Observed Data:", date()))
        message(paste("Computing gene-specific statistics:", 
            date()))
        if (is.null(strat.var)) 
            strat <- NULL
        else strat <- phtype[, strat.var]
        gene.res1 <- PROMISE::promise.genestat(Y, phtype, ph.pattern, strat)
        gene.res2 <- PROMISE::promise.genestat(Y2, phtype, ph.pattern, strat)
        dimnames(gene.res1)[[2]]<-paste(dimnames(gene.res1)[[2]], "_1", sep="")
        dimnames(gene.res2)[[2]]<-paste(dimnames(gene.res2)[[2]], "_2", sep="")
        gene.resf<-cbind(gene.res1, gene.res2, 
                    PROMISE.stat=apply(cbind(gene.res1[, ncol(gene.res1)], gene.res2[, ncol(gene.res2)]), 1, function(prb){return(mean(prb, na.rm=T))}))
        
        m.nph <- dim(gene.resf)
        gene.pvals <- matrix(0, m.nph[1], m.nph[2])
        message(paste("Set seed: seed = ", seed))
        set.seed(seed)
        if(!nbperm){
           perms <- matrix(NA, nperms, nsamp.inc)
           if (is.null(strat.var)) 
               for (i in 1:nperms) perms[i, ] <- sample(nsamp.inc, replace = FALSE)
           else {
               strat <- endpt.data[, strat.var]
               ustrat <- sort(unique(strat))
               nstrat <- length(ustrat)
               ind <- length(strat)
               for (j in 1:nstrat) {
                   this.strat <- (strat == ustrat[j])
                   for (i in 1:nperms) perms[i, this.strat] <- sample((1:ind)[this.strat], replace = FALSE)
               }
           }
           message(paste("Computing Stats for Permuted Data: ", date()))
           for (i in 1:nperms) {
               perm.ind <- unlist(perms[i, ])
               perm.phtype <- phtype[perm.ind, ]
               if (is.null(strat.var)) 
                   perm.strat <- NULL
               else perm.strat <- perm.phtype[, strat.var]
               gene.temp1 <- PROMISE::promise.genestat(Y, perm.phtype, ph.pattern, perm.strat)
               gene.temp2 <- PROMISE::promise.genestat(Y2, perm.phtype, ph.pattern, perm.strat)
               
               dimnames(gene.temp1)[[2]]<-paste(dimnames(gene.temp1)[[2]], "_1", sep="")
               dimnames(gene.temp2)[[2]]<-paste(dimnames(gene.temp2)[[2]], "_2", sep="")
               gene.tempf<-cbind(gene.temp1, gene.temp2, 
                        PROMISE.stat=apply(cbind(gene.temp1[, ncol(gene.temp1)], gene.temp2[, ncol(gene.temp2)]), 1, function(prb){return(mean(prb, na.rm=T))}))
                                                   
               gene.pvals <- gene.pvals + (abs(gene.tempf) >= abs(gene.resf))
           }
           message(paste("Generating Result List: ", date()))
           gene.pvals <- gene.pvals/nperms
           generes.tab <- cbind.data.frame(probeid = probes, gene.resf, perm.p = gene.pvals)
           message(paste("Finished at: ", date()))
           res <- list(generes = generes.tab, setres = NULL)
       }
       if (nbperm){
           #set up tracking
           probe.keep <- rep(T, nprobes)
           probe.done <- rep(F, nprobes)
           probe.ntail <- matrix(0, m.nph[1], m.nph[2])
           probe.apt.ntail <- rep(0, m.nph[1])
           num.perms <- rep(NA, nprobes)
           i <- 0
           while ((i < nperms) && (any(probe.keep))) {
               i <- i + 1
               if (is.null(strat.var)) 
                  perm.ind <- sample(nsamp.inc, replace = FALSE)
               else {
                  perm.ind<-1:nsamp.inc 
                  strat <- X[, strat.var]
                  ustrat <- sort(unique(strat))
                  nstrat <- length(ustrat)
                  ind <- length(strat)
                  for (j in 1:nstrat) {
                    this.strat <- (strat == ustrat[j])
                    perm.ind[this.strat] <- sample((1:ind)[this.strat], replace = FALSE)
                  }
                }
                perm.phtype <- phtype[perm.ind, ]
                if (is.null(strat.var)) 
                  perm.strat <- NULL
                else perm.strat <- perm.phtype[, strat.var]
                gene.temp1 <- PROMISE::promise.genestat(Y[probe.keep, ], perm.phtype, ph.pattern, perm.strat)
                gene.temp2 <- PROMISE::promise.genestat(Y2[probe.keep, ], perm.phtype, ph.pattern, perm.strat)
                dimnames(gene.temp1)[[2]]<-paste(dimnames(gene.temp1)[[2]], "_1", sep="")
                dimnames(gene.temp2)[[2]]<-paste(dimnames(gene.temp2)[[2]], "_2", sep="")
                gene.tempf<-cbind(gene.temp1, gene.temp2, 
                        PROMISE.stat=apply(cbind(gene.temp1[, ncol(gene.temp1)], gene.temp2[, ncol(gene.temp2)]), 1, function(prb){return(mean(prb, na.rm=T))}))
 
                probe.ntail[probe.keep, ] <- probe.ntail[probe.keep, 
                  ] + (abs(gene.tempf) >= abs(gene.resf[probe.keep, 
                  ]))
                probe.apt.ntail[probe.keep] <- probe.ntail[probe.keep, 
                  m.nph[2]]
                probe.done[probe.keep] <- probe.apt.ntail[probe.keep] == 
                  max.ntail
                num.perms[probe.done & probe.keep] <- i
                gene.pvals[probe.done & probe.keep, ] <- probe.ntail[probe.done & 
                  probe.keep, ]/i
                probe.keep[probe.done] <- F
            }
            if (any(probe.keep)) {
                gene.pvals[probe.keep, ] <- probe.ntail[probe.keep, 
                  ]/nperms
                num.perms[probe.keep] <- nperms
            }
            dimnames(gene.pvals)[[2]] <- dimnames(gene.resf)[[2]]
            generes.tab <- cbind.data.frame(probeid = probes, 
                gene.resf, perm.p = gene.pvals, nperms = num.perms)
            message(paste("Finished at: ", date()))
            res <- list(generes = generes.tab, setres = NULL)
     }     
    }
    return(res)
}
#> \details{The function performs PROMISE analysis for two set genomic data with a prefined phenotypic pattern. It is intermediate function called by \emph{CCPROMISE} 
#>  to perform PROMISE analysis with canonical correlation}
#> \value{
#>       \item{$generes}{individual genes' test statistics and p-values for each individual endpoint and PROMISE analysis.}  
#>       \item{$setres}{Gene set level analysis is not implemented with value \emph{NULL}} 
#> }
#> \author{Xueyuan Cao \email{Xueyuan.cao@stjude.org}, Stanley Pounds \email{stanley.pounds@stjude.org}}
#> \seealso{\code{\link{CCPROMISE}}}
#> END PROMISE2

##############################################################################
#> BEGIN CCPROMISE
#> \title{PROMISE Analysis with Canonical Correlation for Two Forms of Genomic Data}
#> \description{PROMISE analysis of two genomic sets with multiple phenotypes under a predefined association pattern at gene level.}
                                                #> \arguments{ 
CCPROMISE<-function(geneSet,                    #> \item{geneSet}{a gene set collection to annotate probes to gene}
                  ESet,                         #> \item{ESet}{an ExpressionSet class contains minimum of exprs (expression matrix) of first form of genomic data such as gene expression and phenoData (AnnotatedDataFrame of end point data). Please refer to Biobase for details on how to create such an ExpressionSet expression set.}
                  MSet,                         #> \item{MSet}{an ExpressionSet class of second form of genomic data such as methylation levels, the subject id of MSet and ESet should be exactly same}
                  promise.pattern,              #> \item{promise.pattern}{PROMISE pattern}
                  strat.var=NULL,               #> \item{strat.var}{stratum variable}
                  prlbl=NULL,                   #> \item{prlbl}{labels}
                  EMlbl=c('Expr', 'Mthyl'),     #> \item{EMlbl}{lablel of the genomic data, default=c('Expr', 'Methyl') for ESet and MSet}
                  nbperm=FALSE,                 #> \item{nbperm}{indicator of fast permuation using negative binomial strategy, taking two valid values: FALSE or TRUE. The default is FALSE.}
                  max.ntail=100,                #> \item{max.ntail}{number of sucess if nbperm = T. Further permutation will not be performed for gene(s) or gene set(s) which max.ntail permutated statistics are greater or equal to the observed statistics, The default is 100.}
                  nperms=10000,                 #> \item{nperms}{number of permutation, default = 10,000}
                  seed=13)                      #> \item{seed}{initial seed of random number generator. The default is 13.}
                                                #> }
{
  Edat<-exprs(ESet)
  Mdat<-exprs(MSet)
  pheDat<-phenoData(ESet)
  #pheDatM<-phenoData(Mdat)
  phdat<-pData(pheDat)
  Edat<-Edat[, is.element(dimnames(Edat)[[2]], dimnames(Mdat)[[2]]) & is.element(dimnames(Edat)[[2]], dimnames(phdat)[[1]])]
  Mdat<-Mdat[, is.element(dimnames(Mdat)[[2]], dimnames(Edat)[[2]]) & is.element(dimnames(Mdat)[[2]], dimnames(phdat)[[1]])]
  phdat<-phdat[is.element(dimnames(phdat)[[1]], dimnames(Edat)[[2]]) & is.element(dimnames(phdat)[[1]], dimnames(Mdat)[[2]]), ]
  
  #Edat<-Edat[, order(dimnames(Edat)[[2]])]
  #Mdat<-Mdat[, order(dimnames(Mdat)[[2]])]
  #phdat<-phdat[order(dimnames(phdat)[[1]]), ]
  if(any(dimnames(Edat)[[2]]!=dimnames(Mdat)[[2]] | dimnames(Edat)[[2]]!=dimnames(phdat)[[1]] ))
  {
    message('Data Error!!!!!!!!!!!!!!!!!!!!!!!!!')
    break()
  }
  
  canres<-CANN(geneSet=geneSet, Edat=Edat, Mdat=Mdat, EMlbl=EMlbl, phdat=phdat)
  canndat<-canres$gen
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
  if (is.null(prlbl ))  prlbl<-c(promise.pattern[, 'endpt.vars'], paste('PR', nrow(promise.pattern), sep=""))
  if (!nbperm)
  dimnames(PRres)[[2]]<-c('Gene', paste(paste(EMlbl[1], prlbl, sep="_"), 'Stat', sep="."), 
                            paste(paste(EMlbl[2], prlbl, sep="_"), 'Stat', sep="."), paste(prlbl[length(prlbl)], 'Stat', sep="."),
                            paste(paste(EMlbl[1], prlbl, sep="_"), 'Pval', sep="."), 
                            paste(paste(EMlbl[2], prlbl, sep="_"), 'Pval', sep="."), paste(prlbl[length(prlbl)], 'Pval', sep="."))
  if (nbperm) 
  dimnames(PRres)[[2]]<-c('Gene', paste(paste(EMlbl[1], prlbl, sep="_"), 'Stat', sep="."), 
                            paste(paste(EMlbl[2], prlbl, sep="_"), 'Stat', sep="."), paste(prlbl[length(prlbl)], 'Stat', sep="."),
                            paste(paste(EMlbl[1], prlbl, sep="_"), 'Pval', sep="."), 
                            paste(paste(EMlbl[2], prlbl, sep="_"), 'Pval', sep="."), paste(prlbl[length(prlbl)], 'Pval', sep="."), 'nperm')                                                              
  return(list(PRres=PRres, CANres=canres$res, FSTgen=canres$gen, CANload=canres$candtl))
} 
#> \details{The function performs PROMISE analysis for two forms of genomic data in minimal expression set format with a prefined phenotypic pattern. It calls two external function 
#> \emph{CANN} and \emph{PROMISE2}}
#> \value{
#>       \item{$PRres}{PROMISE result for the first component of canonical correlation between two forms of geneomic data.
#>           individual genes' test statistics and p-values for each individual endpoint and PROMISE analysis}
#>       \item{$CANres}{result of canonical correlation analysis with six columns: Gene: Gene names; n.EMlbl[1]: number of probe set in the first form data; 
#>           n.EMlbl[2]: number of probe set in the second form data; CanonicalCR: Canonical correlation of first components; WilksPermPval: permuatation p value of Wilks' Lambda;
#>           WilksAsymPval: p value of F-approximations of Wilks' Lambda.}  
#>       \item{$FSTgen}{load of first component: a data frame of loading (each row is for a gene, first column is gene names, 
#>           second column is the probeset ids of first form genomic data seperated by '|', third column is the load for each probeset in first form genomic data seperated by '|',  
#>           fourth column is the probeset ids of second form genomic data seperated by '|', fifth column is the load for each probeset in second form genomic data seperated by '|')}
#>       \item{$CANload}{a data frame of loading (each row is for a gene, first column is gene names, 
#>           second column is the probeset ids of first form genomic data seperated by '|', third column is the load for each probeset in first form genomic data seperated by '|',  
#>           fourth column is the probeset ids of second form genomic data seperated by '|', fifth column is the load for each probeset in second form genomic data seperated by '|')}
#> }
#> \author{Xueyuan Cao \email{Xueyuan.cao@stjude.org}, Stanley Pounds \email{stanley.pounds@stjude.org}}
#> \seealso{\code{\link{CANN}} \code{\link{PROMISE2}}}
#> \examples{
#>  ## load data
#>  data(exmplESet)
#>  data(exmplMSet)
#>  data(exmplGeneSet)
#>  data(exmplPat)

#>  ## Perform canonical correlation test
#> test<- CCPROMISE(geneSet=exmplGeneSet, 
#>              ESet=exmplESet, 
#>              MSet=exmplMSet, 
#>              promise.pattern=exmplPat,
#>              strat.var=NULL,
#>              prlbl=NULL, 
#>              EMlbl=c("Expr", "Methyl"),
#>              nbperm=FALSE,
#>              max.ntail=10,
#>              nperms=100,
#>              seed=13)   
#> }
#> END CCPROMISE

##############################################################################
#> BEGIN PrbCor
#> \title{Probe Level Correlation of Two Sets of Genomic Data}
#> \description{Compute Spearman correlation of all probe combination between two sets of genomic data within a gene.}
                              #> \arguments{
PrbCor<-function(geneSet,     #> \item{geneSet}{a gene set collection to annotate probes to gene}
               Edat,          #> \item{Edat}{data frame of the first form of genomic data, such as gene expression data with row being probes and column being subjects. The column names should match the row names \emph{phdat}}
               Mdat,          #> \item{Mdat}{data frame of the second form of genomic data, such as methylation data with row being probes and column being subjects. The column names should match the row names \emph{phdat}}
               EMlbl=c("Expr", "Methyl"),      #> \item{EMlbl}{lablel of the genomic data, default=c("Expr", "Methyl") for \emph{Edat} and \emph{Mdat}}
               phdat,         #> \item{phdat}{phenotype data with row being subjects and column being phenotype variables. The row names should match the column names of \emph{Edat} and \emph{Mdat}}
               pcut=0.05)     #> \item{pcut}{p value cutoff to eliminate probe pairs that are not significantly correlated. Default is 0.05}
                              #> }
{
  res<-NULL
  gen<-NULL
  sbjnam<-dimnames(phdat)[[1]]
  for (i in 1:length(geneSet))
  { 
    thisgene<-geneSet[[i]]
    thisPrb<-geneIds(thisgene)
    if (is.null(Mdat)) thisMthyl2<-NULL
    if (!is.null(Mdat)) thisMthyl2<-datSel(Mdat, thisPrb, prbcol=FALSE)   
    if (is.null(Edat)) thisExpr2<-NULL
    if (!is.null(Edat))   thisExpr2<-datSel(Edat, thisPrb, prbcol=FALSE)   
    
    if (is.null(thisMthyl2))
    {
      this.res<-cbind(rep(setName(thisgene), nrow(thisExpr2)), rownames(thisExpr2), rep(NA, nrow(thisExpr2)), rep(NA, nrow(thisExpr2)), rep(NA, nrow(thisExpr2)))
      rownames(this.res)<-rownames(thisExpr2)
      colnames(this.res)<-c('Gene', EMlbl, 'Spearman.rstat', 'Spearman.p')
      this.gen<-cbind(rownames(thisExpr2), thisExpr2, matrix(NA, nrow=nrow(thisExpr2), ncol=ncol(thisExpr2)))
      dimnames(this.gen)[[2]]<-c("ProbPair", rep(sbjnam, 2))
    }
    
    if (is.null(thisExpr2))
    {
      this.res<-cbind(rep(setName(thisgene), nrow(thisMthyl2)), rep(NA, nrow(thisMthyl2)),rownames(thisMthyl2),  rep(NA, nrow(thisMthyl2)), rep(NA, nrow(thisMthyl2)))
      rownames(this.res)<-rownames(thisMthyl2)
      colnames(this.res)<-c('Gene', EMlbl, 'Spearman.rstat', 'Spearman.p')
      this.gen<-cbind(rownames(thisMthyl2), matrix(NA, nrow=nrow(thisMthyl2), ncol=ncol(thisMthyl2)), thisMthyl2)
      dimnames(this.gen)[[2]]<-c("ProbPair", rep(sbjnam, 2))
    }
    
    if (!is.null(thisExpr2) & !is.null(thisMthyl2))
    {
      this.res<-NULL
      this.gen<-NULL
      for (k in 1:nrow(thisExpr2))
      {
         spres<-t(apply(thisMthyl2, 1, function(prb, eprb){
                      tt<-cor.test(prb, as.numeric(eprb), method="spearman")
                      return(c(round(tt$estimate, 4), round(tt$p.value, 10)))
                      }, thisExpr2[k, ]))
         tt.res<-cbind(rep(setName(thisgene), nrow(thisMthyl2)), rep(rownames(thisExpr2)[k], nrow(thisMthyl2)), rownames(spres), spres)
         rownames(tt.res)<-paste(rownames(thisExpr2)[k],rownames(spres), sep='*')
         
         # sign of methylation signal is changed to reflect the correlation with expression. This change make the PROMISE analysis easier.
         msign<-apply(thisMthyl2, 2, function(prb, corsign){
                       return(prb*corsign)}, sign(spres[, 1]))
         if (is.vector(msign)) msign<-matrix(msign, nrow=1, ncol=length(msign))
         tt.gen<-cbind(matrix(rep(as.numeric(thisExpr2[k, ]), nrow(thisMthyl2)), nrow=nrow(thisMthyl2), ncol=ncol(thisExpr2), byrow=T), msign)
         rownames(tt.gen)<-paste(rownames(thisExpr2)[k],rownames(spres), sep='*')
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
#> \details{The function performs Spearman correlation for all probe pairs between two forms genomic data within each gene (Edat and Mdat) defined by \emph{gann}. 
#>           If a gene only has one form of genomic data, the other form is coded as NA.  
#>           The function return a list of two components. See \emph{value} for details.  
#> }
#> \value{
#>       \item{$res}{spearman correlation result: a data frame with row for each probe pair with correlation p value < pcut and five columns Gene: Gene: Gene names; EMlbl[1]: probe id in the first form data; 
#>           EMlbl[2]: probe id in the second form data; Spearman.rstat: Spearman r statistics; Spearman.p: Spearman p value.}
#>       \item{$gen}{Probe level data: a data frame with row for each probe pairs, first half of columns for first form genomic 
#>           data and second half of columns for second form genomic data with sign reflecting the correlation of the probe pair.}
#> }
#> \author{Xueyuan Cao \email{Xueyuan.cao@stjude.org}, Stanley Pounds \email{stanley.pounds@stjude.org}}
#> \seealso{\code{\link{CCPROMISE}}}
#> \examples{
#>  ## load exmplPhDat exmplEdat exmplMdat
#>  data(exmplESet)
#>  data(exmplMSet)
#>  data(exmplGeneSet)
#>  ## Perform canonical correlation test
#> test1<- PrbCor(geneSet=exmplGeneSet, 
#>              Edat=exprs(exmplESet), 
#>              Mdat=exprs(exmplMSet), 
#>              EMlbl=c("Expr", "Methyl"), 
#>              phdat=pData(exmplESet))   
#> }
#> END PrbCor


##############################################################################
#> BEGIN PrbPROMISE
#> \title{PROMISE Analysis with Two Forms of Genomic Data at Probe Level}
#> \description{PROMISE analysis of two genomic sets with multiple phenotypes under a predefined association pattern at probe level.}
                                                #> \arguments{ 
PrbPROMISE<-function(geneSet,                   #> \item{geneSet}{a gene set collection to annotate probes to gene}
                  ESet,                         #> \item{ESet}{an ExpressionSet class contains minimum of exprs (expression matrix) of first form of genomic data such as gene expression and phenoData (AnnotatedDataFrame of end point data). Please refer to Biobase for details on how to create such an ExpressionSet expression set.}
                  MSet,                         #> \item{MSet}{an ExpressionSet class of second form of genomic data such as methylation levels, the subject id of MSet and ESet should be exactly same}
                  promise.pattern,              #> \item{promise.pattern}{PROMISE pattern}
                  strat.var=NULL,               #> \item{strat.var}{stratum variable}
                  prlbl=NULL,                   #> \item{prlbl}{labels}
                  EMlbl=c('Expr', 'Mthyl'),     #> \item{EMlbl}{lablel of the genomic data, default=c('Expr', 'Methyl') for ESet and MSet}
                  pcut=0.05,                    #> \item{pcut}{p value cutoff to eliminate probe pairs that are not significantly correlated. Default is 0.05}
                  nbperm=FALSE,                 #> \item{nbperm}{indicator of fast permuation using negative binomial strategy, taking two valid values: FALSE or TRUE. The default is FALSE.}
                  max.ntail=100,                #> \item{max.ntail}{number of sucess if nbperm = T. Further permutation will not be performed for gene(s) or gene set(s) which max.ntail permutated statistics are greater or equal to the observed statistics, The default is 100.}
                  nperms=10000,                 #> \item{nperms}{number of permutation, default = 10,000}
                  seed=13)                      #> \item{seed}{initial seed of random number generator. The default is 13.}
                                                #> }
{
  options(stringsAsFactors=F)
  Edat<-exprs(ESet)
  Mdat<-exprs(MSet)
  pheDat<-phenoData(ESet)
  #pheDatM<-phenoData(Mdat)
  phdat<-pData(pheDat)
  Edat<-Edat[, is.element(dimnames(Edat)[[2]], dimnames(Mdat)[[2]]) & is.element(dimnames(Edat)[[2]], dimnames(phdat)[[1]])]
  Mdat<-Mdat[, is.element(dimnames(Mdat)[[2]], dimnames(Edat)[[2]]) & is.element(dimnames(Mdat)[[2]], dimnames(phdat)[[1]])]
  phdat<-phdat[is.element(dimnames(phdat)[[1]], dimnames(Edat)[[2]]) & is.element(dimnames(phdat)[[1]], dimnames(Mdat)[[2]]), ]
  
  #Edat<-Edat[, order(dimnames(Edat)[[2]])]
  #Mdat<-Mdat[, order(dimnames(Mdat)[[2]])]
  #phdat<-phdat[order(dimnames(phdat)[[1]]), ]
  if(any(dimnames(Edat)[[2]]!=dimnames(Mdat)[[2]] | dimnames(Edat)[[2]]!=dimnames(phdat)[[1]] ))
  {
    message('Data Error!!!!!!!!!!!!!!!!!!!!!!!!!')
    break()
  }
  
  CORres<-PrbCor(geneSet=geneSet, Edat=Edat, Mdat=Mdat, EMlbl=EMlbl, phdat=phdat, pcut=pcut)
  cordat<-CORres$gen
  cordatf<-apply(cordat[, -1], 2, as.numeric)
  dimnames(cordatf)[[1]]<-cordat[,1]
   
  expr<-as.matrix(cordatf[, 1:(ncol(cordatf)/2)])
  mtyl<-as.matrix(cordatf[, (ncol(cordatf)/2+1):ncol(cordatf)])
  colnames(mtyl)<-colnames(expr)
  
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
  if (is.null(prlbl ))  prlbl<-c(promise.pattern[, 'endpt.vars'], paste('PR', nrow(promise.pattern), sep=""))
  if (!nbperm)
  dimnames(PRres)[[2]]<-c('Gene', paste(paste(EMlbl[1], prlbl, sep="_"), 'Stat', sep="."), 
                            paste(paste(EMlbl[2], prlbl, sep="_"), 'Stat', sep="."), paste(prlbl[length(prlbl)], 'Stat', sep="."),
                            paste(paste(EMlbl[1], prlbl, sep="_"), 'Pval', sep="."), 
                            paste(paste(EMlbl[2], prlbl, sep="_"), 'Pval', sep="."), paste(prlbl[length(prlbl)], 'Pval', sep="."))
  if (nbperm) 
  dimnames(PRres)[[2]]<-c('Gene', paste(paste(EMlbl[1], prlbl, sep="_"), 'Stat', sep="."), 
                            paste(paste(EMlbl[2], prlbl, sep="_"), 'Stat', sep="."), paste(prlbl[length(prlbl)], 'Stat', sep="."),
                            paste(paste(EMlbl[1], prlbl, sep="_"), 'Pval', sep="."), 
                            paste(paste(EMlbl[2], prlbl, sep="_"), 'Pval', sep="."), paste(prlbl[length(prlbl)], 'Pval', sep="."), 'nperm')                                                              
  return(list(PRres=PRres, CORres=CORres$res))
} 
#> \details{The function performs PROMISE analysis for two forms of genomic data in minimal expression set format with a prefined phenotypic pattern. It calls two external function 
#> \emph{PrbCor} and \emph{PROMISE2}}
#> \value{
#>       \item{$PRres}{PROMISE result for the first component of canonical correlation between two forms of geneomic data.
#>           individual genes' test statistics and p-values for each individual endpoint and PROMISE analysis}
#>       \item{$CORres}{result of spearman correlation analysis of probe pairs within a gene with five columns: Gene: Gene names; EMlbl[1]: probe id in the first form data; 
#>           EMlbl[2]: probe id in the second form data; Spearman.rstat: Spearman r statistics; Spearman.p: Spearman p value.}  
#> }
#> \author{Xueyuan Cao \email{Xueyuan.cao@stjude.org}, Stanley Pounds \email{stanley.pounds@stjude.org}}
#> \seealso{\code{\link{PrbCor}} \code{\link{PROMISE2}}}
#> \examples{
#>  ## load data
#>  data(exmplESet)
#>  data(exmplMSet)
#>  data(exmplGeneSet)
#>  data(exmplPat)

#>  ## Perform probe level PROMISE analysis
#>test<-PrbPROMISE(geneSet=exmplGeneSet, 
#>              ESet=exmplESet, 
#>              MSet=exmplMSet, 
#>              promise.pattern=exmplPat,
#>              strat.var=NULL,
#>              prlbl=c('LC50', 'MRD22', 'EFS', 'PR3'), 
#>              EMlbl=c("Expr", "Methyl"),
#>              nbperm=TRUE,
#>              max.ntail=10,
#>              nperms=100,
#>              seed=13)
#> }
#> END PrbPROMISE


