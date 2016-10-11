\name{CCPROMISE-package}
\alias{CCPROMISE-package}
 \Rdversion{3.3.1}
 \docType{package}
 \title{
   PRojection Onto the Most Interesting Statistical Evidence with Canonical Correlation
 }
 \description{
   a tool to indentify genes that are correlated between two set of genomic variables and are associated with a predefined  
   pattern of associations with multiple endpoint variables.
 }
 \details{
   \tabular{ll}{
     Package: \tab CCPROMISE\cr
     Type: \tab Package\cr
     Version: \tab 0.99.2\cr
     Date: \tab 2016-10-10\cr
     License: \tab GPL (>=2)\cr
     LazyLoad: \tab yes\cr
   }
   The CCPROMISE (Canonical correlation  with PROMISE analysis) is performed by calling function CCPROMISE.
   The two forms of genomic data such as gene expression and methylation are passed through minimal ExpressionSet;
   the gene annotation (defining relationship between a gene and the two forms of genomic data), phenotypic data and definition of 
   R routines for calculating association statistics with individual endpoint variable are same as in \emph{PROMISE} package. 
   Please refer to \emph{PROMISE} package for writing user defined routines.
 }
 \author{Xueyuan Cao \email{Xueyuan.cao@stjude.org}, Stanley Pounds \email{stanley.pounds@stjude.org}

      Maintainer: Xueyuan Cao \email{xueyuan.cao@stjude.org}
 }
 \references{
  Cao X, Crews KR, Downing J, Lamba J and Pounds XB (2016) 
  CC-PROMISE effectively integrates two forms of molecular data with multiple biologically related endpoint. 
  BMC Bioinformatics 17(Suppl 13):382
             
  Hotelling H. (1936). Relations between two sets of variables. Biometrika, 28, 321-327
  
  Pounds S, Cheng C, Cao X, Crews KR, Plunkett W, Gandhi V, Rubnitz J, Ribeiro RC, Downing JR, and Lamba J (2009) 
  PROMISE: a tool to identify genomic features with a specific biologically interesting pattern of associations 
  with multiple endpoint variables. Bioinformatics 25: 2013-2019
  
  Wilks, S. S. (1935) On the independence of k sets of normally distributed statistical variables. Econometrica, 3 309-326.
 }
 \keyword{ package }
 \examples{
 ## load data
  data(exmplESet)
  data(exmplMSet)
  data(exmplGeneSet)
  data(exmplPat)
  ## Perform CCPROMISE test
 test<- CCPROMISE(geneSet=exmplGeneSet, 
              ESet=exmplESet, 
              MSet=exmplMSet, 
              promise.pattern=exmplPat,
              strat.var=NULL,
              prlbl=NULL, 
              EMlbl=c("Expr", "Methyl"),
              nbperm=TRUE,
              max.ntail=10,
              nperms=100,
              seed=13)    
 }
