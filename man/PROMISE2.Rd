\name{PROMISE2}
\alias{PROMISE2}
 \title{PROMISE Analysis of Two Genomic Sets}
 \description{PROMISE analysis of two genomic sets with multiple phenotypes.}
\usage{
PROMISE2 (exprSet, exprSet2, geneSet = NULL, promise.pattern, 
    strat.var = NULL, nbperm = FALSE, max.ntail = 100, nperms = 10000, 
    seed = 13) 
}
 \arguments{ 
 \item{exprSet}{expression set of first genomic data}
 \item{exprSet2}{expression set of second genomic data}
 \item{geneSet}{geneSet should be NULL.}
 \item{promise.pattern}{PROMISE pattern}
 \item{strat.var}{stratum variable}
 \item{nbperm}{indicator of fast permuation using negative binomial strategy, taking two valid values: FALSE or TRUE. The default is FALSE.}
 \item{max.ntail}{number of sucess if nbperm = T. Further permutation will not be performed for gene(s) or gene set(s) which max.ntail permutated statistics are greater or equal to the observed statistics, The default is 100.}
 \item{nperms}{number of permutation, default = 10,000}
 \item{seed}{random seed, default = 13}
 }
 \details{The function performs PROMISE analysis for two set genomic data with a prefined phenotypic pattern. It is intermediate function called by \emph{CCPROMISE} 
  to perform PROMISE analysis with canonical correlation}
 \value{The output of the function is a list of length 2. The 2 components are as following:
       \item{generes}{individual genes' test statistics and p-values for each individual endpoint and PROMISE analysis.}  
       \item{setres}{Gene set level analysis is not implemented with value \emph{NULL}} 
 }
 \author{Xueyuan Cao \email{Xueyuan.cao@stjude.org}, Stanley Pounds \email{stanley.pounds@stjude.org}}
 \seealso{\code{\link{CCPROMISE}}}
 \examples{
  ## load data
  data(exmplESet)
  data(exmplMSet)
  data(exmplGeneSet)
  data(exmplPat)
  ## Perform canonical correlation test
 test<- PROMISE2(exmplESet[1:10], 
                 exmplMSet[1:10], 
                 promise.pattern=exmplPat,
                 strat.var=NULL,
                 nbperm=FALSE,
                 max.ntail=10,
                 nperms=100,
                 seed=13)   
 }
