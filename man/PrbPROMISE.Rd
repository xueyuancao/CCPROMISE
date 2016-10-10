\name{PrbPROMISE}
\alias{PrbPROMISE}
 \title{PROMISE Analysis with Two Forms of Genomic Data at Probe Level}
 \description{PROMISE analysis of two genomic sets with multiple phenotypes under a predefined association pattern at probe level.}
\usage{
PrbPROMISE (geneSet, ESet, MSet, promise.pattern, strat.var = NULL, 
    prlbl = NULL, EMlbl = c("Expr", "Mthyl"), pcut = 0.05, nbperm = FALSE, 
    max.ntail = 100, nperms = 10000, seed = 13) 
}
 \arguments{ 
 \item{geneSet}{a gene set collection to annotate probes to gene}
 \item{ESet}{an ExpressionSet class contains minimum of exprs (expression matrix) of first form of genomic data such as gene expression and phenoData (AnnotatedDataFrame of end point data). Please refer to Biobase for details on how to create such an ExpressionSet expression set.}
 \item{MSet}{an ExpressionSet class of second form of genomic data such as methylation levels, the subject id of MSet and ESet should be exactly same}
 \item{promise.pattern}{PROMISE pattern}
 \item{strat.var}{stratum variable}
 \item{prlbl}{labels}
 \item{EMlbl}{lablel of the genomic data, default=c('Expr', 'Methyl') for ESet and MSet}
 \item{pcut}{p value cutoff to eliminate probe pairs that are not significantly correlated. Default is 0.05}
 \item{nbperm}{indicator of fast permuation using negative binomial strategy, taking two valid values: FALSE or TRUE. The default is FALSE.}
 \item{max.ntail}{number of sucess if nbperm = T. Further permutation will not be performed for gene(s) or gene set(s) which max.ntail permutated statistics are greater or equal to the observed statistics, The default is 100.}
 \item{nperms}{number of permutation, default = 10,000}
 \item{seed}{initial seed of random number generator. The default is 13.}
 }
 \details{The function performs PROMISE analysis for two forms of genomic data in minimal expression set format with a prefined phenotypic pattern. It calls two external function 
 \emph{PrbCor} and \emph{PROMISE2}}
<<<<<<< HEAD
 \value{The output of the function is a list of length 2. The 2 components are as following:
       \item{PRres}{PROMISE result for the first component of canonical correlation between two forms of geneomic data.
           individual genes' test statistics and p-values for each individual endpoint and PROMISE analysis}
       \item{CORres}{result of spearman correlation analysis of probe pairs within a gene with five columns: Gene: Gene names; EMlbl[1]: probe id in the first form data; 
=======
 \value{
       \item{$PRres}{PROMISE result for the first component of canonical correlation between two forms of geneomic data.
           individual genes' test statistics and p-values for each individual endpoint and PROMISE analysis}
       \item{$CORres}{result of spearman correlation analysis of probe pairs within a gene with five columns: Gene: Gene names; EMlbl[1]: probe id in the first form data; 
>>>>>>> 5b4c40047f362e5ec59f764ac963f67c77446a2d
           EMlbl[2]: probe id in the second form data; Spearman.rstat: Spearman r statistics; Spearman.p: Spearman p value.}  
 }
 \author{Xueyuan Cao \email{Xueyuan.cao@stjude.org}, Stanley Pounds \email{stanley.pounds@stjude.org}}
 \seealso{\code{\link{PrbCor}} \code{\link{PROMISE2}}}
 \examples{
  ## load data
  data(exmplESet)
  data(exmplMSet)
  data(exmplGeneSet)
  data(exmplPat)
  ## Perform probe level PROMISE analysis
test<-PrbPROMISE(geneSet=exmplGeneSet, 
              ESet=exmplESet, 
              MSet=exmplMSet, 
              promise.pattern=exmplPat,
              strat.var=NULL,
              prlbl=c('LC50', 'MRD22', 'EFS', 'PR3'), 
              EMlbl=c("Expr", "Methyl"),
              nbperm=TRUE,
              max.ntail=10,
              nperms=100,
              seed=13)
 }
