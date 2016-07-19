\name{CCPROMISE}
\alias{CCPROMISE}
 \title{PROMISE Analysis with Canonical Correlation for Two Forms of Genomic Data}
 \description{PROMISE analysis of two genomic sets with multiple phenotypes under a predefined association pattern at gene level.}
\usage{
CCPROMISE (geneSet, ESet, MSet, promise.pattern, strat.var = NULL, 
    prlbl = NULL, EMlbl = c("Expr", "Mthyl"), nbperm = FALSE, 
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
 \item{nbperm}{indicator of fast permuation using negative binomial strategy, taking two valid values: FALSE or TRUE. The default is FALSE.}
 \item{max.ntail}{number of sucess if nbperm = T. Further permutation will not be performed for gene(s) or gene set(s) which max.ntail permutated statistics are greater or equal to the observed statistics, The default is 100.}
 \item{nperms}{number of permutation, default = 10,000}
 \item{seed}{initial seed of random number generator. The default is 13.}
 }
 \details{The function performs PROMISE analysis for two forms of genomic data in minimal expression set format with a prefined phenotypic pattern. It calls two external function 
 \emph{CANN} and \emph{PROMISE2}}
 \value{
       \item{$PRres}{PROMISE result for the first component of canonical correlation between two forms of geneomic data.
           individual genes' test statistics and p-values for each individual endpoint and PROMISE analysis}
       \item{$CANres}{result of canonical correlation analysis with six columns: Gene: Gene names; n.EMlbl[1]: number of probe set in the first form data; 
           n.EMlbl[2]: number of probe set in the second form data; CanonicalCR: Canonical correlation of first components; WilksPermPval: permuatation p value of Wilks' Lambda;
           WilksAsymPval: p value of F-approximations of Wilks' Lambda.}  
       \item{$FSTgen}{load of first component: a data frame of loading (each row is for a gene, first column is gene names, 
           second column is the probeset ids of first form genomic data seperated by '|', third column is the load for each probeset in first form genomic data seperated by '|',  
           fourth column is the probeset ids of second form genomic data seperated by '|', fifth column is the load for each probeset in second form genomic data seperated by '|')}
       \item{$CANload}{a data frame of loading (each row is for a gene, first column is gene names, 
           second column is the probeset ids of first form genomic data seperated by '|', third column is the load for each probeset in first form genomic data seperated by '|',  
           fourth column is the probeset ids of second form genomic data seperated by '|', fifth column is the load for each probeset in second form genomic data seperated by '|')}
 }
 \author{Xueyuan Cao \email{Xueyuan.cao@stjude.org}, Stanley Pounds \email{stanley.pounds@stjude.org}}
 \seealso{\code{\link{CANN}} \code{\link{PROMISE2}}}
 \examples{
  ## load data
  data(exmplESet)
  data(exmplMSet)
  data(exmplGeneSet)
  data(exmplPat)
  ## Perform canonical correlation test
 test<- CCPROMISE(geneSet=exmplGeneSet, 
              ESet=exmplESet, 
              MSet=exmplMSet, 
              promise.pattern=exmplPat,
              strat.var=NULL,
              prlbl=NULL, 
              EMlbl=c("Expr", "Methyl"),
              nbperm=FALSE,
              max.ntail=10,
              nperms=100,
              seed=13)   
 }
