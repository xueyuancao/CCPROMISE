\name{PrbCor}
\alias{PrbCor}
 \title{Probe Level Correlation of Two Sets of Genomic Data}
 \description{Compute Spearman correlation of all probe combination between two sets of genomic data within a gene.}
\usage{
PrbCor (geneSet, Edat, Mdat, EMlbl = c("Expr", "Methyl"), phdat, 
    pcut = 0.05) 
}
 \arguments{
 \item{geneSet}{a gene set collection to annotate probes to gene}
 \item{Edat}{data frame of the first form of genomic data, such as gene expression data with row being probes and column being subjects. The column names should match the row names \emph{phdat}}
 \item{Mdat}{data frame of the second form of genomic data, such as methylation data with row being probes and column being subjects. The column names should match the row names \emph{phdat}}
 \item{EMlbl}{lablel of the genomic data, default=c("Expr", "Methyl") for \emph{Edat} and \emph{Mdat}}
 \item{phdat}{phenotype data with row being subjects and column being phenotype variables. The row names should match the column names of \emph{Edat} and \emph{Mdat}}
 \item{pcut}{p value cutoff to eliminate probe pairs that are not significantly correlated. Default is 0.05}
 }
 \details{The function performs Spearman correlation for all probe pairs between two forms genomic data within each gene (Edat and Mdat) defined by \emph{gann}. 
           If a gene only has one form of genomic data, the other form is coded as NA.  
           The function return a list of two components. See \emph{value} for details.  
 }
<<<<<<< HEAD
 \value{The output of the function is a list of length 2. The 2 components are as following:
       \item{res}{spearman correlation result: a data frame with row for each probe pair with correlation p value < pcut and five columns Gene: Gene: Gene names; EMlbl[1]: probe id in the first form data; 
           EMlbl[2]: probe id in the second form data; Spearman.rstat: Spearman r statistics; Spearman.p: Spearman p value.}
       \item{gen}{Probe level data: a data frame with row for each probe pairs, first half of columns for first form genomic 
=======
 \value{
       \item{$res}{spearman correlation result: a data frame with row for each probe pair with correlation p value < pcut and five columns Gene: Gene: Gene names; EMlbl[1]: probe id in the first form data; 
           EMlbl[2]: probe id in the second form data; Spearman.rstat: Spearman r statistics; Spearman.p: Spearman p value.}
       \item{$gen}{Probe level data: a data frame with row for each probe pairs, first half of columns for first form genomic 
>>>>>>> 5b4c40047f362e5ec59f764ac963f67c77446a2d
           data and second half of columns for second form genomic data with sign reflecting the correlation of the probe pair.}
 }
 \author{Xueyuan Cao \email{Xueyuan.cao@stjude.org}, Stanley Pounds \email{stanley.pounds@stjude.org}}
 \seealso{\code{\link{CCPROMISE}}}
 \examples{
  ## load exmplPhDat exmplEdat exmplMdat
  data(exmplESet)
  data(exmplMSet)
  data(exmplGeneSet)
  ## Perform canonical correlation test
 test1<- PrbCor(geneSet=exmplGeneSet, 
              Edat=exprs(exmplESet), 
              Mdat=exprs(exmplMSet), 
              EMlbl=c("Expr", "Methyl"), 
              phdat=pData(exmplESet))   
 }
