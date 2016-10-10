\name{CANN}
\alias{CANN}
 \title{Canonical Correlation of Two Sets of Genomic Data}
 \description{Compute canonical correlation between two sets of genomic data.}
\usage{
CANN (geneSet, Edat, Mdat, EMlbl = c("Expr", "Methyl"), phdat) 
}
 \arguments{
 \item{geneSet}{a gene set collection to annotate probes to gene}
 \item{Edat}{data frame of the first form of genomic data, such as gene expression data with row being probes and column being subjects. The column names should match the row names \emph{phdat}}
 \item{Mdat}{data frame of the second form of genomic data, such as methylation data with row being probes and column being subjects. The column names should match the row names \emph{phdat}}
 \item{EMlbl}{lablel of the genomic data, default=c("Expr", "Methyl") for \emph{Edat} and \emph{Mdat}}
 \item{phdat}{phenotype data with row being subjects and column being phenotype variables. The row names should match the column names of \emph{Edat} and \emph{Mdat}}
 }
 \details{The function performs Canonical correlation between two forms genomic data for each gene (Edat and Mdat) defined by \emph{gann}. If a gene only has one form of genomic data, the first 
           prrcinpal component is used; If one form of data has numberof probesets exceeding the number of subjects, the first number of subjects probesets are used. The function return a list 
           of three components. See \emph{value} for details.  
 }
<<<<<<< HEAD
 \value{ The output of the function is a list of length 3 with thee components:
       \item{CCres}{canonical correlation result: a data frame with row for each each gene and six columns (Gene: gene names; n.EMlbl[1]: number of probes of first form genomic data;
           n.EMlbl[2]: number of probes of second form genomic data; CanonicalCR: Canonical correlation of first components; WilksPermPval: permuatation p value of Wilks' Lambda;
           WilksAsymPval: p value of F-approximations of Wilks' Lambda).}  
       \item{FSTccscore}{the first component of canonical correlation: a data frame with row for each gene, first half of columns for first component of first form genomic 
           data and second half of columns for first component of second form genomic data.}
       \item{CCload}{a data frame of loading (each row is for a gene, first column is gene names, 
=======
 \value{
       \item{$res}{canonical correlation result: a data frame with row for each each gene and six columns (Gene: gene names; n.EMlbl[1]: number of probes of first form genomic data;
           n.EMlbl[2]: number of probes of second form genomic data; CanonicalCR: Canonical correlation of first components; WilksPermPval: permuatation p value of Wilks' Lambda;
           WilksAsymPval: p value of F-approximations of Wilks' Lambda).}  
       \item{$gen}{the first component of canonical correlation: a data frame with row for each gene, first half of columns for first component of first form genomic 
           data and second half of columns for first component of second form genomic data.}
       \item{$candtl}{a data frame of loading (each row is for a gene, first column is gene names, 
>>>>>>> 5b4c40047f362e5ec59f764ac963f67c77446a2d
           second column is the probeset ids of first form genomic data seperated by '|', third column is the load for each probeset in first form genomic data seperated by '|',  
           fourth column is the probeset ids of second form genomic data seperated by '|', fifth column is the load for each probeset in second form genomic data seperated by '|')}
 }
 \author{Xueyuan Cao \email{Xueyuan.cao@stjude.org}, Stanley Pounds \email{stanley.pounds@stjude.org}}
 \references{Hotelling H. (1936). Relations between two sets of variables. Biometrika, 28, 321-327} 
 \seealso{\code{\link{CCPROMISE}}}
 \examples{
  ## load  exmplEdat exmplMdat
  data(exmplESet)
  data(exmplMSet)
  data(exmplGeneSet)
  ## Perform canonical correlation test
 test1<- CANN(geneSet=exmplGeneSet, 
              Edat=exprs(exmplESet), 
              Mdat=exprs(exmplMSet), 
              EMlbl=c("Expr", "Methyl"), 
              phdat=pData(exmplESet))   
 }
