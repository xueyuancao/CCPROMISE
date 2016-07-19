\name{datSel}
\alias{datSel}
 \title{Subset of Genetic Data}
 \description{Subset of genetic data for a gene according to probe names}
\usage{
datSel (dat, thisPrb, prbcol = TRUE) 
}
 \arguments{
 \item{dat}{a data frame of genetic data with row as probes and column as subjects}
 \item{thisPrb}{a vector of probe names that are annotated to a gene}
 \item{prbcol}{indicator to return selected data: TRUE, columns are probes; FALSE, rows are probes}
 }
 \details{The function extract the genetic data for a gene according to probe names that are annotated the gene. 
  The function is internally called by \emph{CANN}.}
 \value{a data frame}
 \examples{
  ## load  exmplMSet exmplGeneSet
  data(exmplMSet)
  data(exmplGeneSet)
  ## sel Methylation data corresponding to first gene set
  seldata<-datSel(exprs(exmplMSet), geneIds(exmplGeneSet[[1]])) 
 }
