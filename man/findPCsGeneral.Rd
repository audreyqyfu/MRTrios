% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findPCsGeneral.R
\name{findPCsGeneral}
\alias{findPCsGeneral}
\title{A general function to derive Principal Components (PCs)
and identify significantly associated PCs with the input variables.}
\usage{
findPCsGeneral(data, startCol = 3, GeneNameCol = 1, com.ind, bsize)
}
\arguments{
\item{data}{A data matrix; probes or genes in rows, individuals in columns}

\item{startCol}{The column where numeric values begin in the data}

\item{GeneNameCol}{The column where the gene names are located}

\item{com.ind}{A vector of common individuals between the datasets you want to find PCs for}

\item{bsize}{The number of columns of data to use in each block of correlation calculations}
}
\value{
A list of length 3:

pca.type

\if{html}{\out{<div class="sourceCode">}}\preformatted{     A PC score matrix.
}\if{html}{\out{</div>}}

data.with.conf

\if{html}{\out{<div class="sourceCode">}}\preformatted{     A list of 6 elements:
             sig.asso.covs (list of associated PCs for each gene), 
             pvalues (data matrix), 
             qvalues (list of one data frame), 
             cors (data matrix), 
             sig (list of one data frame containing TRUE/FALSE),
             adj.p (data matrix).
}\if{html}{\out{</div>}}

new.data.no.var

\if{html}{\out{<div class="sourceCode">}}\preformatted{     A updated data matrix with individuals in rows and probes/genes in columns.
}\if{html}{\out{</div>}}
}
\description{
This function performs principal component (PC) analysis on the input data
and identifies scores that are significantly associated with the variables in the input.
}
\examples{
#load the datasets
data(gene)
data(meth)

#Find common individuals between the methylation and gene expression dataset
com.ind = intersect(colnames(gene)[3:ncol(gene)], colnames(meth)[5:ncol(meth)])

#Use the function to get PC score matrix and significantly associated PCs.
pc.gene = findPCsGeneral(as.data.frame(gene), 3, 1, com.ind, 1)

#The PC matrix
dim (pc.gene[[1]])
pc.gene[[1]]

#List of significantly associated PCs
pc.gene[[2]]$sig.asso.covs

}
\seealso{
\link{prcomp} used to calculate the PC score matrix;
\link{get.conf.matrix} used to identify PCs that are significantly associated with the input;
\link{findPCs} used to derive and identify PCs for the breast cancer cohort with two subtypes (ER+ and ER-).
}
