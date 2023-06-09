% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findPCs.R
\name{findPCs}
\alias{findPCs}
\title{A function to calculate Principal Components (PCs) and significantly associated PCs with the data.}
\usage{
findPCs(data, startCol, GeneNameCol, type.ind, com.ind, type, bsize)
}
\arguments{
\item{data}{A data matrix; probes or genes in rows, individuals in columns}

\item{startCol}{The column where numeric values begin in the data}

\item{GeneNameCol}{The column where the gene names are located}

\item{type.ind}{A vector of positive or negative ER individuals ID}

\item{com.ind}{A vector of common individuals between the datasets you want to find PCs for}

\item{type}{A string of individual type ("ER pos" or "ER neg")}

\item{bsize}{The number of columns of data to use in each block of correlation calculations}
}
\value{
A list of length 4:

pca.type

\if{html}{\out{<div class="sourceCode">}}\preformatted{     A PC score matrix
}\if{html}{\out{</div>}}

data.with.conf

\if{html}{\out{<div class="sourceCode">}}\preformatted{     A list of 7 elements with the first element being a list of significantly associated PCs to each gene or probe
}\if{html}{\out{</div>}}

new.data.no.var

\if{html}{\out{<div class="sourceCode">}}\preformatted{     A updated data matrix used for analysis with individuals in rows and probes or genes in column
}\if{html}{\out{</div>}}
}
\description{
This function uses Principal Component Analysis which can be used to analyze the high dimensional data while retaining the information that is important in the data. The PC score matrix is the coefficients of the linear combination of the initial variables that form the Principle Components and we find the PCs that are highly associated with the data.
}
\examples{
#load the datasets
data(clinical.neg)
data(clinical.pos)

#Find common individuals between the methylation and gene expression dataset
com.ind = intersect(colnames(gene)[3:ncol(gene)], colnames(meth)[5:ncol(meth)])

#Use the function to get PC score matrix and significantly associated PCs.
pc.gene = findPCs(gene, 3, 1, clinical.pos[,1], com.ind, "Pos", 1)

#The PC matrix
pc.gene[[1]]

#List of significantly associated PCs
pc.gene[[2]]$sig.asso.covs

}
\seealso{
\code{\link[=prcomp]{prcomp()}} used to calculate the PC score matrix; \code{\link[=get.conf.matrix]{get.conf.matrix()}} used to calculate the significantly associated PCs
}
