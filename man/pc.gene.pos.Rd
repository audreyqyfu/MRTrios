\name{pc.gene.pos}
\alias{pc.gene.pos}
\docType{data}
\title{
Principal component (PC) score matrix for the gene expression data in ER+ breast cancer patients
}
\description{
 The PC score matrix is derived from the genomewide gene expression datasets for each individual.  The individuals are matched between the gene expression and methylation datasets so their numbers are lower compared to the dimension of the respective datasets.
}
\usage{data(pc.gene.pos)}
\format{
  A data frame with 572 scores for 572 PCs.
}
\examples{
data(pc.gene.pos)
}

\keyword{datasets}

\seealso{
\code{\link[=findPCs]{findPCs}}
\code{\link[=findPCsGeneral]{findPCsGeneral}}
\code{\link{pc.gene.pos}}
}
