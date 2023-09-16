\name{biomart}
\alias{biomart}
\docType{data}
\title{
 Biomart data from Ensembl
}
\description{
 The biomart data extracted from Ensembl that contained information of the genes in human.
}
\usage{data("biomart")}
\format{
  A data frame with 5 rows and 6 columns.
  \describe{
    \item{\code{Gene.stable.ID}}{The identifiers that refer to the same genomic features or species}
    \item{\code{Gene.stable.ID.version}}{The increments when the set of transcripts linked to a gene changes}
    \item{\code{Gene.start..bp.}}{The location where the gene begins}
    \item{\code{Gene.end..bp.}}{The location where the gene ends}
    \item{\code{Gene.name}}{The gene name}
    \item{\code{Gene...GC.content}}{The content of Guanine-Cytosin for the gene}
  }
}
\examples{
data(biomart)
}

\keyword{datasets}

\seealso{
\code{\link[=extractHumanMethProbeInfo]{extractHumanMethProbeInfo}}
}
