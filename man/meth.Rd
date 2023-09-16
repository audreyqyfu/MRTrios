\name{meth}
\alias{meth}
\docType{data}
\title{
  DNA methylation dataset
}
\description{
 A dataset containing the methylation probes, gene names, and methylation values for individuals. It is a subset of a whole genome methylation dataset.
}
\usage{data("meth")}
\format{
  A data frame with 24 observations on the following 899 variables.  The first four variables are Row.names (i.e., methylation probe IDs), Gene_Symbol (the nearby genes where the methylation probe is located), Chromosome, and Genomic_Coordinate.  The remaining variables are IDs of individuals.
}
\examples{
 print(meth[1:5,1:5])
}
\keyword{datasets}

\seealso{
\code{\link[=trios]{trios}}
}

