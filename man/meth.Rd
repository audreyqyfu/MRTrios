\name{meth}
\alias{meth}
\docType{data}
\title{
  DNA methylation dataset
}
\description{
 A dataset containing the methylation probes, gene names, and methylation values for individuals. The variables are as follows:
}
\usage{data("meth")}
\format{
  A data frame with 24 observations on the following 899 variables.
  \describe{
    \item{\code{Row.names}}{Methylation Probe IDs}
    \item{\code{Gene_Symbol}}{Gene names consisting of the methylation probe}
    \item{\code{Chromosome}}{The number of chromosome (any number between 1 and 22 or X or Y)}
    \item{\code{Genomic_Coordinate}}{The position of the probe in the genome}
    \item{\samp{TCGA-xx-xxxx-01}}{Individual ID}
  }
}
\details{
The data provided in the package is a sample from the original data.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
 print(meth[1:5,1:5])
}
\keyword{datasets}

\seealso{
\code{\link[=trios]{trios()}}
}

