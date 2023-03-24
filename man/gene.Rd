\name{gene}
\alias{gene}
\docType{data}
\title{
 Gene Expression dataset
}
\description{
 A dataset containing the gene names, Entrez ID and expression values for individuals. The variables are as follows:
}
\usage{data("gene")}
\format{
  A data frame with 22 observations on the following 1102 variables.
  \describe{
    \item{\code{Hugo_Symbol}}{Gene names}
    \item{\code{Entrez_Gene_Id}}{Unique identifiers for the gene}
    \item{\samp{TCGA-3C-AAAU-01}}{Individual ID}
  }
}
\details{
The data provided in the package is a sample from the original data
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(gene)
}
\keyword{datasets}

\seealso{
\code{\link[=trios]{trios()}}
}

