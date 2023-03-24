\name{cna}
\alias{cna}
\docType{data}
\title{
 Copy Number Alteration dataset
}
\description{
 A dataset containing the gene names, Entrez ID and CNA values for individuals. The variables are as follows:
}
\usage{data("cna")}
\format{
  A data frame with 24 observations on the following 1082 variables.
  \describe{
    \item{\code{Hugo_Symbol}}{Gene names}
    \item{\code{Entrez_Gene_Id}}{Unique identifiers for the gene}
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
data(cna)
}
\keyword{datasets}

\seealso{
\code{\link[=trios]{trios()}}
}

