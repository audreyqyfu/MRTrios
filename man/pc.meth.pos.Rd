\name{pc.meth.pos}
\alias{pc.meth.pos}
\docType{data}
\title{
 Score matrix for ER + individuals in methylation dataset
}
\description{
 A data with the Prinicipal Component score matrix calculated from the genomewide gene expression datasets for each individual.
}
\usage{data("pc.meth.pos")}
\format{
  A data frame with 572 observations on the following 572 variables.
  \describe{
    \item{\code{PCx}}{PC score for each individual}
  }
}
\details{
 The individuals are matched between the gene expression and methylation datasets so their numbers are lower compared to the dimension of the respective datasets.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(pc.meth.pos)
}
\keyword{datasets}

\seealso{
\code{\link[=findPCs]{findPCs()}}
}
