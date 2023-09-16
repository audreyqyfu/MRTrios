\name{clinical.pos}
\alias{clinical.pos}
\docType{data}
\title{
 Clinical dataset for ER - individuals
}
\description{
 A dataset containing the ID and demographic information of the individuals indentified with Estrogen receptor positive (ER +) breast cancer.
}
\usage{data("clinical.pos")}
\format{
  A data frame with 650 observations on the following 3 variables.
  \describe{
    \item{\code{ID}}{Individual ID}
    \item{\code{Age}}{Age of the individuals}
    \item{\code{Race}}{Race of the individuals}
  }
}
\examples{
data(clinical.pos)
}
\keyword{datasets}

\seealso{
\code{\link[=findPCs]{findPCs()}}
}
