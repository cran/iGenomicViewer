\name{writeExFiles}
\alias{writeExFiles}

\title{Writes Sample Files}

\description{
  This function writes example annotation, mapping, and band information files.
}

\usage{
writeExFiles(direct="./")
}

\arguments{
\item{direct}{directory path where files will be written}
}

\details{
  This function writes out sample annotation, mapping, and band
  information files.
  see individual help files for more details
  \code{CancerGenes}
  \code{cytoBand}
  \code{DiseaseGenes}
  \code{DNArepairgenes}
  \code{HB19Kv2.HG18}
}

\value{
  Five files are written to the given directory
}

\author{Lori A. Shepherd}

\seealso{\code{\link{annObj}},\code{\link{Band.Info}}, \code{\link{mapping.info}} }

\examples{

library("iGenomicViewer")

# writes to current working directory
writeExFiles()

}

\keyword{methods}
