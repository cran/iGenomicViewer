\name{makeTrack}
\alias{makeTrack}
\alias{trackRegions}
\alias{trackRegion}

\title{Creates a trackRegion object}

\description{
  This function creates a track region object. This object is
  utilized in makeGGV to store known genomic regions of interest.
}

\usage{

makeTrack(Broad.Band=NA,
          Fine.Band=NA,
          genomicLoc=NA,
          geneName=NA)

}

\arguments{
\item{Broad.Band}{character vector of broad band regions. These regions
  should match to an associated mapObj\$band.info\$Broad.Band}
\item{Fine.Band}{character vector of fine band regions. These regions
  should match to an associated mapObj\$band.info\$Fine.Band}
\item{genomicLoc}{numeric vector or matrix of genomic locations. see
  details}
\item{geneName}{character vector of geneNames, or name in label column
  of an associated annObj's annotation data.frame.}
  
}

\details{
  The makeTrack function makes a trackRegion object to store known
  regions of interest. There are four ways to speicify a region: broad
  band, fine band, genomic start and stop locations, and geneName.
  
  Broad.Band and Fine.Band are character vectors containing names of
  broad band or fine band regions. They should match an entry in the
  broad.band or fine.band column in a mapObj\$band.info data frames
  object.
  The geneName character vector is associates with an annObj. This
  character vector should contain entries that match the labels in the
  columns of an annObj's annotation  data.frame. The list may contain
  entries from any and all objects of the annObj.
  The genomicLoc entry may either be a numeric vector or matrix. If it
  is a vector it assumes the entries are listed in order of start loc,
  stop loc, start loc, stop loc, for given regions. If it is a matrix,
  it assumes the first column is start locations and the second is stop
  locations.


  Any of the four entries may be left as NA. 
  
  A trackRegion object is returned. This object is a list of known
  regions, with entries for Broad.Band, Fine.Band, genomicLoc, and
  geneName.
  
}

\value{trackRegion object, list containing known regions of interest,
  see details}

\note{see vignette for more details}

\author{Lori A. Shepherd}

\seealso{\code{\link{updateGGV}}, \code{\link{makeGGV}}}

\examples{

library("iGenomicViewer")

trackRegion = makeTrack(Fine.Band = c("8p11.22","18p11.21"),
                        genomicLoc = NA,
                        geneName = "FANCE")

}

\keyword{methods}
