\name{mapobj-class}
\docType{class}

\alias{mapobj}
\alias{mapObj}

\title{mapobj 'mapping object' class}

\description{
The mapping object contains all mapping information which includes but
is not limited to spotIDs, chromosome locations, and genomic
locations. The mapping object will also have a band.info object
associated with it in order to map spot.IDs to chromosome, arm, broad
and fine bands.  
}

\section{mapobj object structure}{
  \describe{
  \item{band.info}{a band.info object. see \code{makeBandInfo}}
  \item{mapping.info}{A data.frame of mapping information. The
    information included depends on the mapping file used and the user
    settins when creating the object with the mappingObj function. The
    data.frame must have Spot.ID, Chrom, as well as one of the options
    for genomic locations. see \code{mappingObj} for more details.}
  \item{links}{a data frame of web address information for the
    mapping. It should have the same row dimensions as the mapping.info
    data.frame. }
  \item{images}{a data frame of paths to images relating to the
    mapping. It should have the same row dimensions as the mapping.info
    data.frame. }
  }
  
}

\author{Lori A. Shepherd}

\seealso{\code{\link{mappingObj}}}

\keyword{classes}
