\name{anninfo-class}

\docType{class}

\alias{anninfo}

\title{anninfo 'annotation information' class}

\description{
The anninfo object contains all information for a given annoation
set. This includes but is not limited to names, chromosome locations,
and genomic locations.
}

\section{anninfo object structure}{

 \describe{
	
  \item{annotation}{a data frame of information for the annotation
    set. The information included depends on the annotation file used
    and the user settings when creating the object with the
    makeAnnotation function. The data.frame must have Label and Chrom,
    as well as one of the options for genomic locations. see
    \code{makeAnnotation} for more details. }
  \item{links}{a data frame of web address information for the
    annotation set. It should have the same row dimension as the
    annotation data.frame.}
  \item{images}{a data frame of paths to images for the annotation
    set. It should have the same row dimensions as the annotation
    data.frame}
  }
}

\author{Lori A. Shepherd}

\seealso{\code{\link{annotation}}}

\keyword{classes}
