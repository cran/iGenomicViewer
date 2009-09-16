\name{bandinfo-class}
\docType{class}

\alias{bandinfo}
\alias{BandInfo}

\title{bandinfo 'band information' object}

\description{
  The bandinfo object contains useful data.frames of starting and
  ending genomic locations for chromosomes, arms, broad.bands, and
  fine.bands. Each data.frame has five columns: the region (chrom, arm),
  label, lower, center, and upper. lower, center, and upper refer to
  genomic locations. label is the label that will be used for the region
  when plotting. The bandinfo object will also contain an offset
  object. This is a numeric vector of length equal to the number of
  chromsomes. It contains the numeric buffer that should be added to
  genomic locations given within chromosome to get the location with
  respect to the entire genome.
  
}

\section{mapobj object structure}{
  \describe{
  \item{offset}{offset is a numeric vector containing buffers to convert
    genomic locations with respect to chromosome to genomic locations with
    respect to the entire genome}
  \item{Chrom}{A data.frame with 5 columns: Chrom, Label, lower, center,
    and upper}
  \item{Arm}{A data.frame with 5 columns: Arm, Label, lower, center,
    and upper}
  \item{Broad.Band}{A data.frame with 5 columns: Broad.Band, Label,
    lower, center,and upper}
  \item{Fine.Band}{A data.frame with 5 columns: Fine.Band, Label, lower,
    center,and upper}
  }
}
\author{Lori A. Shepherd}

\seealso{\code{\link{makeBandInfo}}, \code{\link{Band.Info}}, \code{\link{cytoBand}}}

\keyword{classes}
