\name{TIplot-class}
\docType{class}
\alias{TIplot}

\title{TIplot 'tiled image' class}

\description{
  The tiled image object contains all necessary information to plot a
  tiled image heatmap. This includes but is not limited to matrix of
  values, axis ranges, axis labels, size, and mapping information. 
}

\section{TIplot structure}{
  \describe{
  \item{tractBound}{The tractBound object contains all information for
    breaking up the yaxis spot list into different tracks. For detailed
    information please see vignette}
  \item{vls}{The vls list contains information on values for
    plotting. There is a matrix or start and stop genomic locations
    for the y axix, Y. There is also the matrix of values, Z. The object
    H is a numeric indicating the number of tracks to break each sample
    into. Xcoords is a numeric list indicating each tracks x
    value. Finally nsmp and nBAC indicate the number of samples and the
    number of spots used. For more detail see vignette.}
  \item{lims}{The lims list contains vectors of the range for the xaxis,
    yaxis and Z matrix [xlim, ylim, zlim]}
  \item{labels}{The labels list contains values for individual x labels
    and ylabels, as well as labels for title, xaxis and yaxis. [ttl,
    xlab, ylab]}
  \item{Zcol}{Zcol is a vector of colors for the Z matrix scale}
  \item{cex}{The cex, or size, list consists of the x axis and y axis
    values for point size, given as xcex and ycex respectively}
  \item{map}{The map list consists of the spot ID index [bacDX] and the
    associated mapping object [mapObj]. Please see mapobj for details on
    this structure}
  }
}

\note{Users are encouraged to view the vignette for more detail}

\author{Lori A. Shepherd, Daniel P. Gaile}

\seealso{\code{\link{initTiled}}}

\keyword{classes}
